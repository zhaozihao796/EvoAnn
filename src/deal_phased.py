import logging
import os
from typing import Dict, Tuple, Optional, List
from collections import defaultdict
from . annotate_utils import fasta_to_dict, parse_attributes,rev_comp

try:
    from cyvcf2 import VCF
except ImportError:
    raise ImportError("Please install cyvcf2 first: pip install cyvcf2")

logger = logging.getLogger(__name__)

class PhasedExtractor:
    def __init__(self, ref_fasta_file: str, vcf_file: str, gff_file: str, output_dir: str = "phased_output",
                 cds_only: bool = False, all_transcripts: bool = False):
        self.ref_fasta_file = ref_fasta_file
        self.gff_file = gff_file
        self.vcf_file = vcf_file
        self.output_dir = output_dir
        self.ref_dict = fasta_to_dict(ref_fasta_file)
        self.cds_only = cds_only
        self.all_transcripts = all_transcripts

        v = VCF(vcf_file)
        self.samples = v.samples
        v.close()

        self.gene_info: Dict[str, Tuple[str, int, int, str]] = {}
        self.gene_transcripts = defaultdict(dict)
        self._parse_gff()
        
        os.makedirs(self.output_dir, exist_ok=True)

    def _parse_gff(self):
        logger.info("Parsing GFF structure...")
        gene_count = 0
        mrna_map = defaultdict(lambda: {'gene_id': '', 'strand': '', 'cds': []})

        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                chrom, _, feature, start_str, end_str, _, strand, _, attr_field = fields
                attrs = parse_attributes(attr_field)
                start, end = int(start_str), int(end_str)

                if feature == 'gene':
                    gid = attrs.get('ID')
                    if gid:
                        self.gene_info[gid] = (chrom, start, end, strand)
                        gene_count += 1
                elif feature == 'mRNA':
                    mid = attrs.get('ID')
                    pid = attrs.get('Parent')
                    if mid and pid:
                        mrna_map[mid]['gene_id'] = pid
                        mrna_map[mid]['strand'] = strand
                elif feature == 'CDS':
                    pid = attrs.get('Parent')
                    if pid and pid in mrna_map:
                        mrna_map[pid]['cds'].append((start, end))

        for mid, info in mrna_map.items():
            if info['cds']:
                gid = info['gene_id']
                info['cds'].sort(key=lambda x: x[0])
                self.gene_transcripts[gid][mid] = info

        logger.info(f"Successfully parsed {gene_count} genes")

    def prepare_gene_list(self, single_gene: str = None, gene_file: str = None) -> List[str]:
        genes = []

        if gene_file:
            if os.path.exists(gene_file):
                with open(gene_file, 'r') as f:
                    for line in f:
                        g = line.strip()
                        if g and not g.startswith('#'):
                            genes.append(g)
            else:
                logger.error(f"Gene list file does not exist: {gene_file}")

        elif single_gene:
            genes.append(single_gene)

        else:
            genes = list(self.gene_info.keys())
        return genes
        
    def get_gene_info(self, gene_name: str) -> Optional[Tuple[str, int, int, str]]:
        return self.gene_info.get(gene_name)

    def extract_region_phased(self, chrom: str, start: int, end: int) -> Dict[str, Tuple[str, str]]:
        if chrom not in self.ref_dict:
            logger.error(f"Chromosome {chrom} is not present in the reference genome")
            return {}

        full_ref_seq = self.ref_dict[chrom]
        # sample -> [hap1_parts, hap2_parts]
        sample_builders = {s: [[], []] for s in self.samples}

        try:
            if not os.path.exists(self.vcf_file + ".tbi") and not os.path.exists(self.vcf_file + ".csi"):
                 logger.error(f"VCF index missing: {self.vcf_file}.tbi does not exist. Please run 'tabix -p vcf' first.")
                 return {}

            region_str = f"{chrom}:{start}-{end}"
            vcf = VCF(self.vcf_file)
            current_ref_pos = start

            for v in vcf(region_str):
                if v.POS > current_ref_pos:
                    chunk = full_ref_seq[current_ref_pos - 1: v.POS - 1]
                    for s in self.samples:
                        sample_builders[s][0].append(chunk)
                        sample_builders[s][1].append(chunk)

                if v.POS < current_ref_pos:
                    continue

                ref = v.REF
                alts = v.ALT if v.ALT else []
                all_alleles = [ref] + alts

                max_len = max(len(a) for a in all_alleles)
                
                gts = v.genotypes
                
                for i, sample in enumerate(self.samples):
                    gt = gts[i]
                    idx1, idx2 = gt[0], gt[1]

                    seq1 = all_alleles[idx1] if idx1 != -1 else 'N' * max_len # 缺失处理为N
                    seq2 = all_alleles[idx2] if idx2 != -1 else 'N' * max_len

                    if len(seq1) < max_len:
                        seq1 += '-' * (max_len - len(seq1))
                    if len(seq2) < max_len:
                        seq2 += '-' * (max_len - len(seq2))
                        
                    sample_builders[sample][0].append(seq1)
                    sample_builders[sample][1].append(seq2)

                current_ref_pos = v.POS + len(ref)

            if current_ref_pos <= end:
                chunk = full_ref_seq[current_ref_pos - 1: end]
                for s in self.samples:
                    sample_builders[s][0].append(chunk)
                    sample_builders[s][1].append(chunk)

        except Exception as e:
            logger.error(f"VCF file reading or processing error: {e}")
            import traceback
            traceback.print_exc()
            return {}

        final_seqs = {}
        for s in self.samples:
            hap1 = "".join(sample_builders[s][0])
            hap2 = "".join(sample_builders[s][1])
            final_seqs[s] = (hap1, hap2)
            
        return final_seqs

    def _extract_single_transcript_cds_phased(self, chrom: str, t_info: dict) -> Dict[str, Tuple[str, str]]:
        cds_segments = t_info['cds']
        strand = t_info['strand']
        full_ref_seq = self.ref_dict[chrom]

        sample_builders = {s: [[], []] for s in self.samples}
        vcf = VCF(self.vcf_file)

        for seg_start, seg_end in cds_segments:
            region_str = f"{chrom}:{seg_start}-{seg_end}"
            current_ref_pos = seg_start

            try:
                records = vcf(region_str)
            except:
                records = []

            for v in records:
                if v.POS < current_ref_pos: continue

                if v.POS > current_ref_pos:
                    chunk = full_ref_seq[current_ref_pos - 1: v.POS - 1]
                    for s in self.samples:
                        sample_builders[s][0].append(chunk)
                        sample_builders[s][1].append(chunk)

                ref = v.REF
                alts = v.ALT if v.ALT else []
                all_alleles = [ref] + alts
                max_len = max(len(a) for a in all_alleles)
                gts = v.genotypes

                for i, sample in enumerate(self.samples):
                    gt = gts[i]
                    idx1, idx2 = gt[0], gt[1]

                    seq1 = all_alleles[idx1] if idx1 != -1 else 'N' * max_len
                    seq2 = all_alleles[idx2] if idx2 != -1 else 'N' * max_len

                    if len(seq1) < max_len:
                        seq1 += '-' * (max_len - len(seq1))
                    if len(seq2) < max_len:
                        seq2 += '-' * (max_len - len(seq2))

                    sample_builders[sample][0].append(seq1)
                    sample_builders[sample][1].append(seq2)

                current_ref_pos = v.POS + len(ref)

            if current_ref_pos <= seg_end:
                chunk = full_ref_seq[current_ref_pos - 1: seg_end]
                for s in self.samples:
                    sample_builders[s][0].append(chunk)
                    sample_builders[s][1].append(chunk)

        final_seqs = {}
        for s, parts in sample_builders.items():
            hap1 = "".join(parts[0])
            hap2 = "".join(parts[1])

            if strand == '-':
                hap1 = rev_comp(hap1)
                hap2 = rev_comp(hap2)

            final_seqs[s] = (hap1, hap2)

        return final_seqs

    def extract_gene(self, gene_name: str) -> bool:
        gene_info = self.get_gene_info(gene_name)
        
        if not gene_info:
            logger.warning(f"Gene {gene_name} not found in GFF file")
            return False
        
        chrom, start, end, strand = gene_info
        if self.cds_only:
            if gene_name not in self.gene_transcripts:
                logger.warning(f"Gene {gene_name} don't have CDS")
                return False

            transcripts = self.gene_transcripts[gene_name]

            if not self.all_transcripts:
                best_mid = max(transcripts.keys(), key=lambda k: sum(e - s + 1 for s, e in transcripts[k]['cds']))
                transcripts = {best_mid: transcripts[best_mid]}

            all_final_seqs = {}

            for mrna_id, t_info in transcripts.items():
                seqs_dict = self._extract_single_transcript_cds_phased(chrom, t_info)
                for sample, (hap1, hap2) in seqs_dict.items():
                    all_final_seqs[f"{sample}_{mrna_id}"] = (hap1, hap2)

            output_file = os.path.join(self.output_dir, f"{gene_name}_cds_phased.fasta")
            self._write_fasta(gene_name, all_final_seqs, output_file)
            logger.info(f"Gene {gene_name} CDS has been extracted (including {len(transcripts)} mRNA)")
        else:
            final_seqs = self.extract_region_phased(chrom, start, end)
            output_file = os.path.join(self.output_dir, f"{gene_name}_phased.fasta")
            self._write_fasta(gene_name, final_seqs, output_file)
            logger.info(f"Sequence for gene {gene_name} saved to: {output_file}")
        
        return True
    
    def extract_genes(self, gene_list: List[str]) -> Tuple[int, int]:
        success_count = 0
        fail_count = 0
        
        for gene_name in gene_list:
            if self.extract_gene(gene_name):
                success_count += 1
            else:
                fail_count += 1
        
        return success_count, fail_count

    def _write_fasta(self, gene_name: str, final_seqs: Dict[str, Tuple[str, str]], out_file: str):
        with open(out_file, 'w') as f:
            for header_suffix, (hap1, hap2) in final_seqs.items():
                f.write(f">{gene_name}_{header_suffix}_h1\n{hap1}\n")
                f.write(f">{gene_name}_{header_suffix}_h2\n{hap2}\n")

    def extract_all(self, single_gene=None, gene_file=None):
        target_genes = self.prepare_gene_list(single_gene, gene_file)

        if not target_genes:
            logger.warning("No genes to extract, task terminated.")
            return 0, 0
        return self.extract_genes(target_genes)


def run(args):
    extractor = PhasedExtractor(
        ref_fasta_file=args.fasta,
        vcf_file=args.vcf,
        gff_file=args.gff,
        output_dir=args.output_dir,
        cds_only=getattr(args, 'cds_only', False),
        all_transcripts=getattr(args, 'all_transcripts', False)
    )

    success_count, fail_count = extractor.extract_all(single_gene=args.gene, gene_file=args.gene_file)

    logger.info(f"Extraction completed! Success: {success_count}, Failed: {fail_count}")
