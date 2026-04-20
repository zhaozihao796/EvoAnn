import logging
import sys
import time
from . annotate_utils import parse_attributes, rev_comp, parse_cds_fasta, parse_vcf_all2

logger = logging.getLogger(__name__)

CODON2AA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
START_CODONS = {'ATG', 'CTG', 'GTG', 'TTG', 'ATT', 'ATC', 'ATA'}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


class CDSAnnotator:
    def __init__(self, gff_file, vcf_file, cds_file, vcf_dict=None):
        self.gff_file = gff_file
        self.vcf_file = vcf_file
        self.cds_file = cds_file
        self.vcf_dict = vcf_dict if vcf_dict is not None else {}

        self.cds_by_gene = {}  # {gene: {mrna: [cds_tuple]}}
        self.chrom_name = {}  # {gene: chrom}
        # self.vcf_dict = {}
        self.transcripts = {}  # {(gene, mrna): seq}
        self.nsmutation = {}  # {(gene, mrna): [mutations]}

    def _parse_gff_structure(self):
        cds_by_mrna = {}
        mrna_by_gene = {}

        with open(self.gff_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\t')
                if len(fields) < 9: continue

                chrom, _, feature, start, end, _, strand, phase, attrs_str = fields
                attrs = parse_attributes(attrs_str)

                if feature == 'gene':
                    self.chrom_name[attrs.get('ID')] = chrom
                elif feature == 'mRNA':
                    mrna_by_gene.setdefault(attrs.get('Parent'), []).append(attrs.get('ID'))
                elif feature == 'CDS':
                    parent = attrs.get('Parent')
                    phase_int = int(phase) if phase.isdigit() else 0
                    cds_by_mrna.setdefault(parent, []).append(
                        (int(start), int(end), strand, chrom, phase_int)
                    )

        for gene_id, mrna_list in mrna_by_gene.items():
            valid_mrnas = {}
            for mid in mrna_list:
                if mid in cds_by_mrna:
                    valid_mrnas[mid] = cds_by_mrna[mid]

            if valid_mrnas:
                self.cds_by_gene[gene_id] = valid_mrnas

    def _get_longest_mrna_id(self, mrna_map):
        return max(mrna_map.items(), key=lambda kv: sum(e - s + 1 for s, e, _, _, _ in kv[1]))[0]

    def _translate(self, codon):
        if len(codon) != 3: return None
        c = codon.upper()
        if any(x in c for x in ['N', '-', 'B', 'D', 'H', 'V']): return None
        return CODON2AA.get(c)

    def _determine_variant_type(self, original, mutated, pos_in_transcript):
        orig_U = original.upper()
        mut_U = mutated.upper()

        if pos_in_transcript == 0:
            if orig_U in START_CODONS and mut_U not in START_CODONS:
                return 'start_lost'

        orig_stop = orig_U in STOP_CODONS
        mut_stop = mut_U in STOP_CODONS

        if orig_stop and not mut_stop: return 'stop_lost'
        if not orig_stop and mut_stop: return 'stop_gained'

        return 'missense_variant'

    def run_analysis(self):
        logger.info("Starting to parse GFF file...")
        self._parse_gff_structure()
        logger.info(f"GFF parsing completed | Associated with {len(self.cds_by_gene)} genes")

        if not self.vcf_dict:
            logger.info("Starting to parse VCF file...")
            self.vcf_dict = parse_vcf_all2(self.vcf_file)
            for chrom in self.vcf_dict:
                self.vcf_dict[chrom].sort(key=lambda x: x[0])
            logger.info(f"VCF parsing completed | Contains {sum(len(v) for v in self.vcf_dict.values())} variant sites")

        logger.info("Starting to match CDS sequences with structure...")
        raw_cds_dict, _ = parse_cds_fasta(self.cds_file, validate=True)

        cds_mapping_filtered = {}

        for header, seq in raw_cds_dict.items():
            if '|' in header:
                gene, mrna_id = header.split('|', 1)
                if gene in self.cds_by_gene and mrna_id in self.cds_by_gene[gene]:
                    self.transcripts[(gene, mrna_id)] = seq
                    cds_mapping_filtered.setdefault(gene, {})[mrna_id] = self.cds_by_gene[gene][mrna_id]
            else:
                gene = header
                if gene in self.cds_by_gene:
                    best_mrna = self._get_longest_mrna_id(self.cds_by_gene[gene])
                    self.transcripts[gene] = seq
                    cds_mapping_filtered[gene] = self.cds_by_gene[gene][best_mrna]

        logger.info(f"CDS sequence matching completed | Matched {len(self.transcripts)} transcripts")
        logger.info("Starting to identify non-synonymous mutations...")
        self._annotate_variants(cds_mapping_filtered)

    def _annotate_variants(self, cds_mapping_filtered):
        total = len(self.transcripts)
        count = 0
        start_time = time.time()

        for key, seq in self.transcripts.items():
            count += 1

            if count % 100 == 0 or count == total:
                elapsed = time.time() - start_time
                avg = elapsed / count
                rem_sec = (total - count) * avg
                if rem_sec < 60:
                    rem_str = f"{rem_sec:.0f}s"
                elif rem_sec < 3600:
                    rem_str = f"{rem_sec/60:.1f}min"
                else:
                    rem_str = f"{rem_sec/3600:.1f}h"
                percent = (count / total) * 100
                sys.stdout.write(f"\rIdentification progress: {count}/{total} [{percent:3.0f}%] | Speed:{avg:.2f}s/each | Remaining:{rem_str}   ")
                sys.stdout.flush()

            if isinstance(key, tuple):
                gene, mrna_id = key
            else:
                gene = key
                mrna_id = None

            chrom = self.chrom_name.get(gene)
            if not chrom or chrom not in self.vcf_dict: continue
            if mrna_id:
                cds_list = cds_mapping_filtered.get(gene, {}).get(mrna_id)
            else:
                cds_list = cds_mapping_filtered.get(gene)

            if not cds_list: continue

            cds_sorted = sorted(cds_list, key=lambda x: x[0])

            vlist = self.vcf_dict[chrom]

            transcript_offset = 0
            transcript_len = sum(c[1] - c[0] + 1 for c in cds_list)
            snp_idx = 0

            for start, end, strand, _, _ in cds_sorted:
                while snp_idx < len(vlist) and vlist[snp_idx][0] < start:
                    snp_idx += 1

                curr = snp_idx
                while curr < len(vlist) and vlist[curr][0] <= end:
                    pos, ref, alt = vlist[curr]

                    if strand == '+':
                        tpos = pos - start + transcript_offset
                    else:
                        pre_index = transcript_offset + (pos - start)
                        tpos = transcript_len - pre_index - 1

                        ref = rev_comp(ref)
                        alt = rev_comp(alt)

                    if 0 <= tpos < len(seq) and seq[tpos] == ref:
                        frame = tpos % 3
                        codon_start = tpos - frame
                        codon_end = codon_start + 3

                        if codon_end <= len(seq):
                            original_codon = seq[codon_start:codon_end]

                            prefix = original_codon[:frame]
                            suffix = original_codon[frame + 1:]
                            mutated_codon = prefix + alt + suffix

                            if len(mutated_codon) == 3:
                                aa_ref = self._translate(original_codon)
                                aa_alt = self._translate(mutated_codon)

                                if aa_ref and aa_alt and aa_ref != aa_alt:
                                    aa_pos = codon_start // 3 + 1
                                    var_type = self._determine_variant_type(original_codon, mutated_codon, aa_pos - 1)
                                    change_str = f"{aa_ref}{aa_pos}{aa_alt}"

                                    res_key = (gene, mrna_id) if mrna_id else gene
                                    raw_pos, raw_ref, raw_alt = vlist[curr]

                                    self.nsmutation.setdefault(res_key, []).append(
                                        (raw_pos, raw_ref, raw_alt, var_type, change_str)
                                    )

                    curr += 1
                transcript_offset += (end - start + 1)
        
        sys.stdout.write("\n")
        sys.stdout.flush()

    def write_output(self, output_file):
        logger.info(f"Starting to write result file: {output_file}")
        with open(output_file, 'w') as out:
            for k, muts in self.nsmutation.items():
                gene = k[0] if isinstance(k, tuple) else k
                chrom = self.chrom_name.get(gene, "Unknown")

                for pos, ref, alt, vtype, change in muts:
                    out.write(f">{chrom}\t{pos}\t{ref}\t{alt}\t{vtype}\t{gene}\t{change}\n")

        total_mutations = sum(len(m) for m in self.nsmutation.values())
        logger.info(f"Analysis completed | Found {total_mutations} non-synonymous mutation sites | Results saved to {output_file}")


def run(args):
    annotator = CDSAnnotator(args.gff, args.vcf, args.cds)
    annotator.run_analysis()
    annotator.write_output(args.output)
