import logging
import time
import sys
from . annotate_utils import parse_attributes, parse_vcf_all2
from . GFFCorrector import repair_gff_missing_exons_and_phases

logger = logging.getLogger(__name__)


class UTRAnnotator:
    def __init__(self, gff_file, vcf_file, vcf_dict=None):
        self.gff_file = gff_file
        self.vcf_file = vcf_file
        self.vcf_dict = vcf_dict if vcf_dict is not None else {}

        self.chrom_name = {}  # {gene_id: chrom}
        self.mrna_by_gene = {}  # {gene: [(mrna, strand), ...]}
        self.utr_by_mrna = {}  # {mrna: [(start, end, tag), ...]}
        # self.vcf_dict = {}
        self.annotated_snps = {}  # {(gene, mrna, tag): [snps]}

    def _parse_and_infer_utrs(self):
        logger.info("Preprocessing GFF...")
        repaired_gff = repair_gff_missing_exons_and_phases(self.gff_file)

        exon_by_mrna = {}
        cds_by_mrna = {}

        logger.info("Parsing GFF structure...")
        with open(repaired_gff, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\t')
                if len(fields) < 9: continue

                chrom, _, feature, start, end, _, strand, _, attr_field = fields
                attrs = parse_attributes(attr_field)

                if feature == 'gene':
                    self.chrom_name[attrs.get('ID')] = chrom
                elif feature == 'mRNA':
                    parent = attrs.get('Parent')
                    self.mrna_by_gene.setdefault(parent, []).append((attrs.get('ID'), strand))
                elif feature == 'exon':
                    exon_by_mrna.setdefault(attrs.get('Parent'), []).append((int(start), int(end)))
                elif feature == 'CDS':
                    cds_by_mrna.setdefault(attrs.get('Parent'), []).append((int(start), int(end)))

        for gene, mrna_list in self.mrna_by_gene.items():
            for mrna_id, strand in mrna_list:
                exons = sorted(exon_by_mrna.get(mrna_id, []), key=lambda x: x[0])
                cdss = sorted(cds_by_mrna.get(mrna_id, []), key=lambda x: x[0])

                if not exons or not cdss: continue

                if strand == '-':
                    tag_5, tag_3 = '3_prime_UTR_variant', '5_prime_UTR_variant'
                else:
                    tag_5, tag_3 = '5_prime_UTR_variant', '3_prime_UTR_variant'

                i, j = 0, 0
                while i < len(exons) and j < len(cdss):
                    e_start, e_end = exons[i]
                    c_start, c_end = cdss[j]

                    if e_end < c_start:
                        self.utr_by_mrna.setdefault(mrna_id, []).append((e_start, e_end, tag_5))
                        i += 1
                    elif e_start > c_end:
                        self.utr_by_mrna.setdefault(mrna_id, []).append((e_start, e_end, tag_3))
                        i += 1

                    else:
                        if e_start < c_start:
                            self.utr_by_mrna.setdefault(mrna_id, []).append((e_start, c_start - 1, tag_5))

                        if e_end > c_end:
                            self.utr_by_mrna.setdefault(mrna_id, []).append((c_end + 1, e_end, tag_3))

                        if e_end >= c_end:
                            j += 1
                        i += 1

                while i < len(exons):
                    e_start, e_end = exons[i]
                    if j >= len(cdss):
                        self.utr_by_mrna.setdefault(mrna_id, []).append((e_start, e_end, tag_3))
                    i += 1

        count = sum(len(u) for u in self.utr_by_mrna.values())
        logger.info(f"GFF parsing completed | Identified {count} UTR segments in total")

    def run_analysis(self):
        self._parse_and_infer_utrs()

        if not self.vcf_dict:
            logger.info("Starting to parse VCF file...")
            self.vcf_dict = parse_vcf_all2(self.vcf_file)
            for chrom in self.vcf_dict:
                self.vcf_dict[chrom].sort(key=lambda x: x[0])

        self._annotate_utrs()

    def _annotate_utrs(self):
        utr_by_chrom = {}
        for gene, mrna_list in self.mrna_by_gene.items():
            chrom = self.chrom_name.get(gene)
            if not chrom: continue

            for mrna_id, _ in mrna_list:
                utrs = self.utr_by_mrna.get(mrna_id, [])
                for start, end, tag in utrs:
                    utr_by_chrom.setdefault(chrom, []).append(
                        (start, end, tag, gene, mrna_id)
                    )

        total_hits = 0
        total_chroms = len(utr_by_chrom)
        processed_chroms = 0
        start_time = time.time()

        logger.info(f"Starting to annotate UTR variants...")

        for chrom, utr_regions in utr_by_chrom.items():
            processed_chroms += 1
            
            if chrom not in self.vcf_dict:
                sys.stdout.write(f"\rAnnotation progress: Chromosome {processed_chroms}/{total_chroms} | Hits:{total_hits}   ")
                sys.stdout.flush()
                continue

            utr_regions.sort(key=lambda x: x[0])
            snp_list = self.vcf_dict[chrom]

            n_snps = len(snp_list)
            snp_idx = 0

            for u_start, u_end, u_tag, gene, mrna in utr_regions:
                while snp_idx < n_snps and snp_list[snp_idx][0] < u_start:
                    snp_idx += 1

                k = snp_idx
                while k < n_snps and snp_list[k][0] <= u_end:
                    pos, ref, alt = snp_list[k]

                    alt_val = alt[0] if isinstance(alt, list) else alt
                    self.annotated_snps.setdefault((gene, mrna, u_tag), []).append(
                        (chrom, pos, ref, alt_val)
                    )
                    total_hits += 1
                    k += 1

            sys.stdout.write(f"\rAnnotation progress: Chromosome {processed_chroms}/{total_chroms} | Hits:{total_hits}   ")
            sys.stdout.flush()

        sys.stdout.write("\n")
        sys.stdout.flush()
        logger.info(f"Annotation completed | Total hits {total_hits} UTR variants")

    def write_output(self, output_file):
        logger.info(f"Starting to write result file: {output_file}")
        with open(output_file, 'w') as f:
            for (gene, mrna, tag), snps in self.annotated_snps.items():
                for chrom, pos, ref, alt in snps:
                    f.write(f"{gene}\t{mrna}\t{tag}\t{chrom}\t{pos}\t{ref}\t{alt}\n")
        
        total_results = sum(len(snps) for snps in self.annotated_snps.values())
        logger.info(f"Result writing completed | Saved {total_results} UTR variants to {output_file}")


def run(args):
    annotator = UTRAnnotator(args.gff, args.vcf)
    annotator.run_analysis()
    annotator.write_output(f"{args.output}.fa")