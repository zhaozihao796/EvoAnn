import time
import logging
import sys
from . annotate_utils import parse_attributes, parse_vcf_all2
from . GFFCorrector import repair_gff_missing_exons_and_phases

logger = logging.getLogger(__name__)


class IntronAnnotator:
    def __init__(self, gff_file, vcf_file, vcf_dict=None):
        self.gff_file = gff_file
        self.vcf_file = vcf_file
        self.vcf_dict = vcf_dict if vcf_dict is not None else {}

        self.chrom_name = {}  # {gene_id: chrom}
        self.mrna_by_gene = {}  # {gene_id: [(mrna_id, start, end), ...]}
        self.intron_by_mrna = {}  # {mrna_id: [(intron_id, start, end), ...]}
        # self.vcf_dict = {}
        self.annotated_snps = {}  # {(gene, mrna, intron): [(chrom, pos, ref, alt), ...]}
        self.mrna_info = {}

    def _get_intron_info(self):
        exon_by_mrna = {}

        logger.info("Preprocessing GFF file")
        repaired_gff = repair_gff_missing_exons_and_phases(self.gff_file)

        logger.info("Parsing GFF file...")
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
                    self.mrna_info[attrs.get('ID')] = strand
                    self.mrna_by_gene.setdefault(parent, []).append(
                        (attrs.get('ID'), int(start), int(end))
                    )
                elif feature == 'exon':
                    parent = attrs.get('Parent')
                    exon_by_mrna.setdefault(parent, []).append(
                        (attrs.get('ID'), int(start), int(end))
                    )

        for gene, mrna_list in self.mrna_by_gene.items():
            for mrna_id, _, _ in mrna_list:
                exon_list = exon_by_mrna.get(mrna_id, [])
                if not exon_list: continue

                exon_sorted = sorted(exon_list, key=lambda x: x[1])

                merged_exons = []
                for _, s, e in exon_sorted:
                    if not merged_exons or s > merged_exons[-1][1]:
                        merged_exons.append([s, e])
                    else:
                        merged_exons[-1][1] = max(merged_exons[-1][1], e)

                for i in range(len(merged_exons) - 1):
                    intron_start = merged_exons[i][1] + 1
                    intron_end = merged_exons[i + 1][0] - 1

                    if intron_start <= intron_end:
                        intron_id = f"{mrna_id}_intron_{i + 1}"
                        self.intron_by_mrna.setdefault(mrna_id, []).append(
                            (intron_id, intron_start, intron_end)
                        )

    def run_analysis(self):
        self._get_intron_info()
        logger.info(
            f"GFF parsing completed | {len(self.mrna_by_gene)} genes | {sum(len(v) for v in self.intron_by_mrna.values())} introns")

        if not self.vcf_dict:
            logger.info("Starting to parse VCF file...")
            self.vcf_dict = parse_vcf_all2(self.vcf_file)
            for chrom in self.vcf_dict:
                self.vcf_dict[chrom].sort(key=lambda x: x[0])

        self._annotate_introns()

    def _annotate_introns(self):
        total_genes = len(self.mrna_by_gene)
        processed = 0
        total_hits = 0
        start_time = time.time()

        logger.info(f"Starting to annotate intron SNPs...")

        for gene, mrna_list in self.mrna_by_gene.items():
            processed += 1

            if processed % 100 == 0 or processed == total_genes:
                elapsed = time.time() - start_time
                avg = elapsed / processed
                rem_sec = (total_genes - processed) * avg
                if rem_sec < 60:
                    rem_str = f"{rem_sec:.0f}s"
                elif rem_sec < 3600:
                    rem_str = f"{rem_sec/60:.1f}min"
                else:
                    rem_str = f"{rem_sec/3600:.1f}h"
                percent = (processed / total_genes) * 100
                sys.stdout.write(f"\rAnnotation progress: {processed}/{total_genes} [{percent:3.0f}%] | Hits:{total_hits} | Speed:{avg:.2f}s/each | Remaining:{rem_str}   ")
                sys.stdout.flush()

            chrom = self.chrom_name.get(gene)
            if not chrom or chrom not in self.vcf_dict: continue

            vlist = self.vcf_dict[chrom]
            n_snps = len(vlist)

            for mrna_id, _, _ in mrna_list:
                strand = self.mrna_info.get(mrna_id)
                introns = self.intron_by_mrna.get(mrna_id, [])
                if not introns: continue

                introns_sorted = sorted(introns, key=lambda x: x[1])

                snp_idx = 0
                for intron_id, start, end in introns_sorted:
                    while snp_idx < n_snps and vlist[snp_idx][0] < start:
                        snp_idx += 1

                    curr = snp_idx
                    while curr < n_snps and vlist[curr][0] <= end:
                        pos, ref, alt = vlist[curr]
                        alt_val = alt if isinstance(alt, str) else alt[0]

                        tag = "intron_variant"
                        intron_len = end - start + 1
                        if intron_len >= 4:
                            if strand == '+':
                                if pos <= start + 1:
                                    tag = "splice_donor_variant"
                                elif pos >= end - 1:
                                    tag = "splice_acceptor_variant"
                            else:
                                if pos >= end - 1:
                                    tag = "splice_donor_variant"
                                elif pos <= start + 1:
                                    tag = "splice_acceptor_variant"

                        key = (gene, mrna_id, intron_id, tag)
                        self.annotated_snps.setdefault(key, []).append(
                            (chrom, pos, ref, alt_val)
                        )
                        total_hits += 1
                        curr += 1
                    snp_idx = curr

        sys.stdout.write("\n")
        sys.stdout.flush()
        logger.info(f"Annotation completed | Total hits {total_hits} intron SNPs")

    def write_output(self, output_file):
        logger.info(f"Starting to write result file: {output_file}")
        with open(output_file, 'w') as f:
            for (gene, mrna, intron), snps in self.annotated_snps.items():
                for chrom, pos, ref, alt in snps:
                    f.write(f"{gene}\t{mrna}\t{intron}\t{chrom}\t{pos}\t{ref}\t{alt}\n")
        
        total_results = sum(len(snps) for snps in self.annotated_snps.values())
        logger.info(f"Result writing completed | Saved {total_results} intron SNPs to {output_file}")


def run(args):
    annotator = IntronAnnotator(args.gff, args.vcf)
    annotator.run_analysis()
    annotator.write_output(f"{args.output}.fa")