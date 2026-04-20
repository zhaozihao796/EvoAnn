import time
import logging
import sys
from . annotate_utils import parse_attributes, parse_vcf_all2

logger = logging.getLogger(__name__)


class StreamAnnotator:
    def __init__(self, gff_file, vcf_file, flank=2000, vcf_dict=None):
        self.gff_file = gff_file
        self.vcf_file = vcf_file
        self.flank = flank
        self.vcf_dict = vcf_dict if vcf_dict is not None else {}

        self.chrom_name = {}  # {gene_id: chrom}
        self.mrna_by_gene = {}  # {gene_id: [(mrna_id, start, end, strand), ...]}
        # self.vcf_dict = {}
        self.annotated_snps = {}  # {(gene, mrna, tag): [snps]}

    def _parse_gff(self):
        logger.info("Parsing GFF file...")
        with open(self.gff_file, 'r') as f:
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
                    self.mrna_by_gene.setdefault(parent, []).append(
                        (attrs.get('ID'), int(start), int(end), strand)
                    )

    def run_analysis(self):
        self._parse_gff()
        logger.info(f"GFF parsing completed | Found {len(self.mrna_by_gene)} genes in total")

        if not self.vcf_dict:
            logger.info("Starting to parse VCF file...")
            self.vcf_dict = parse_vcf_all2(self.vcf_file)
            for chrom in self.vcf_dict:
                self.vcf_dict[chrom].sort(key=lambda x: x[0])

        self._annotate_streams()

    def _annotate_streams(self):
        total_genes = len(self.mrna_by_gene)
        processed = 0
        total_hits = 0
        start_time = time.time()

        logger.info(f"Starting to annotate upstream and downstream regions | Range: {self.flank}bp...")

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

            mrna_list_sorted = sorted(mrna_list, key=lambda x: x[1])

            j = 0

            for mrna_id, start, end, strand in mrna_list_sorted:
                up_s = max(1, start - self.flank)
                up_e = start - 1
                down_s = end + 1
                down_e = end + self.flank

                if strand == '+':
                    range1_tag, range1_s, range1_e = 'upstream', up_s, up_e
                    range2_tag, range2_s, range2_e = 'downstream', down_s, down_e
                else:
                    range1_tag, range1_s, range1_e = 'downstream', up_s, up_e
                    range2_tag, range2_s, range2_e = 'upstream', down_s, down_e

                while j < n_snps and vlist[j][0] < range1_s:
                    j += 1

                curr = j
                while curr < n_snps and vlist[curr][0] <= range1_e:
                    self._add_snp(gene, mrna_id, range1_tag, vlist[curr], chrom)
                    total_hits += 1
                    curr += 1

                search_idx = curr
                while search_idx < n_snps and vlist[search_idx][0] < range2_s:
                    search_idx += 1

                while search_idx < n_snps and vlist[search_idx][0] <= range2_e:
                    self._add_snp(gene, mrna_id, range2_tag, vlist[search_idx], chrom)
                    total_hits += 1
                    search_idx += 1

        sys.stdout.write("\n")
        sys.stdout.flush()
        logger.info(f"Annotation completed | Total hits {total_hits} sites")

    def _add_snp(self, gene, mrna, tag, vcf_record, chrom):
        pos, ref, alt = vcf_record
        alt_val = alt if isinstance(alt, str) else alt[0]
        self.annotated_snps.setdefault((gene, mrna, tag), []).append(
            (chrom, pos, ref, alt_val)
        )

    def write_output(self, output_file):
        logger.info(f"Starting to write result file: {output_file}")
        with open(output_file, 'w') as f:
            for (gene, mrna, tag), snps in self.annotated_snps.items():
                for chrom, pos, ref, alt in snps:
                    f.write(f"{gene}\t{mrna}\t{tag}\t{chrom}\t{pos}\t{ref}\t{alt}\n")
        
        total_results = sum(len(snps) for snps in self.annotated_snps.values())
        logger.info(f"Result writing completed | Saved {total_results} upstream and downstream SNPs to {output_file}")


def run(args):
    annotator = StreamAnnotator(args.gff, args.vcf, args.flank)
    annotator.run_analysis()
    annotator.write_output(f"{args.output}.fa")