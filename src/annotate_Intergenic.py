import sys
import time
import logging
import bisect
from . annotate_utils import parse_attributes, parse_vcf_all2

logger = logging.getLogger(__name__)

class IntergenicAnnotator:
    def __init__(self, gff_file, vcf_file, output_file, vcf_dict=None):
        self.gff_file = gff_file
        self.vcf_file = vcf_file
        self.output_file = output_file
        self.vcf_dict = vcf_dict if vcf_dict is not None else {}
        self.genes_by_chrom = {}
        self.precomputed_max_ends = {}
        
        self.vcf_dict = {}
        self.results = []

    def _parse_gff(self):
        logger.info("Parsing GFF file for genes...")
        count = 0
        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\t')
                if len(fields) < 9: continue

                chrom, _, feature, start, end, _, strand, _, attr_field = fields
                
                if feature == 'gene':
                    attrs = parse_attributes(attr_field)
                    gene_id = attrs.get('ID')
                    if gene_id:
                        self.genes_by_chrom.setdefault(chrom, []).append({
                            'start': int(start),
                            'end': int(end),
                            'id': gene_id
                        })
                        count += 1

        logger.info(f"Loaded {count} genes. Sorting and precomputing indices...")
        for chrom, genes in self.genes_by_chrom.items():
            genes.sort(key=lambda x: x['start'])

            max_ends = []
            current_max_end = -1
            current_max_gene = None
            
            for gene in genes:
                if gene['end'] > current_max_end:
                    current_max_end = gene['end']
                    current_max_gene = gene
                
                max_ends.append({
                    'end': current_max_end,
                    'id': current_max_gene['id']
                })
            self.precomputed_max_ends[chrom] = max_ends

    def run_analysis(self):
        self._parse_gff()
        
        if not self.vcf_dict:
            logger.info("Parsing VCF file...")
            self.vcf_dict = parse_vcf_all2(self.vcf_file)
        
        logger.info("Identifying intergenic variants...")
        total_snps = 0
        intergenic_snps = 0
        
        start_time = time.time()
        
        for chrom, snps in self.vcf_dict.items():
            if chrom not in self.genes_by_chrom:
                continue
                
            genes = self.genes_by_chrom[chrom]
            max_ends = self.precomputed_max_ends[chrom]
            gene_starts = [g['start'] for g in genes]
            
            for pos, ref, alt in snps:
                total_snps += 1

                idx = bisect.bisect_right(gene_starts, pos)
                
                left_gene = None
                right_gene = None
                is_intergenic = True

                if idx > 0:
                    left_max_info = max_ends[idx - 1]
                    
                    if left_max_info['end'] >= pos:
                        is_intergenic = False
                    else:
                        left_gene = {
                            'id': left_max_info['id'],
                            'end': left_max_info['end']
                        }

                if not is_intergenic:
                    continue

                if idx < len(genes):
                    right_gene = {
                        'id': genes[idx]['id'],
                        'start': genes[idx]['start']
                    }
                
                self.results.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'left_gene': left_gene['id'] if left_gene else '.',
                    'left_end': str(left_gene['end']) if left_gene else '.',
                    'right_gene': right_gene['id'] if right_gene else '.',
                    'right_start': str(right_gene['start']) if right_gene else '.'
                })
                intergenic_snps += 1
                
        elapsed = time.time() - start_time
        logger.info(f"Analysis finished in {elapsed:.2f}s. Found {intergenic_snps} intergenic SNPs out of {total_snps} total.")

    def write_output(self):
        logger.info(f"Writing results to {self.output_file}")
        with open(self.output_file, 'w') as f:
            f.write("Chrom\tPos\tRef\tAlt\tLeft_Gene\tLeft_Gene_End\tRight_Gene\tRight_Gene_Start\n")
            for res in self.results:
                f.write(f"{res['chrom']}\t{res['pos']}\t{res['ref']}\t{res['alt']}\t"
                        f"{res['left_gene']}\t{res['left_end']}\t"
                        f"{res['right_gene']}\t{res['right_start']}\n")

def run(args):
    annotator = IntergenicAnnotator(args.gff, args.vcf, args.output)
    annotator.run_analysis()
    annotator.write_output()
