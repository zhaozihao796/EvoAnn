import os
import logging
from . annotate_utils import parse_attributes, rev_comp, _check_cds, fasta_to_dict

logger = logging.getLogger(__name__)


class CDSExtractor:
    def __init__(self, gff_file, fasta_file):
        self.gff_file = gff_file
        self.fasta_file = fasta_file
        self.fasta_dict = {}

        self.cds_by_mrna = {}  # {mrna_id: [(start, end, strand, chrom, phase), ...]}
        self.mrna_by_gene = {}  # {gene_id: [mrna_id1, mrna_id2, ...]}
        self.chrom_map = {}  # {gene_id: chrom}

        self.transcripts = {}  # {(gene_id, mrna_id): seq}

    def load_data(self):
        logger.info("Parsing GFF file...")
        self._parse_gff()

        logger.info("Parsing FASTA file...")
        self.fasta_dict = fasta_to_dict(self.fasta_file)
        logger.info(f"Parsing completed: {len(self.fasta_dict)} chromosome sequences")

    def _parse_gff(self):
        with open(self.gff_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line or line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                chrom, _, feature, start, end, _, strand, phase, attrs_str = fields
                attrs = parse_attributes(attrs_str)

                if feature == 'gene':
                    gene_id = attrs.get('ID')
                    self.chrom_map[gene_id] = chrom

                elif feature == 'mRNA':
                    mrna_id = attrs.get('ID')
                    gene_id = attrs.get('Parent')
                    self.mrna_by_gene.setdefault(gene_id, []).append(mrna_id)

                elif feature == 'CDS':
                    mrna_id = attrs.get('Parent')
                    phase_int = int(phase) if phase.isdigit() else 0
                    self.cds_by_mrna.setdefault(mrna_id, []).append(
                        (int(start), int(end), strand, chrom, phase_int)
                    )

    def extract(self, longest_only=True):
        self.transcripts.clear()
        logger.info(f"Starting to extract CDS: ({'longest transcript' if longest_only else 'all transcripts'})...")

        for gene_id, mrna_list in self.mrna_by_gene.items():
            chrom = self.chrom_map.get(gene_id)
            if not chrom or chrom not in self.fasta_dict:
                continue

            target_mrnas = []
            if longest_only:
                max_len = -1
                best_mrna = None
                for mid in mrna_list:
                    if mid in self.cds_by_mrna:
                        current_len = sum(c[1] - c[0] + 1 for c in self.cds_by_mrna[mid])
                        if current_len > max_len:
                            max_len = current_len
                            best_mrna = mid
                if best_mrna:
                    target_mrnas = [best_mrna]
            else:
                target_mrnas = [mid for mid in mrna_list if mid in self.cds_by_mrna]

            for mrna_id in target_mrnas:
                seq = self._assemble_cds(mrna_id, chrom)
                if seq:
                    has_atg, has_stop, in_frame = _check_cds(seq)
                    if not has_atg: logger.warning(f"{gene_id}|{mrna_id} start is not ATG")
                    if not has_stop: logger.warning(f"{gene_id}|{mrna_id} no stop codon")
                    if not in_frame: logger.warning(f"{gene_id}|{mrna_id} length not multiple of 3: {len(seq)}")

                    self.transcripts[(gene_id, mrna_id)] = seq

    def _assemble_cds(self, mrna_id, chrom):
        cds_list = self.cds_by_mrna.get(mrna_id)
        if not cds_list:
            return ""

        strand = cds_list[0][2]
        if strand == '+':
            sorted_cds = sorted(cds_list, key=lambda x: x[0])
        else:
            sorted_cds = sorted(cds_list, key=lambda x: x[0], reverse=True)

        chrom_seq = self.fasta_dict[chrom]
        seq_parts = []

        for start, end, _, _, _ in sorted_cds:
            sub_seq = chrom_seq[max(0, start - 1): end]
            if strand == '-':
                sub_seq = rev_comp(sub_seq)
            seq_parts.append(sub_seq)

        return "".join(seq_parts)

    def write_fasta(self, output_file):
        logger.info(f"Writing output file: {output_file}")
        with open(output_file, 'w') as out:
            for (gene_id, mrna_id), seq in self.transcripts.items():
                header = f"{gene_id}|{mrna_id}"
                out.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")
        logger.info("Output completed.")


def run(args):
    extractor = CDSExtractor(gff_file=args.gff, fasta_file=args.fasta)

    extractor.load_data()

    longest_only = not getattr(args, 'all_transcripts', False)
    extractor.extract(longest_only=longest_only)

    output_path = f"{args.output}.fa"
    if not longest_only:
        output_path = f"{args.output}_all.fa"

    extractor.write_fasta(output_path)