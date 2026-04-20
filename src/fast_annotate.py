import sys
import os
import time
import logging
import gzip
import bisect
from collections import defaultdict
from .annotate_utils import parse_attributes, rev_comp, fasta_to_dict
from .GFFCorrector import repair_gff_missing_exons_and_phases

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
START_CODONS = {'ATG', 'CTG', 'TTG'}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

BIN_SIZE = 5000


class Transcript:
    __slots__ = ('id', 'gene_id', 'chrom', 'strand', 'start', 'end', 'exons', 'cds', 'utrs', 'cds_seq', 'cds_segments',
                 '_cds_min', '_cds_max')

    def __init__(self, tid, gid, chrom, strand, start, end):
        self.id = tid
        self.gene_id = gid
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end
        self.exons = []
        self.cds = []
        self.utrs = []
        self.cds_seq = None
        self.cds_segments = []
        self._cds_min = None
        self._cds_max = None

    def build_structure(self, ref_seq_chrom):
        if not self.cds:
            return

        if self.strand == '-':
            sorted_cds = sorted(self.cds, key=lambda x: x[0], reverse=True)
        else:
            sorted_cds = sorted(self.cds, key=lambda x: x[0])

        seq_parts = []
        cum_offset = 0

        for start, end, phase in sorted_cds:
            seg_len = end - start + 1

            part = ref_seq_chrom[start - 1:end]
            if self.strand == '-':
                part = rev_comp(part)

            seq_parts.append(part)

            self.cds_segments.append((start, end, cum_offset))

            cum_offset += seg_len

        self.cds_seq = "".join(seq_parts)


class GenomeIndex:
    def __init__(self):
        self.bins = defaultdict(lambda: defaultdict(list))
        self.genes = {}
        self.transcripts = {}
        self.chrom_genes = defaultdict(list)

    def add_feature(self, chrom, start, end, ftype, data):
        start_bin = start // BIN_SIZE
        end_bin = end // BIN_SIZE

        for b in range(start_bin, end_bin + 1):
            self.bins[chrom][b].append((start, end, ftype, data))

    def query(self, chrom, pos):
        bin_id = pos // BIN_SIZE
        hits = []
        if chrom in self.bins and bin_id in self.bins[chrom]:
            for start, end, ftype, data in self.bins[chrom][bin_id]:
                if start <= pos <= end:
                    hits.append((ftype, data))
        return hits

    def finalize(self):
        self.chrom_gene_starts = {}
        self.chrom_gene_max_ends = {}

        for chrom in self.chrom_genes:
            genes = self.chrom_genes[chrom]
            genes.sort(key=lambda x: x[0])

            starts = []
            max_ends = []
            current_max_end = -1
            current_max_gid = None

            for start, end, gid in genes:
                starts.append(start)
                if end > current_max_end:
                    current_max_end = end
                    current_max_gid = gid
                max_ends.append((current_max_end, current_max_gid))

            self.chrom_gene_starts[chrom] = starts
            self.chrom_gene_max_ends[chrom] = max_ends


class FastAnnotator:
    def __init__(self, vcf_file, gff_file, fasta_file, output_dir, flank=2000):
        self.vcf_file = vcf_file
        self.gff_file = gff_file
        self.fasta_file = fasta_file
        self.output_dir = output_dir
        self.flank = flank

        self.index = GenomeIndex()
        self.ref_genome = {}

        self.handles = {}
        self.counts = defaultdict(int)
        self.total_processed_variants = 0

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def load_genome(self):
        logger.info("Loading reference genome...")
        self.ref_genome = fasta_to_dict(self.fasta_file)

    def parse_gff(self):
        logger.info("Preprocessing GFF (repair missing exons and fix phases)...")
        repaired_gff = repair_gff_missing_exons_and_phases(self.gff_file)

        logger.info("Parsing GFF and building index...")

        gene_map = {}
        mrna_map = {}

        with open(repaired_gff, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\t')
                if len(fields) < 9: continue

                chrom, _, feature, start, end, _, strand, phase, attr_str = fields
                start, end = int(start), int(end)
                attrs = parse_attributes(attr_str)

                if feature == 'gene':
                    gid = attrs.get('ID')
                    if gid:
                        self.index.chrom_genes[chrom].append((start, end, gid))
                        gene_map[gid] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}

                elif feature == 'mRNA':
                    mid = attrs.get('ID')
                    pid = attrs.get('Parent')
                    if mid and pid:
                        t = Transcript(mid, pid, chrom, strand, start, end)
                        mrna_map[mid] = t
                        self.index.transcripts[mid] = t

                        if strand == '-':
                            self.index.add_feature(chrom, end + 1, end + self.flank, 'upstream', (pid, mid))
                            self.index.add_feature(chrom, max(1, start - self.flank), start - 1, 'downstream',
                                                   (pid, mid))
                        else:
                            self.index.add_feature(chrom, max(1, start - self.flank), start - 1, 'upstream', (pid, mid))
                            self.index.add_feature(chrom, end + 1, end + self.flank, 'downstream', (pid, mid))

                elif feature == 'CDS':
                    pid = attrs.get('Parent')
                    if pid and pid in mrna_map:
                        p_int = int(phase) if phase.isdigit() else 0
                        mrna_map[pid].cds.append((start, end, p_int))
                        self.index.add_feature(chrom, start, end, 'CDS', pid)

                elif feature == 'exon':
                    pid = attrs.get('Parent')
                    if pid and pid in mrna_map:
                        mrna_map[pid].exons.append((start, end))

        logger.info("Building transcript structures...")
        for mid, t in mrna_map.items():
            if t.chrom not in self.ref_genome: continue

            t.build_structure(self.ref_genome[t.chrom])

            # ---- Intron ----
            t.exons.sort(key=lambda x: x[0])
            merged_exons = []
            for es, ee in t.exons:
                if not merged_exons or es > merged_exons[-1][1]:
                    merged_exons.append([es, ee])
                else:
                    merged_exons[-1][1] = max(merged_exons[-1][1], ee)

            for i in range(len(merged_exons) - 1):
                intron_start = merged_exons[i][1] + 1
                intron_end = merged_exons[i + 1][0] - 1
                if intron_start <= intron_end:
                    intron_name = f"{mid}_intron_{i + 1}"
                    self.index.add_feature(t.chrom, intron_start, intron_end, 'intron', (mid, intron_name))
                    if intron_end - intron_start + 1 >= 4:
                        if t.strand == '+':
                            self.index.add_feature(t.chrom, intron_start, intron_start + 1, 'splice_donor', mid)
                            self.index.add_feature(t.chrom, intron_end - 1, intron_end, 'splice_acceptor', mid)
                        else:
                            self.index.add_feature(t.chrom, intron_end - 1, intron_end, 'splice_donor', mid)
                            self.index.add_feature(t.chrom, intron_start, intron_start + 1, 'splice_acceptor', mid)

            # ---- UTR ----
            for estart, eend in t.exons:
                self.index.add_feature(t.chrom, estart, eend, 'exon_region', mid)

            if t.cds:
                t._cds_min = min(c[0] for c in t.cds)
                t._cds_max = max(c[1] for c in t.cds)
            else:
                t._cds_min = t._cds_max = None

        self.index.finalize()
        logger.info("Index built.")

    def _open_output_files(self):
        files = {
            'cds': 'cds_annotation.txt',
            'intron': 'intron_annotation.txt',
            'utr': 'utr_annotation.txt',
            'stream': 'flank_annotation.txt',
            'intergenic': 'intergenic_annotation.txt'
        }

        headers = {
            'cds': "#CHROM\tPOS\tREF\tALT\tVARIANT_TYPE\tGENE|MRNA\tAA_CHANGE\n",
            'intron': "#CHROM\tPOS\tREF\tALT\tGENE\tMRNA\tFEATURE_ID\n",
            'utr': "#CHROM\tPOS\tREF\tALT\tGENE\tMRNA\tUTR_TYPE\n",
            'stream': "#CHROM\tPOS\tREF\tALT\tGENE\tMRNA\tSTREAM_TYPE\tDISTANCE\n",
            'intergenic': "#CHROM\tPOS\tREF\tALT\tLEFT_GENE\tDIST_LEFT\tRIGHT_GENE\tDIST_RIGHT\n"
        }

        for k, v in files.items():
            f = open(os.path.join(self.output_dir, v), 'w')
            self.handles[k] = f

            if k in headers:
                f.write(headers[k])

        # for k, v in files.items():
        #     self.handles[k] = open(os.path.join(self.output_dir, v), 'w')

    def _close_output_files(self):
        for h in self.handles.values():
            h.close()

    def process_vcf(self):
        logger.info("Streaming VCF and annotating...")
        self._open_output_files()

        def is_gzip(f):
            with open(f, 'rb') as h: return h.read(2) == b'\x1f\x8b'

        opener = gzip.open if is_gzip(self.vcf_file) else open

        t0 = time.time()
        processed = 0

        with opener(self.vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue

                fields = line.strip().split('\t')
                chrom = fields[0]
                try:
                    pos = int(fields[1])
                except:
                    continue

                ref = fields[3].upper()
                alt_field = fields[4].upper()

                if len(ref) != 1 or alt_field == '.':
                    continue

                if ',' in alt_field:
                    alts = alt_field.split(',')
                else:
                    alts = [alt_field]

                for alt_allele in alts:
                    if alt_allele in ('.', '*') or len(alt_allele) != 1:
                        continue
                    self._annotate_one_variant(chrom, pos, ref, alt_allele)

                processed += 1
                if processed % 100000 == 0:
                    sys.stdout.write(f"\rProcessed {processed} variants...")
                    sys.stdout.flush()

        self.total_processed_variants = processed
        sys.stdout.write("\n")
        self._close_output_files()
        logger.info(f"Finished processing {processed} variants in {time.time() - t0:.2f}s")

    def _annotate_one_variant(self, chrom, pos, ref, alt):
        hits = self.index.query(chrom, pos)

        # is_intergenic = True
        hit_types = set()

        cds_hits = []
        exon_hits = []
        intron_hits = []
        stream_hits = []
        splice_donor_mrnas = set()
        splice_acceptor_mrnas = set()

        for ftype, data in hits:
            if ftype == 'CDS':
                cds_hits.append(data)
            elif ftype == 'exon_region':
                exon_hits.append(data)
            elif ftype == 'intron':
                intron_hits.append(data)
            elif ftype in ['upstream', 'downstream']:
                stream_hits.append((ftype, data))
            elif ftype == 'splice_donor':
                splice_donor_mrnas.add(data)
            elif ftype == 'splice_acceptor':
                splice_acceptor_mrnas.add(data)

        # Logic:
        # If CDS hit -> Check NS/Syn. (It is also an exon hit, but CDS is more specific)
        # If Exon hit but NOT CDS hit -> UTR.
        # If Intron hit -> Intron.
        # If Stream hit -> Stream.
        # If NO hit -> Intergenic.

        # CDS Annotation
        if cds_hits:
            # is_intergenic = False
            for mrna_id in cds_hits:
                self._analyze_cds_variant(chrom, pos, ref, alt, mrna_id)
                hit_types.add('cds')

        # UTR Annotation
        if exon_hits:
            # is_intergenic = False
            for mrna_id in exon_hits:
                if mrna_id in cds_hits:
                    continue  # already handled as CDS for this mRNA
                t = self.index.transcripts[mrna_id]
                if t._cds_min is None:
                    continue  # non-coding transcript, skip

                if t.strand == '+':
                    if pos < t._cds_min:
                        tag = '5_prime_UTR_variant'
                    elif pos > t._cds_max:
                        tag = '3_prime_UTR_variant'
                    else:
                        continue
                else:
                    if pos > t._cds_max:
                        tag = '5_prime_UTR_variant'
                    elif pos < t._cds_min:
                        tag = '3_prime_UTR_variant'
                    else:
                        continue

                self.handles['utr'].write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{tag}\t{t.gene_id}\t{mrna_id}\n")
                self.counts['UTR'] += 1
                self.counts['tag'] += 1
                hit_types.add('utr')

        # Intron Annotation
        if intron_hits:
            # is_intergenic = False
            for mid, intron_name in intron_hits:  # data is now (mrna_id, intron_name) tuple
                t = self.index.transcripts.get(mid)
                if t:
                    tag = 'intron_variant'
                    if mid in splice_donor_mrnas:
                        tag = 'splice_donor_variant'
                        self.counts['Splice Donor'] += 1
                    elif mid in splice_acceptor_mrnas:
                        tag = 'splice_acceptor_variant'
                        self.counts['Splice Acceptor'] += 1
                    else:
                        self.counts[tag] += 1
                    self.handles['intron'].write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{tag}\t{t.gene_id}\t{mid}\t{intron_name}\n")
                    self.counts['Intron_Total'] += 1
                    hit_types.add('intron')

        # Stream Annotation
        if stream_hits:
            # is_intergenic = False
            for tag, (gid, mrna_id) in stream_hits:
                t = self.index.transcripts.get(mrna_id)
                if not t: continue
                if t.strand == '+':
                    dist = (t.start - pos) if tag == 'upstream' else (pos - t.end)
                else:
                    dist = (pos - t.end) if tag == 'upstream' else (t.start - pos)
                dist = abs(dist)
                self.handles['stream'].write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{tag}\t{gid}\t{mrna_id}\t{dist}\n")
                self.counts['Upstream/Downstream'] += 1
                self.counts[tag] += 1
                hit_types.add('stream')

        # Intergenic Annotation
        if self._analyze_intergenic(chrom, pos, ref, alt):
            self.counts['Intergenic'] += 1

    def _analyze_cds_variant(self, chrom, pos, ref, alt, mrna_id):
        t = self.index.transcripts[mrna_id]
        if not t.cds_segments or not t.cds_seq:
            return

        try:
            cds_pos = self._genomic_to_cds_pos(t, pos)
            if cds_pos is None:
                return

            frame = cds_pos % 3
            codon_start = cds_pos - frame

            if codon_start + 3 > len(t.cds_seq):
                return

            original_codon = t.cds_seq[codon_start:codon_start + 3].upper()

            if t.strand == '-':
                eff_ref = rev_comp(ref)
                eff_alt = rev_comp(alt)
            else:
                eff_ref = ref
                eff_alt = alt

            if len(eff_ref) != 1 or original_codon[frame] != eff_ref:
                return

            mutated_codon = original_codon[:frame] + eff_alt + original_codon[frame + 1:]

            if len(mutated_codon) != 3:
                return

            aa_ref = CODON2AA.get(original_codon)
            aa_alt = CODON2AA.get(mutated_codon)

            if aa_ref and aa_alt:
                aa_pos = (codon_start // 3) + 1
                vtype = None
                if aa_pos == 1:
                    orig_U = original_codon.upper()
                    mut_U = mutated_codon.upper()
                    if orig_U in START_CODONS and mut_U not in START_CODONS:
                        vtype = 'start_lost'
                    if orig_U in START_CODONS and mut_U in START_CODONS:
                        aa_alt = aa_ref

                if not vtype and aa_ref != aa_alt:
                    vtype = 'missense_variant'
                    if aa_ref == '*' and aa_alt != '*':
                        vtype = 'stop_lost'
                    elif aa_ref != '*' and aa_alt == '*':
                        vtype = 'stop_gained'

                if vtype:
                    change = f"{aa_ref}{aa_pos}{aa_alt}"
                    self.handles['cds'].write(
                        f"{chrom}\t{pos}\t{ref}\t{alt}\t{vtype}\t{t.gene_id}|{mrna_id}\t{change}\n"
                    )
                    self.counts['CDS Non-synonymous'] += 1
                    self.counts[vtype] += 1
                else:
                    self.counts['CDS Synonymous'] += 1

        except Exception as e:
            logger.error(f"Error analyzing CDS variant at {chrom}:{pos} for {mrna_id}: {e}")

    @staticmethod
    def _genomic_to_cds_pos(t, pos):
        import bisect

        segments = t.cds_segments

        if t.strand == '+':
            starts = [s[0] for s in segments]
            idx = bisect.bisect_right(starts, pos) - 1
            if idx < 0:
                return None
            seg_start, seg_end, cum_offset = segments[idx]
            if pos > seg_end:
                return None
            return cum_offset + (pos - seg_start)

        else:
            starts = [s[0] for s in segments]
            ends = [s[1] for s in segments]
            neg_ends = [-e for e in ends]
            idx = bisect.bisect_right(neg_ends, -pos) - 1
            if idx < 0:
                return None
            seg_start, seg_end, cum_offset = segments[idx]
            if pos < seg_start:
                return None
            return cum_offset + (seg_end - pos)

    def _analyze_intergenic(self, chrom, pos, ref, alt):
        if chrom not in self.index.chrom_genes:
            return False

        starts = self.index.chrom_gene_starts[chrom]
        max_ends = self.index.chrom_gene_max_ends[chrom]
        genes = self.index.chrom_genes[chrom]

        import bisect
        idx = bisect.bisect_right(starts, pos)

        left_gene = '.'
        left_dist = '.'
        right_gene = '.'
        right_dist = '.'

        if idx > 0:
            left_max_end, left_max_gid = max_ends[idx - 1]
            if left_max_end >= pos:
                return False

            left_gene = left_max_gid
            left_dist = str(pos - left_max_end)

        if idx < len(genes):
            right_gene = genes[idx][2]
            right_start = genes[idx][0]
            right_dist = str(right_start - pos)

        self.handles['intergenic'].write(
            f"{chrom}\t{pos}\t{ref}\t{alt}\t{left_gene}\t{left_dist}\t{right_gene}\t{right_dist}\n"
        )
        return True

    def generate_report(self):
        import time
        import os
        html_file = os.path.join(self.output_dir, "fast_annotation_summary.html")

        total_snps = self.total_processed_variants

        categories = ['CDS Non-synonymous', 'Intron_Total', 'UTR', 'Upstream/Downstream', 'Intergenic',
                      'CDS Synonymous', 'Splice Donor', 'Splice Acceptor']
        counts = [self.counts.get(c, 0) for c in categories]
        chart_colors = [
            '#ff6384',  # CDS Non-syn (Red)
            '#ff9f40',  # CDS Syn (Orange)
            '#d9534f',  # Splice Donor (Dark Red)
            '#c92a2a',  # Splice Acceptor (Deeper Red)
            '#36a2eb',  # Intron (Blue)
            '#ffce56',  # UTR (Yellow)
            '#4bc0c0',  # Stream (Teal)
            '#9966ff'  # Intergenic (Purple)
        ]

        hierarchy = [
            ("1. CDS Non-synonymous", "CDS Non-synonymous", [
                ("1.1 Missense Variant", "missense_variant"),
                ("1.2 Start Lost", "start_lost"),
                ("1.3 Stop Gained", "stop_gained"),
                ("1.4 Stop Lost", "stop_lost")
            ]),
            ("2. CDS Synonymous", "CDS Synonymous", []),
            ("3. Intronic Region", "Intron_Total", [
                ("3.1 Intron Variant", "intron_variant"),
                ("3.2 Splice Donor", "Splice Donor"),
                ("3.3 Splice Acceptor", "Splice Acceptor")
            ]),
            ("4. UTR Region", "UTR", [
                ("4.1 5' UTR Variant", "5_prime_UTR_variant"),
                ("4.2 3' UTR Variant", "3_prime_UTR_variant")
            ]),
            ("5. Upstream and Downstream", "Upstream/Downstream", [
                ("5.1 Upstream", "upstream"),
                ("5.2 Downstream", "downstream")
            ]),
            ("6. Intergenic Region", "Intergenic", [])
        ]

        table_rows = ""
        for main_label, main_key, sub_list in hierarchy:
            main_count = self.counts.get(main_key, 0)
            main_pct = (main_count / total_snps * 100) if total_snps > 0 else 0

            table_rows += f"""
            <tr class="main-row">
                <td><strong>{main_label}</strong></td>
                <td><strong>{main_count:,}</strong></td>
                <td><strong>{main_pct:.2f}%</strong></td>
            </tr>
            """

            for sub_label, sub_key in sub_list:
                sub_count = self.counts.get(sub_key, 0)
                sub_pct = (sub_count / total_snps * 100) if total_snps > 0 else 0
                rel_pct = (sub_count / main_count * 100) if main_count > 0 else 0

                if sub_count > 0:
                    table_rows += f"""
                    <tr class="sub-row">
                        <td class="indent">{sub_label}</td>
                        <td>{sub_count:,}</td>
                        <td>{sub_pct:.3f}% <span class="rel-pct">( {rel_pct:.1f}% of region )</span></td>
                    </tr>
                    """

        html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>EvoAnn's Annotation Report</title>
        <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
        <link href="https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap" rel="stylesheet">
        <style>
            body {{ font-family: 'Roboto', sans-serif; background-color: #f4f6f8; color: #333; margin: 0; padding: 0; }}
            .header {{ background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%); color: white; padding: 40px 20px; text-align: center; }}
            .header h1 {{ margin: 0; font-size: 2.5em; }}
            .header p {{ opacity: 0.8; margin-top: 10px; }}

            .container {{ max-width: 1200px; margin: -30px auto 40px; padding: 0 20px; }}

            .card {{ background: white; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); padding: 25px; margin-bottom: 25px; display: flex; flex-direction: column; }}
            .card h2 {{ border-bottom: 2px solid #f0f0f0; padding-bottom: 10px; margin-top: 0; color: #2c3e50; font-size: 1.2em; flex-shrink: 0; }}

            .summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; text-align: center; }}
            .stat-box {{ padding: 15px; background: #f8f9fa; border-radius: 8px; border: 1px solid #e9ecef; }}
            .stat-value {{ font-size: 2em; font-weight: 700; color: #2a5298; }}
            .stat-label {{ color: #6c757d; font-size: 0.9em; text-transform: uppercase; letter-spacing: 1px; margin-top: 5px; }}

            .charts-row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; align-items: stretch; }}
            @media (max-width: 900px) {{ .charts-row {{ grid-template-columns: 1fr; }} }}

            .chart-wrapper {{ position: relative; height: 400px; width: 100%; flex-grow: 1; margin-top: 10px; }}
            .table-responsive {{ height: 400px; overflow-y: auto; flex-grow: 1; margin-top: 10px; }}

            table {{ width: 100%; border-collapse: collapse; position: relative; }}
            th, td {{ padding: 12px 15px; text-align: left; border-bottom: 1px solid #e1e4e8; }}
            th {{ background-color: #f1f3f5; font-weight: 600; color: #495057; position: sticky; top: 0; z-index: 10; box-shadow: 0 1px 2px rgba(0,0,0,0.05); }}
            tr:last-child td {{ border-bottom: none; }}
            table th:nth-child(2), table td:nth-child(2),
            table th:nth-child(3), table td:nth-child(3) {{ white-space: nowrap; }}
            table th:nth-child(1), table td:nth-child(1) {{ min-width: 200px; }}

            .main-row {{ background-color: #ffffff; transition: background 0.2s; }}
            .main-row:hover {{ background-color: #f8f9fa; }}
            .sub-row {{ background-color: #fafbfc; }}
            .sub-row:hover {{ background-color: #f1f3f5; }}
            .indent {{ padding-left: 45px; color: #495057; font-size: 0.95em; position: relative; }}
            .indent::before {{ content: ""; position: absolute; left: 25px; top: 50%; width: 5px; height: 5px; background-color: #adb5bd; border-radius: 50%; transform: translateY(-50%); }}
            .rel-pct {{ font-size: 0.85em; color: #868e96; margin-left: 8px; font-weight: normal; }}

            ::-webkit-scrollbar {{ width: 8px; }}
            ::-webkit-scrollbar-track {{ background: #f8f9fa; }}
            ::-webkit-scrollbar-thumb {{ background: #c1c1c1; border-radius: 4px; }}
            ::-webkit-scrollbar-thumb:hover {{ background: #a8a8a8; }}

            .footer {{ text-align: center; color: #adb5bd; margin-top: 20px; padding-bottom: 40px; }}
            a {{ color: #2a5298; text-decoration: none; font-weight: 500; }}
            a:hover {{ text-decoration: underline; }}
        </style>
    </head>
    <body>

        <div class="header">
            <h1>Annotation Report</h1>
            <p>Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>

        <div class="container">

            <div class="card">
                <h2>Overview</h2>
                <div class="summary-grid">
                    <div class="stat-box">
                        <div class="stat-value">{total_snps:,}</div>
                        <div class="stat-label">Total Variants</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{self.counts.get('CDS Non-synonymous', 0):,}</div>
                        <div class="stat-label">Non-synonymous</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{len(self.index.chrom_genes):,}</div>
                        <div class="stat-label">Chromosomes</div>
                    </div>
                </div>
            </div>

            <div class="charts-row">
                <div class="card">
                    <h2>Distribution by Region</h2>
                    <div class="chart-wrapper">
                        <canvas id="regionChart"></canvas>
                    </div>
                </div>
                <div class="card">
                    <h2>Number of annotaitons and region counts</h2>
                    <div class="table-responsive">
                        <table>
                            <thead>
                                <tr><th>Region / Variant Type</th><th>Count</th><th>Percentage</th></tr>
                            </thead>
                            <tbody>
                                {table_rows}
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>

            <div class="card" style="display: block;">
                <h2>Output Files</h2>
                <ul style="list-style-type: none; padding: 0;">
                    <li style="margin-bottom: 10px;">📄 <a href="cds_annotation.txt" download>CDS Annotation (Non-synonymous)</a> - <span style="color:#666">Detailed list of protein-altering variants</span></li>
                    <li style="margin-bottom: 10px;">📄 <a href="intron_annotation.txt" download>Intron Annotation</a></li>
                    <li style="margin-bottom: 10px;">📄 <a href="utr_annotation.txt" download>UTR Annotation</a></li>
                    <li style="margin-bottom: 10px;">📄 <a href="flank_annotation.txt" download>Upstream/Downstream Annotation</a></li>
                    <li style="margin-bottom: 10px;">📄 <a href="intergenic_annotation.txt" download>Intergenic Annotation</a></li>
                </ul>
            </div>

        </div>

        <div class="footer">
            Powered by <strong>Evoann</strong> | An ultra-fast variant annotation toolkit
        </div>

        <script>
            const ctx = document.getElementById('regionChart').getContext('2d');
            new Chart(ctx, {{
                type: 'doughnut',
                data: {{
                    labels: {categories},
                    datasets: [{{
                        data: {counts},
                        backgroundColor: {chart_colors},
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            position: 'right',
                            labels: {{ boxWidth: 12, font: {{size: 11}} }}
                        }}
                    }},
                    layout: {{ padding: 10 }}
                }}
            }});
        </script>
    </body>
    </html>
    """

        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        logger.info(f"Report generated: {html_file}")


def run(args):
    ann = FastAnnotator(args.vcf, args.gff, args.fasta, args.output, args.flank)
    ann.load_genome()
    ann.parse_gff()
    ann.process_vcf()
    ann.generate_report()
