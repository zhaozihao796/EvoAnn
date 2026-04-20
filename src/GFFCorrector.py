import logging
from typing import Dict, List, Tuple
from .annotate_utils import parse_attributes

logger = logging.getLogger(__name__)


class GFFFeature:

    def __init__(self, line: str, start: int, end: int, strand: str, feature_type: str, attributes: dict):
        self.line = line.strip()
        self.start = start
        self.end = end
        self.strand = strand
        self.feature_type = feature_type
        self.attributes = attributes
        self.children = []

    def update_phase(self, new_phase: int):
        fields = self.line.split('\t')
        if len(fields) >= 8:
            fields[7] = str(new_phase)
            self.line = '\t'.join(fields)

    def update_id(self, new_id: str):
        self.attributes['ID'] = new_id

        new_attr_parts = []

        new_attr_parts.append(f"ID={new_id}")

        for k, v in self.attributes.items():
            if k == 'ID': continue
            new_attr_parts.append(f"{k}={v}")

        new_attr_str = ";".join(new_attr_parts)

        fields = self.line.split('\t')
        if len(fields) >= 9:
            fields[8] = new_attr_str
            self.line = '\t'.join(fields)


class GFFCorrector:
    def __init__(self, input_gff: str, output_gff: str):
        self.input_gff = input_gff
        self.output_gff = output_gff

        self.headers = []
        self.genes: Dict[str, GFFFeature] = {}  # gene_id -> GeneFeature
        self.mrnas: Dict[str, GFFFeature] = {}  # mrna_id -> mRNAFeature

        self.mrna_children: Dict[str, Dict[str, List[GFFFeature]]] = {}

    def run_pipeline(self):
        logger.info(f"Initiating GFF file repair: {self.input_gff}")

        self._load_gff()

        total_repaired_exons = 0
        total_corrected_phases = 0

        for mrna_id, children in self.mrna_children.items():
            repaired_count = self._repair_exons(mrna_id, children)
            total_repaired_exons += repaired_count

            corrected_count = self._correct_phases(children['CDS'])
            total_corrected_phases += corrected_count

        logger.info(
            f"Repair completed: {total_repaired_exons} Exons supplemented, {total_corrected_phases} Phases corrected")

        self._write_gff()
        logger.info(f"Results saved to: {self.output_gff}")

    def _load_gff(self):
        logger.info("Parsing GFF structure...")
        count = 0
        with open(self.input_gff, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith('#'):
                    self.headers.append(line)
                    continue

                fields = line.split('\t')
                if len(fields) < 9: continue

                chrom, _, feature, start, end, _, strand, phase, attr_str = fields
                attrs = parse_attributes(attr_str)
                obj = GFFFeature(line, int(start), int(end), strand, feature, attrs)

                if feature == 'gene':
                    gid = attrs.get('ID')
                    if gid: self.genes[gid] = obj

                elif feature == 'mRNA':
                    mid = attrs.get('ID')
                    parent = attrs.get('Parent')
                    if mid:
                        self.mrnas[mid] = obj
                        if parent and parent in self.genes:
                            self.genes[parent].children.append(mid)

                else:
                    parent = attrs.get('Parent')
                    if parent:
                        if parent not in self.mrna_children:
                            self.mrna_children[parent] = {'exon': [], 'CDS': [], 'other': []}

                        if feature in ['exon', 'CDS']:
                            self.mrna_children[parent][feature].append(obj)
                        else:
                            self.mrna_children[parent]['other'].append(obj)

                count += 1
        logger.info(f"Parsing completed, {count} records processed in total")

    def _repair_exons(self, mrna_id, children_map) -> int:
        exons = children_map['exon']
        cdss = children_map['CDS']
        others = children_map.get('other', [])

        if not cdss: return 0
        if mrna_id not in self.mrnas: return 0
        mrna_obj = self.mrnas[mrna_id]

        new_exons_count = 0

        if not exons:
            utrs = [o for o in others if 'UTR' in o.feature_type.upper()]

            parts = [(c.start, c.end) for c in cdss] + [(u.start, u.end) for u in utrs]
            if not parts: return 0

            parts.sort(key=lambda x: x[0])
            merged_exons = []

            for s, e in parts:
                if not merged_exons or s > merged_exons[-1][1] + 1:
                    merged_exons.append([s, e])
                else:
                    merged_exons[-1][1] = max(merged_exons[-1][1], e)

            for s, e in merged_exons:
                fields = mrna_obj.line.split('\t')
                fields[2] = 'exon'
                fields[3] = str(s)
                fields[4] = str(e)
                fields[7] = '.'
                fields[8] = f"ID={mrna_id}_temp_exon;Parent={mrna_id}"
                new_line = '\t'.join(fields)

                new_exon = GFFFeature(new_line, s, e, mrna_obj.strand, 'exon', {})
                new_exon.attributes['Parent'] = mrna_id
                children_map['exon'].append(new_exon)
                new_exons_count += 1

            return new_exons_count

        for cds in cdss:
            is_covered = False
            for exon in exons:
                if cds.start >= exon.start and cds.end <= exon.end:
                    is_covered = True
                    break

            if not is_covered:
                fields = mrna_obj.line.split('\t')
                fields[2] = 'exon'
                fields[3] = str(cds.start)
                fields[4] = str(cds.end)
                fields[7] = '.'
                fields[8] = f"ID={mrna_id}_temp_exon;Parent={mrna_id}"
                new_line = '\t'.join(fields)

                new_exon = GFFFeature(new_line, cds.start, cds.end, cds.strand, 'exon', {})
                new_exon.attributes['Parent'] = mrna_id
                children_map['exon'].append(new_exon)
                new_exons_count += 1

        return new_exons_count

    def _correct_phases(self, cds_list: List[GFFFeature]) -> int:
        if not cds_list: return 0

        strand = cds_list[0].strand
        if strand == '+':
            sorted_cds = sorted(cds_list, key=lambda x: x.start)
        else:
            sorted_cds = sorted(cds_list, key=lambda x: x.start, reverse=True)

        corrections = 0
        prev_length = 0
        prev_phase = 0

        for i, cds in enumerate(sorted_cds):
            fields = cds.line.split('\t')
            try:
                orig_phase = int(fields[7]) if fields[7] != '.' else 0
            except:
                orig_phase = 0

            if i == 0:
                expected = 0
            else:
                expected = (3 - ((prev_length - prev_phase) % 3)) % 3

            if orig_phase != expected:
                cds.update_phase(expected)
                corrections += 1
                prev_phase = expected
            else:
                prev_phase = orig_phase

            prev_length = cds.end - cds.start + 1

        return corrections

    def _write_gff(self):
        with open(self.output_gff, 'w') as f:
            for h in self.headers:
                f.write(h + '\n')

            sorted_genes = sorted(self.genes.values(), key=lambda x: x.start)

            for gene in sorted_genes:
                f.write(gene.line + '\n')

                for mrna_id in gene.children:
                    if mrna_id in self.mrnas:
                        mrna = self.mrnas[mrna_id]
                        f.write(mrna.line + '\n')

                        if mrna_id in self.mrna_children:
                            children = self.mrna_children[mrna_id]
                            exons = sorted(children['exon'], key=lambda x: x.start)
                            for i, ex in enumerate(exons, 1):
                                new_id = f"{mrna_id}.exon{i}"
                                ex.update_id(new_id)
                                f.write(ex.line + '\n')

                            cdss = sorted(children['CDS'], key=lambda x: x.start)
                            for i, cds in enumerate(cdss, 1):
                                new_id = f"{mrna_id}.cds{i}"
                                cds.update_id(new_id)
                                f.write(cds.line + '\n')

                            others = sorted(self.mrna_children[mrna_id].get('other', []), key=lambda x: x.start)
                            utr_counters = {}
                            for o in others:
                                ftype = o.feature_type
                                utr_counters[ftype] = utr_counters.get(ftype, 0) + 1
                                new_id = f"{mrna_id}.{ftype}.{utr_counters[ftype]}"
                                o.update_id(new_id)
                                f.write(o.line + '\n')


def repair_gff_missing_exons_and_phases(input_gff, output_gff=None):
    if output_gff is None:
        output_gff = input_gff.replace('.gff', '_fixed.gff')

    corrector = GFFCorrector(input_gff, output_gff)
    corrector.run_pipeline()
    return output_gff


def run(args):
    input_gff = args.gff
    output_gff = args.output if args.output else input_gff.replace('.gff', '_fixed.gff')

    repair_gff_missing_exons_and_phases(input_gff, output_gff)
