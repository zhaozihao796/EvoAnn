import gzip
import time
import sys
from Bio import SeqIO

# 工具函数

def parse_attributes(attr_field):
    attrs = {}
    for attr in attr_field.strip().split(';'):
        if not attr:
            continue
        attr = attr.strip()
        if '=' in attr:
            key, value = attr.split('=')
            attrs[key] = value
        else:
            parts = attr.split()
            if len(parts) >= 2:
                attrs[parts[0].strip()] = parts[1].strip('"')
    return attrs


def fasta_to_dict(fasta_file):
    fasta_dict = {}
    fasta_record = SeqIO.parse(fasta_file, 'fasta')
    for record in fasta_record:
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict

def parse_vcf_all2(vcf_file):
    vcf_dict = {}

    def is_gzip_file(filename):
        try:
            with open(filename, 'rb') as f:
                return f.read(2) == b'\x1f\x8b'
        except:
            return False

    opener = gzip.open if is_gzip_file(vcf_file) else open
    mode = 'rt' if is_gzip_file(vcf_file) else 'r'

    kept = 0
    t0 = time.time()
    
    try:
        f_obj = opener(vcf_file, mode)
    except Exception as e:
        raise IOError(f"Failed to open VCF file {vcf_file}: {e}")

    with f_obj as f:
        for i, line in enumerate(f, 1):
            if i % 100000 == 0:
                elapsed = time.time() - t0
                speed = i / elapsed if elapsed > 0 else 0
                sys.stdout.write(f"\r[VCF] Line {i} | Collected {kept} SNPs | Speed {speed:.0f} lines/s | Time elapsed {elapsed:.1f}s   ")
                sys.stdout.flush()
            if not line or line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 7:
                continue
            chrom = fields[0]
            try:
                pos = int(fields[1])
            except:
                continue
            ref = fields[3].upper()
            alt_field = fields[4].upper()
            filt = fields[6]

            if len(ref) != 1 or alt_field == '.':
                continue
            if ',' in alt_field:
                for a in alt_field.split(','):
                    if a == '*':
                        continue
                    if len(a) == 1:
                        vcf_dict.setdefault(chrom, []).append((pos, ref, a))
                        kept += 1
            else:
                if alt_field == '*':
                    continue
                if len(alt_field) == 1:
                    vcf_dict.setdefault(chrom, []).append((pos, ref, alt_field))
                    kept += 1

    for chrom in vcf_dict:
        vcf_dict[chrom].sort(key=lambda x: x[0])

    elapsed = time.time() - t0
    sys.stdout.write("\r[VCF] Parsing completed | Total SNPs: {} | Total time: {:.1f}s                    \n".format(kept, elapsed))
    sys.stdout.flush()
    return vcf_dict

def rev_comp(seq):
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def _check_cds(seq: str):
    seq_u = (seq or '').upper()
    stop_set = {'TAA', 'TAG', 'TGA'}
    has_atg = seq_u.startswith('ATG') if len(seq_u) >= 3 else False
    has_stop = (len(seq_u) >= 3 and seq_u[-3:] in stop_set)
    in_frame = (len(seq_u) % 3 == 0)
    return has_atg, has_stop, in_frame


def parse_cds_fasta(cds_fasta_file: str, validate: bool = True):
    cds_dict = {}
    report = {} if validate else None
    count = 0
    t0 = time.time()

    for record in SeqIO.parse(cds_fasta_file, 'fasta'):
        gid = record.id
        seq = str(record.seq).upper()
        cds_dict[gid] = seq
        count += 1
        if validate:
            has_atg, has_stop, in_frame = _check_cds(seq)
            report[gid] = {
                'has_atg': has_atg,
                'has_stop': has_stop,
                'in_frame': in_frame,
                'length': len(seq),
            }

    print(f"[CDS-FASTA] Parsing completed: {count} records")
    return cds_dict, report


def repair_gff_missing_exons_enhanced(input_gff, output_gff=None):
    if output_gff is None:
        output_gff = input_gff.replace('.gff', '_repaired.gff')

    genes = {}
    mrnas = {}
    existing_exons = {}  # mRNA ID -> [exon intervals]
    cds_by_mrna = {}  # mRNA ID -> [CDS intervals]
    other_features = []

    print("Parsing GFF file...")

    with open(input_gff, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if not line.endswith('\n'):
                line = line + '\n'

            if line.startswith('#'):
                other_features.append(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                other_features.append(line)
                continue

            chrom, source, feature, start, end, score, strand, phase, attr_field = fields
            attrs = parse_attributes(attr_field)

            if feature == 'gene':
                gene_id = attrs.get('ID')
                if gene_id:
                    genes[gene_id] = line
            elif feature == 'mRNA':
                mrna_id = attrs.get('ID')
                if mrna_id:
                    mrnas[mrna_id] = line
            elif feature == 'exon':
                parent = attrs.get('Parent')
                if parent:
                    start_int, end_int = int(start), int(end)
                    existing_exons.setdefault(parent, []).append((start_int, end_int, line))
            elif feature == 'CDS':
                parent = attrs.get('Parent')
                if parent:
                    start_int, end_int = int(start), int(end)
                    cds_by_mrna.setdefault(parent, []).append((start_int, end_int, line))
            else:
                other_features.append(line)

    print(f"Parsing completed: {len(genes)} genes, {len(mrnas)} mRNA")

    with open(output_gff, 'w') as f_out:
        for line in other_features:
            f_out.write(line)

        for gene_line in genes.values():
            f_out.write(gene_line)

        for mrna_id, mrna_line in mrnas.items():
            f_out.write(mrna_line)

            current_exons = existing_exons.get(mrna_id, [])
            cds_list = cds_by_mrna.get(mrna_id, [])

            for _, _, exon_line in current_exons:
                f_out.write(exon_line)

            i = 1

            for start, end, cds_line in cds_list:
                if not cds_line.endswith('\n'):
                    cds_line = cds_line + '\n'
                f_out.write(cds_line)
                pi = False
                for exon_start, exon_end, _ in current_exons:
                    if start >= exon_start and end <= exon_end:
                        pi = True
                        i += 1
                        break
                if not pi:
                    f_out.write(create_exon_line(mrna_line, start, end, mrna_id, i))
                    i += 1
    print(f"\nOutput file: {output_gff}")
    return output_gff


def create_exon_line(mrna_line, exon_start, exon_end, mrna_id, exon_num):
    mrna_fields = mrna_line.strip().split('\t')
    chrom, source = mrna_fields[0], mrna_fields[1]
    strand = mrna_fields[6]

    exon_id = f"{mrna_id}.exon_{exon_num}"
    attr_field = f"ID={exon_id};Parent={mrna_id}"

    return f"{chrom}\t{source}\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{attr_field}\n"
