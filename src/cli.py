import argparse
from . import annotate_CDS
from . import deal_cds
from . import annotate_Intron
from . import annotate_UTR
from . import annotate_Stream
from . import deal_consensus
from . import GFFCorrector
from . import annotate_Intergenic
from . import deal_phased
from . import fast_annotate


class HelpFormatter(argparse.HelpFormatter):
    def __init__(self, prog):
        super().__init__(prog, max_help_position=30, width=100)

    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)

        parts = []
        if action.option_strings:
            parts.extend(action.option_strings)
        return ', '.join(parts)


def setup_cds_ann(subparsers):
    p = subparsers.add_parser('cds-ann', help='Standalone: Annotate non-synonymous variants in CDS.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-v', '--vcf', required=True, metavar='FILE', help='Input VCF.')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF.')
    req.add_argument('-c', '--cds', required=True, metavar='FILE', help='CDS FASTA file.')
    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-o', '--output', default='cds_ann_out', metavar='FILE', help='Output file(default:cds_ann_out).')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=annotate_CDS.run)


def setup_intron_ann(subparsers):
    p = subparsers.add_parser('intron-ann', help='Standalone: Annotate non-synonymous variants in introns.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-v', '--vcf', required=True, metavar='FILE', help='Input VCF.')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF.')
    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-o', '--output', default='intron_out', metavar='FILE', help='Output file(default:intron_out).')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=annotate_Intron.run)


def setup_utr_ann(subparsers):
    p = subparsers.add_parser('utr-ann', help='Standalone: Annotate non-synonymous variants in UTRs.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-v', '--vcf', required=True, metavar='FILE', help='Input VCF.')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF.')
    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-all', '--all_transcripts', action='store_true',
                     help='Extract all transcripts (default: longest transcript only).')
    opt.add_argument('-o', '--output', default='utr_out', metavar='FILE', help='Output file(default:utr_out).')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=annotate_UTR.run)


def setup_stream_ann(subparsers):
    p = subparsers.add_parser('flank-ann', help='Standalone: Annotate upstream/downstream variants.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-v', '--vcf', required=True, metavar='FILE', help='Input VCF.')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF.')
    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-o', '--output', default='stream_out', metavar='FILE', help='Output file(default:flank_out).')
    opt.add_argument('--flank', type=int, default=2000, metavar='INT',
                     help='Upstream/downstream range (bp), default: 2000.')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=annotate_Stream.run)


def setup_cds_extract(subparsers):
    p = subparsers.add_parser('cds-extract', help='Standalone: Extract CDS sequences from genome.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF.')
    req.add_argument('-f', '--fasta', required=True, metavar='FILE', help='Input FASTA.')
    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-all', '--all_transcripts', action='store_true',
                     help='Extract all transcripts (default: longest transcript only).')
    opt.add_argument('-o', '--output', default='cds_out', metavar='FILE', help='Output file(default:cds_out).')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=deal_cds.run)


def setup_consensus_extract(subparsers):
    p = subparsers.add_parser('consensus-extract',
                              help='Standalone: Extract consensus sequences for genes from multi-sample VCF.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-v', '--vcf', required=True, metavar='FILE', help='Input VCF file (supports .vcf or .vcf.gz).')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF annotation file.')
    req.add_argument('-f', '--fasta', required=True, metavar='FILE', help='Input FASTA.')

    gene = p.add_argument_group('Gene Selection')
    gene.add_argument('--gene', metavar='NAME', help='Single gene name to extract.')
    gene.add_argument('--gene-file', metavar='FILE', help='File containing gene names (one per line).')
    gene.add_argument('--cds-only', action='store_true',
                      help='Extract ONLY the coding sequence (CDS) and reverse complement negative strands.')
    gene.add_argument('--all-transcripts', action='store_true',
                      help='Extract CDS for ALL transcripts instead of just the longest one (Requires --cds-only).')

    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-o', '--output-dir', default='consensus_output', metavar='DIR',
                     help='Output directory (default: consensus_output).')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=deal_consensus.run)


def setup_gff_correct(subparsers):
    p = subparsers.add_parser('correct-gff',
                              help='Standalone: Fix GFF files by adding missing exons and correcting CDS phases.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF.')

    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-o', '--output', metavar='FILE',
                     help='Output GFF file (default: input_fixed.gff).')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=GFFCorrector.run)


def setup_intergenic_ann(subparsers):
    p = subparsers.add_parser('intergenic-ann', help='Standalone: Annotate intergenic variants.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-v', '--vcf', required=True, metavar='FILE', help='Input VCF.')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF.')
    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-o', '--output', default='intergenic_out.txt', metavar='FILE',
                     help='Output file(default:intergenic_out.txt).')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=annotate_Intergenic.run)


def setup_phased_extract(subparsers):
    p = subparsers.add_parser('phased-extract',
                              help='Standalone: Extract phased haplotype sequences for genes.',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-v', '--vcf', required=True, metavar='FILE', help='Input VCF file (supports .vcf or .vcf.gz).')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF annotation file.')
    req.add_argument('-f', '--fasta', required=True, metavar='FILE', help='Input FASTA.')

    gene = p.add_argument_group('Gene Selection (at least one required)')
    gene.add_argument('--gene', metavar='NAME', help='Single gene name to extract.')
    gene.add_argument('--gene-file', metavar='FILE', help='File containing gene names (one per line).')
    gene.add_argument('--cds-only', action='store_true',
                      help='Extract ONLY the coding sequence (CDS) and reverse complement negative strands.')
    gene.add_argument('--all-transcripts', action='store_true',
                      help='Extract CDS for ALL transcripts instead of just the longest one (Requires --cds-only).')

    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('-o', '--output-dir', default='phased_output', metavar='DIR',
                     help='Output directory (default: phased_output).')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=deal_phased.run)


def setup_fast_ann(subparsers):
    p = subparsers.add_parser('all-ann',
                              help='Optimized one-pass annotation pipeline (Faster, Less Memory).',
                              formatter_class=HelpFormatter, add_help=False)
    req = p.add_argument_group('Required Arguments')
    req.add_argument('-v', '--vcf', required=True, metavar='FILE', help='Input VCF file.')
    req.add_argument('-g', '--gff', required=True, metavar='FILE', help='Input GFF file.')
    req.add_argument('-f', '--fasta', required=True, metavar='FILE', help='Input Genome FASTA file.')
    req.add_argument('-o', '--output', required=True, metavar='DIR', help='Output directory.')

    opt = p.add_argument_group('Optional Arguments')
    opt.add_argument('--flank', type=int, default=2000, metavar='INT',
                     help='Upstream/downstream range (bp), default: 2000.')
    opt.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
    p.set_defaults(func=fast_annotate.run)


def create_parser():
    menu_text = """
        EvoAnn: Evolutionary Annotation Toolkit v1.0
        ----------------------------------------------------
        A high-performance, GFF-based variant annotation and sequence 
        processing suite tailored for evolutionary genomics.
    
        evoann <command> [options]
        example: evoann fast-ann -v <VCF> -g <GFF> -f <Reference> -o <OutFile> [args...]
        or evoann cds-ann -v <VCF> -g <GFF> -c <CDS> -o <OutFile>

        Commands:
          cds-ann               Annotate CDS variants.
          intron-ann            Annotate intron variants.
          utr-ann               Annotate UTR variants.
          stream-ann            Annotate up/downstream variants.
          intergenic-ann        Annotate intergenic variants and nearest genes.
          cds-extract           Extract CDS sequences from genome based on GFF models.
          consensus-extract     Generate consensus sequences for specific samples/genes (handles Indels).     
          phased-extract        Extract haplotype sequences (Hap1/Hap2) from phased VCFs.
          correct-gff           Repair GFF: fix missing exons and correct CDS phases.
          all-ann              High-speed, one-pass annotation for CDS/Intron/UTR/Intergenic.
        """

    parser = argparse.ArgumentParser(
        prog='evoann',
        usage='evoann <command> [options]\n',
        description=menu_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='For more details on a specific command, run: ann <command> --help',
        add_help=False
    )

    subparsers = parser.add_subparsers(dest='command', metavar='COMMAND', help=argparse.SUPPRESS)

    # common = parser.add_argument_group('Common Arguments (used in subcommands)')
    # common.add_argument('-v, --vcf', help='Input VCF file (bgzipped & indexed).', action='store_true')
    # common.add_argument('-g, --gff', help='Reference GFF3 annotation file.', action='store_true')
    # common.add_argument('-o, --output', help='Output directory or file.', action='store_true')
    # common.add_argument('-f', '--fasta', help='Reference genome FASTA file.', action='store_true')
    # common.add_argument('-c', '--cds', help='Input CDS FILE.', action='store_true')
    # common.add_argument('--flank', help='Upstream/downstream range (bp), default: 2000.', action='store_true')

    setup_cds_ann(subparsers)
    setup_intron_ann(subparsers)
    setup_utr_ann(subparsers)
    setup_stream_ann(subparsers)
    setup_cds_extract(subparsers)
    setup_consensus_extract(subparsers)
    setup_gff_correct(subparsers)
    setup_intergenic_ann(subparsers)
    setup_phased_extract(subparsers)
    setup_fast_ann(subparsers)

    return parser
