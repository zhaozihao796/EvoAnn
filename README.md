================================================================================
EvoAnn: Evolutionary Annotation Toolkit
================================================================================

Welcome to EvoAnn! This directory contains the source code and essential 
resources for running the software. EvoAnn operates as a unified command-line 
tool, seamlessly integrated into your environment.

--------------------------------------------------------------------------------
1. DIRECTORY STRUCTURE
--------------------------------------------------------------------------------
├── environment.yaml      # Conda environment dependencies file
├── setup.py              # Python package installation script
├── src/                  # Source code directory
├── README.md             # This manual file
├── LICENSE               # Terms of use (MIT License)
└── test_data/            # Small empirical dataset for quick verification
    ├── test_ref.fa
    ├── test_anno.gff3
    └── test_pop.vcf.gz

--------------------------------------------------------------------------------
2. INSTALLATION
--------------------------------------------------------------------------------
EvoAnn manages its dependencies robustly via Conda and is installed as a 
standard Python package.

Step 1: Create the isolated environment using the provided yaml file
  conda env create -f environment.yaml

Step 2: Activate the newly created environment (assuming the name is 'evoann')
  conda activate evoann

Step 3: Install the package locally
  pip install -e .

Once installed, the `evoann` command will be globally available in this Conda 
environment.

--------------------------------------------------------------------------------
3. QUICK TEST RUN (Sanity Check)
--------------------------------------------------------------------------------
Before processing your own data, we highly recommend running the provided 
test dataset to ensure EvoAnn is functioning correctly on your system.

Ensure your Conda environment is activated, then run:
  
  evoann all-ann -f test_data/test_ref.fa -g test_data/test_anno.gff3 -v test_data/test_pop.vcf.gz -o test_out

If successful, you will see a `test_out/` directory containing the type-stratified 
results and the interactive HTML summary.

--------------------------------------------------------------------------------
4. QUICK START & EXAMPLES
--------------------------------------------------------------------------------

EvoAnn operates using a unified CLI with specific subcommands. Below are the 
most common use cases.

### 1. One-pass Global Annotation (Faster, Less Memory)
To perform an exhaustive variant annotation across all genomic regions 
(CDS, Intron, UTR, Intergenic, etc.) in a single run:

  evoann all-ann \
    -v population.vcf.gz \
    -g annotation.gff3 \
    -f reference.fa \
    -o ./evoann_results

### 2. Extract Phased Haplotypes for Specific Genes
To extract Hap1/Hap2 sequences utilizing IUPAC ambiguity encoding:

  evoann phased-extract \
    -v population.vcf.gz \
    -g annotation.gff3 \
    -f reference.fa \
    --gene-file target_genes.txt \
    --cds-only \
    -o ./phased_output

### 3. Repair Fragmented GFF3 Files
To automatically fix missing exons and correct CDS phases before downstream analysis:

  evoann correct-gff \
    -g raw_annotation.gff3 \
    -o fixed_annotation.gff3


--------------------------------------------------------------------------------
🛠️ MAIN SUBCOMMANDS & USAGE
--------------------------------------------------------------------------------

You can view the full help menu by running `evoann --help`. 
For detailed options on a specific command, run `evoann <command> --help`.

[A] COMPREHENSIVE PIPELINE
-------------------------------------------
* all-ann : Optimized one-pass annotation pipeline for CDS, Introns, UTRs, and Intergenic regions.
  Usage: evoann all-ann -v <VCF> -g <GFF> -f <FASTA> -o <OUTDIR> [--flank INT]

[B] SEQUENCE EXTRACTION TOOLS
-------------------------------------------
* consensus-extract : Generate consensus sequences for specific samples/genes (handles Indels).
  Usage: evoann consensus-extract -v <VCF> -g <GFF> -f <FASTA> [--gene NAME | --gene-file FILE] [--cds-only] -o <OUTDIR>

* phased-extract    : Extract haplotype sequences (Hap1/Hap2) from phased VCFs.
  Usage: evoann phased-extract -v <VCF> -g <GFF> -f <FASTA> [--gene NAME | --gene-file FILE] [--cds-only] -o <OUTDIR>

* cds-extract       : Extract CDS sequences directly from the reference genome based on GFF models.
  Usage: evoann cds-extract -g <GFF> -f <FASTA> [-all] -o <OUTFILE>

[C] STANDALONE ANNOTATION MODULES
-------------------------------------------
Use these if you only need to process specific genomic features:
* cds-ann        : Annotate non-synonymous variants in CDS. (Requires -c <CDS_FASTA>)
* intron-ann     : Annotate variants in introns.
* utr-ann        : Annotate variants in UTRs. (Supports -all transcripts)
* stream-ann     : Annotate up/downstream variants. (Supports --flank range)
* intergenic-ann : Annotate intergenic variants and nearest genes.

[D] UTILITIES
-------------------------------------------
* correct-gff    : Repair GFF files by fixing missing exons and correcting CDS phases.