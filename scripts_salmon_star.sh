###################################
# select GEO fastq files from ENA with download script on database 

# Script for iterating fastqc over samples. Need to ensure that load module fastqc and parallel. 

# inputs directory where the fastq files are 

#!/bin/bash 

if [ ! -d "$1" ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

mydir="${1}"
dirname=$(basename "$mydir")

# Create an output directory for FastQC reports in the working directory 
output_dir="$(pwd)/${dirname}_fastqc"
mkdir -p "$output_dir"

# Loop over all files matching the pattern *.fastq.* and run fastqc in parallel

ls ${mydir}/*.fastq.* | parallel -j 12 '
  file="{}"  # This is the file passed to parallel
  echo "Processing file: $file"

  # Extract the filename and sample name
  filename=$(basename "$file")
  samplename="${filename%%.*}"

  echo "Current file:               $file"
  echo "Basename of file:           $filename"
  echo "Sample:                     $samplename"

  # Run FastQC
  fastqc "$file" -o "$output_dir"
  echo -e "##########################\n\n"
'
# merge runs of same biological samples and add experiment accession number with below
# - replace with file names

cat RUN_1_1.fastq RUN_2_1.fastq RUN_3_1.fastq > Experiment_Accession_1.fastq

#############################################################################

# script for carrying out trimming with fastp. Ensure to load fastp module. 

# used input as fastq files 

# need to specify the output directory want as well. Specified when running on samples e.g. bulk_ACP_trimmed 

#!/bin/bash

input_dir="$1"
output_dir="$2"
mkdir -p "$output_dir"

# Loop through all SR*_1.fastq.gz files
for r1 in "$input_dir"/SR*_1.fastq.gz; do
    sample=$(basename "$r1" _1.fastq.gz)
    r2="${input_dir}/${sample}_2.fastq.gz"

    if [ ! -f "$r2" ]; then
        echo "Warning: Missing paired file for $sample. Skipping..."
        continue
    fi

    ./bin/fastp \
        -i "$r1" \
        -I "$r2" \
        -o "${output_dir}/${sample}_1.fastq.gz" \
        -O "${output_dir}/${sample}_2.fastq.gz" \
        --detect_adapter_for_pe \
	--trim_poly_x \
        -l 30 -e 20 \
        --overrepresentation_analysis \
        -h "${output_dir}/${sample}_fastp.html" \
        -j "${output_dir}/${sample}_fastp.json" \
        --thread 16
done

#############################################################################

# commands for carring out salmon indexing 

# need to ensure up to date version of salmon is loaded module load salmon/v1.10.1 
# replace INDEX with chosen index name 
# used human transcriptome and the human genome for creating the decoy-aware transcriptome from GENCODE 

gunzip -c GRCh38.primary_assembly.genome.fa.gz | grep "^>" | cut -d " " -f 1 | sed 's/>//' > decoys.txt
cat gencode.v47.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
salmon index -t ../gencode_files/gentrome.fa.gz -d ../gencode_files/decoys.txt -p 12 -i INDEX --gencode

#######################################################################

#script for carrying out quantification on files to produce abundances

# input directory with fastq read files in and specify the output directory 

#!/bin/bash

input_direc="$1"
output_direc="$2"

for r1 in "$input_direc"/*_1.fastq.gz; do
    samp=$(basename "$r1" _1.fastq.gz)
    r2="${input_direc}/${samp}_2.fastq.gz"

    salmon quant -i ./gencode_files/INDEX -l A \
        -1 "$r1" -2 "$r2" \
        --gcBias \
        --seqBias \
        --numBootstraps 20 \
        -p 32 --validateMappings \
        -o "${output_direc}/${samp}_quant" \
        &> "${output_direc}/${samp}_salmon.log"
done

#######################################################################

# STAR ALIGNMENT 

# ran STAR alignment to the genome to determine how the reads map to the genome 

# generate STAR index - check the read lengths data 
# set parameter to 150 for overhang as bulk-rna dataset read lengths 151pb 
# replace INDEX with chosen index name/output directory 

STAR \
  --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir INDEX \
  --genomeFastaFiles ../PATH/TO/FILE/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile ../PATH/TO/FILE/gencode.v47.annotation.gtf \
  --sjdbOverhang 150


# script for aligning to genome

#!/bin/bash

# Set number of threads 
# specify the paths to STAR index directory, name the output directory and directory where fastq files are

THREADS=16
GENOME_DIR="/PATH/TO/INDEX/FILES/INDEX"
OUTPUT_BASE="/PATH/TO/OUTPUT_DIRECTORY"
READ_DIR="/PATH/TO/FASTQ_FILE_DIRECTORY"

module load star/v2.7

# Loop over all R1 FASTQ files with _1.fastq.gz
for R1 in "$READ_DIR"/*_1.fastq.gz
do
    # Derive sample name by removing _1.fastq.gz
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    # Set corresponding R2 file
    R2="${READ_DIR}/${SAMPLE}_2.fastq.gz"

    # Create output directory
    OUT_DIR="${OUTPUT_BASE}/${SAMPLE}"
    mkdir -p "$OUT_DIR"

    # Run STAR alignment
    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${OUT_DIR}/${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM
done

# index BAM file and move to different output folder for visualisation in IGV 

#!/bin/bash

# Find and index BAM files as needed for IGV
find . -type f -name "*sortedByCoord.out.bam" | while read -r bamfile; do
    echo "Indexing: $bamfile"
    samtools index "$bamfile"
done

#!/bin/bash

# Set the target directory
TARGET_DIR= "sorted_bams"

# Create it if it doesn't exist
mkdir -p "$TARGET_DIR"

# Find and move sorted BAM and their BAI index files
find . -type f -name "*sortedByCoord.out.bam" | while read bam; do
  mv "$bam" "$TARGET_DIR/"
  
  # Check and move associated index file if it exists
  bai="${bam}.bai"
  if [ -f "$bai" ]; then
    mv "$bai" "$TARGET_DIR/"
  fi
done

#########################################################################################
# run flagstat tools to check quality of alignment 

# script for going over samples 

#!/bin/bash

# specify directory your sorted BAM files in 
BAM_DIR="PATH/TO/FILE_DIRECTORY/sorted_bams"

# specify output directory for flagstat reports
OUT_DIR= "PATH/TO/FILE_DIRECTORY/flagstat_sorted_bams"
mkdir -p "$OUT_DIR"

# Loop over BAM files in the directory
for bam in "$BAM_DIR"/*.bam; do
    if [[ -f "$bam" ]]; then
        # Extract base name (e.g., sample1.bam → sample1)
        sample_name=$(basename "$bam" .bam)
        
        # Run samtools flagstat and save output
        samtools flagstat "$bam" > "$OUT_DIR/${sample_name}.flagstat.txt"
    fi
done

#######################################################################

# can then look on IGV Desktop with human Gencode annotation RefSeq

