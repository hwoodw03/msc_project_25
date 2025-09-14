# make directory for each sample e.g. CP1

mkdir CP1

# select fastq files from ENA with download script on database 
# select SAMPLE_NAME_barcodes.tsv.gz from GEO 

# ensure one cell barcode on each line and appears without hyphens/spaces  
gunzip SAMPLE_NAME_barcodes.tsv.gz
sed 's/-1$//' SAMPLE_NAME_barcodes.tsv > barcodes_stripped.tsv

# then have directory for each sample with cell barcodes and corresponding FILE_1.fastq.gz that is the cell barcode + UMI 
# and FILE_2.fastq.gz which is the correpsonding reads

###############################################################################################
# Script for iterating fastqc over samples. Need to ensure that load module fastqc and parallel. 

# input directory where the fastq files are

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
echo "FASTQC analysis complete." 

##################################################################
# Commands for generating salmon index and running salmon Alevin 

# need to ensure up to date version of salmon is loaded module load salmon/v1.10.1 
# replace INDEX with chosen index name output directory
# used human transcriptome and the human genome for creating the decoy-aware transcriptome from GENCODE 

gunzip -c GRCh38.primary_assembly.genome.fa.gz | grep "^>" | cut -d " " -f 1 | sed 's/>//' > decoys.txt
cat gencode.v47.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz

salmon index -t ../gencode_files/gentrome.fa.gz -d ../gencode_files/decoys.txt -p 12 -i INDEX --gencode

# as doing on a transcript-level rather than a gene-level use transcript to transcript map 
# rather than a transcript to gene map

# used R with package DTUrtle to create txmap 
tx2gene <- import_gtf(gtf_file = "../gencode.vM24.annotation.gtf")

write.table(x = tx2gene[,c("transcript_id", "transcript_name")], 
            file = "txmap.tsv", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# alevin transcript abundance quantification and repeat for each sample with command below 
# make direcotry for quantification 

mkdir whole_transcriptome_alevin
cd whole_transcriptome_alevin
salmon alevin -l ISR -1 ../SAMPLE_DIRECTORY/FILE_1.fastq.gz -2 ./SAMPLE_DIRECTORY/FILE_2.fastq.gz --chromiumV3 -i ./PATH/TO/SALMON/INDEX -p 16 --dumpFeatures --whitelist ../SAMPLE_NAME/barcodes_stripped.tsv -o SAMPLE_NAME_alevin --tgMap ./PATH/TO/FILE/txmap.tsv


##########################################
