#! /bin/bash

###############################################################################

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace

export PATH=.:/home/edrabek/bin:${PATH:-}
export PYTHONPATH=/home/edrabek/lib/python:${PYTHONPATH:-}
export LC_ALL=C

###############################################################################
    #ABOUT SCRIPT

#August 2018, Kara Moser

#modified from: /home/edrabek/daily/1601/25/do_pilon_corrections
#               /home/edrabek/daily/1509/21/do_eplas_72

#INPUT1: Raw illumina reads, polished assembly (quiver, 2x)
#	-I normally run this on scratch

#OUTPUT: Illumina-corrected assembly
#        Output statistics comparing corrected assembly to 3D7

###############################################################################

#options

sample=$1
country=$2

#This will be created wherever you run the script from
work_d=${sample}_pilon

#data directories
d=/local/projects-t3/p_falciparum/auxiliary_files
data_d=/local/projects-t3/p_falciparum/samples/$country/$sample

###############################################################################

#requires at least java1.6

snpeff=/home/kara.moser/bin/snpEff

###############################################################################

mkdir -p $work_d
cd $work_d

#Get Illumina reads and assembly 
f=ILLUMINA_DATA.correct_files_to_use
if [[ ! -e $f ]]; then
  ln -s $data_d/$f/ 
fi

f=best_pacbio_assembly.fa
if [[ ! -e $f ]]; then
  ln -s $data_d/canu_assembly/polished_assembly.quiver2.fasta best_pacbio_assembly.fa
fi

###############################################################################

fastq_count=`ls ILLUMINA_DATA.correct_files_to_use/*_R1.fastq.gz | wc -l`
if [[ $fastq_count == 1 ]]; then
  lnsreal ILLUMINA_DATA.correct_files_to_use/*_R1.fastq.gz 1.fastq.gz
  lnsreal ILLUMINA_DATA.correct_files_to_use/*_R2.fastq.gz 2.fastq.gz
else
  zcat ILLUMINA_DATA.correct_files_to_use/*_R1.fastq.gz \
    | gzip -v \
    > 1.fastq.gz
#  lnsreal 1.fastq.gz 1.fastq.gz
  zcat ILLUMINA_DATA.correct_files_to_use/*_R2.fastq.gz \
    | gzip -v \
    > 2.fastq.gz
#  lnsreal 2.fastq.gz 2.fastq.gz
fi

#Align reads against assembly to be corrected for Pilon

f=$data_d/canu_assembly/final_assembly.fasta
if [[ ! -e $f ]]; then

  f=reads_vs_best_pacbio_assembly.sorted.bam
    if [[ ! -e $f ]]; then
      bowtie2-build best_pacbio_assembly.fa vs_best_pacbio_assembly

      bowtie2 --no-unal -t --rfg 0,3 -x vs_best_pacbio_assembly -1 *1.fastq.gz -2 *2.fastq.gz \
      | samtools view - -bS -h \
      > reads_vs_best_pacbio_assembly.bam
      samtools sort reads_vs_best_pacbio_assembly.bam reads_vs_best_pacbio_assembly.sorted
      samtools index reads_vs_best_pacbio_assembly.sorted.bam
  fi

  java="/usr/local/packages/jdk1.7.0_40/bin/java"

  #Running Pilon (optimized paramaters + phillipey's --fix bases option)

  f=pilon_output.fasta
  if [[ ! -e $f ]]; then
    $java -Xmx16G -jar /usr/local/packages/pilon-1.13/pilon-1.13.jar \
      --genome best_pacbio_assembly.fa \
      --output pilon_output \
      --frags reads_vs_best_pacbio_assembly.sorted.bam \
      --changes \
      --vcf \
      --fix bases \
      --mindepth 5 \
      --K 85 \
      --minmq 0 \
      --minqual 35 \
      > pilon.log

    cp $f $data_d/canu_assembly/polished_assembly.quiver2.pilon.fasta

  fi

fi

###############################################################################

echo "Pilon run complete!"
echo "Inspect assembly and make changes as required (order, revcomp, etc)"
echo "If nothing else is required, you can annotate your assembly!"
echo "run_annotation.sh"

###############################################################################

echo "Read through script."

###############################################################################
