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

#Here should be the final corrected assembly, final_assembly.fasta
mkdir -p $data_d/canu_assembly/compare_to_3d7
cd $data_d/canu_assembly/compare_to_3d7

f=assembly.fasta
if [[ ! -e $f ]]; then
  ln -s $data_d/canu_assembly/polished_assembly.quiver2.pilon.fasta $data_d/canu_assembly/compare_to_3d7/$f
fi

#Generating whole-genome stats compared to 3D7

tag=against_reference
f=$tag.delta
if [[ ! -e $f ]]; then
  echo 1>&2 f=$f

  nucmer -p $tag $d/reference.fa assembly.fasta
#  gzip -v final_assembly.fasta
fi

f=$tag.r.filtered.delta
if [[ ! -e $f ]]; then
  echo 1>&2 f=$f

  delta-filter -r $tag.delta \
    > $f
fi

f=position_identities.gz
if [[ ! -e $f ]]; then
  echo 1>&2 f=$f

  show-coords -cHlT $tag.r.filtered.delta \
    | /home/kara.moser/codes/collect_identity \
    | sort \
    | cut -f 3 \
    | gzip -v \
    > $f
fi

f=r.genome_coverage
if [[ ! -e $f ]]; then
  echo 1>&2 f=$f

  zcat position_identities.gz \
    | awk -F "\t" -v OFS="\t" '{print $1 != 0.0}' \
    | /home/edrabek/bin/groupby 0 avg \
    > $f
fi

f=r.mean_identity
if [[ ! -e $f ]]; then
  echo 1>&2 f=$f

  zcat position_identities.gz \
    | awk -F "\t" -v OFS="\t" '{print "all", $1; if ($1 != 0.0) {print "covered", $1}}' \
    | sort \
    | /home/edrabek/bin/groupby -c 1 avg \
    > $f
fi

f=bottom_line
if [[ ! -e $f ]]; then
  < r.mean_identity /home/edrabek/bin/eq 1 all \
    | awk -F "\t" -v OFS="\t" '{print 100 - $3}' \
    > $f
fi

###############################################################################

#Generating gene-based stats compared to 3D7

f=$tag.r.filtered.bed
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

        show-coords -cHlT $tag.r.filtered.delta \
          | tabawk '{if ($2 < $1) {temp = $2; $2 = $1; $1 = temp}; print $12, $1 - 1, $2}' \
          > $f
      fi

      f=$tag.r.filtered.bedtools_out
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

        bedtools coverage -a $tag.r.filtered.bed -b $d/reference.gff \
          > $f
      fi

      f=mrna-covered_bases-total_bases-covered_fraction
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

      < $tag.r.filtered.bedtools_out grep -v '^#' \
          | eq 3 CDS \
          | cut -f 9,11,12 \
          | perl -pe 's/^\S+;Parent=(\S+?)(;\S*)?\t/$1\t/ or die $_' \
          | sort \
          | groupby 1 sum sum \
          | tabawk '{print $0, $2 / $3}' \
          > $f
      fi

      f=$tag.r.filtered.show-snps_out
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

        show-snps -CHT $tag.r.filtered.delta \
          > $f
      fi

      f=$tag.r.filtered.snps.vcf
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

        {
          echo '#CHROM        POS     ID      REF     ALT     QUAL    FILTER INFO    FORMAT  NF54'
          < $tag.r.filtered.show-snps_out neq 2 . \
            | neq 3 . \
            | tabawk '{print $9, $1, ".", $2, $3, "100", "PASS", ".", "GT", "1"}'
        } \
            > $f
      fi

      f=$tag.r.filtered.raw
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

        < $tag.r.filtered.show-snps_out /home/edrabek/daily/1604/18/collapse_indels \
          | /home/edrabek/daily/1604/18/get_sillab-122_stats $d/reference.{fa,gff} /dev/null \
          | sort \
          | groupby -c 2 \
          > $f
      fi

      f=$tag.r.filtered.variant_counts
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

        < $tag.r.filtered.raw cut -f 2- \
          | sort \
          | groupby 1 sum \
          > $f
      fi

      f=$tag.r.filtered.gene_table
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

        < $tag.r.filtered.raw perl -pe 's/\b(non)?coding\.//' \
          | perl -pe 's/\.(non)?coding\b//' \
          | sort -k 2,2 \
          | project-table -f mrna -d 0 \
          | neq 1 None \
          > $f
      fi

f=$tag.r.filtered.gene_counts
      if [[ ! -e $f ]]; then
        echo 1>&2 f=$f

        < $tag.r.filtered.gene_table transpose \
          | tail -n +2 \
          | tabawk '{count = 0; for (field = 2; field <= NF; field++) {count += $field >= 1}; print $1, count}' \
          | cat - <(< mrna-covered_bases-total_bases-covered_fraction lt 4 0.99 | wc -l | prefix not_fully_covered) \
          > $f
      fi


f=variants.vcf
if [[ ! -e $f ]]; then
  show-snps -HCT against_reference.r.filtered.delta \
  | /home/edrabek/daily/1605/16/collapse_indels \
  | awk -F "\t" -v OFS="\t" '{print $9, $1, ".", $2, $3, "100", "PASS", ".", "GT", "1"}' \
  | /home/edrabek/bin/label-columns '#CHROM' POS ID REF ALT QUAL FILTER INFO FORMAT ${sample} \
  > variants.vcf
fi

f=variants.allsites.vcf
if [[ ! -e $f ]]; then
  show-coords -cHlT against_reference.r.filtered.delta \
  | /home/edrabek/daily/1605/16/fill_in_nucmer_to_vcf variants.vcf $d/reference_cleaned.fa \
  > $f
fi

###############################################################################

#Annotate variants

tag="variants"
f=$tag.annotated.html
if [[ ! -e $f ]]; then
  ${java} -Xmx8g -jar $snpeff/snpEff.jar \
  -c $snpeff/snpEff.config Pf3D7v24 \
  -v $tag.vcf \
  > $tag.annotated.vcf -stats $tag.annotated.html
fi

f=variants.snps.vcf
if [[ ! -e $f ]]; then
  echo "##fileformat=VCFv4.1" > $f;
  < variants.vcf grep "\#" >> $f;
  
  < variants.vcf awk '$4 != "."' \
  | awk '$5 != "."' \
  >> $f
fi

export PATH=/usr/local/packages/vcftools/bin:$PATH
export PATH=/usr/local/packages/tabix:$PATH
export PERL5LIB=/usr/local/packages/vcftools/perl:$PERL5LIB

f=variants.snps.snpden
if [[ ! -e $f ]]; then
  vcftools --vcf variants.snps.vcf --SNPdensity 1000 --out variants.snps 
fi 

###############################################################################

echo "Pilon and comparison stats complete!"
echo "This is more or less the end of the assembly pipeline."
echo "Inspect assembly and make changes as required (order, revcomp, etc)"
echo "If nothing else is required, you can annotate your assembly!"
echo "run_annotation.sh"

###############################################################################

echo  "Read through script."

###############################################################################











