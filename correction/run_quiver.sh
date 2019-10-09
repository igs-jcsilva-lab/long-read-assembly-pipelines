#! /bin/bash

###############################################################################

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace

###############################################################################
    #ABOUT SCRIPT

#August 8, Kara moser

#INPUT1: Uncorrected assembly from canu
#	-Also requires the .bax.h5 files from the raw PacBio Data
#	 (Listed in Input.fofn)
#	-Quiver is part of the smrtanalysis toolkit
#	-While this is installe din GRC space, it's quicker to run ourselves (normally)
#	-Jeff/Eric set up a machine with smrtanalysis installed
#	-For instructions on logging in and setting up the env, see CVD-125

#OUTPUT: Polished assembly with quiver 2x

#Notes: if the below pipeline does not work (specifically an error occurs in pbalign):
#       Run pbalign on individual cells with separate fofns
#       Merge: cmph5tools.py merge --outFile out_all.cmp.h5 <list of cmp.h5s>
#       Sort: cmph5tools.py sort --outFile out_all.sort.cmp.h5 out_all.cmp.h5

#       Additionally, Quiver is now outdated. Future projects should use blasr. 

###############################################################################

#change as needed

sample=$1
country=$2

#work_d=/usr/local/packages/smrtanalysis/test/${sample}.canu
work_d=$3

#data directories
data_d=/local/projects-t3/p_falciparum/samples/$country/$sample

###############################################################################

mkdir -p $work_d
cd $work_d

#Getting data

f=raw.contigs.fasta
if [[ ! -e $f ]]; then
  ln -s $data_d/canu_assembly/raw.contigs.fasta $work_d/
fi

f=Input.fofn
if [[ ! -e $f ]]; then
  for i in $data_d/PACBIO_DATA/*bax.h5; do echo "${i}" >> $f;
   done
fi

###############################################################################

#Polish with Quiver

f=polished_assembly.quiver1.fasta
if [[ ! -e $f ]]; then
  samtools faidx raw.contigs.fasta

  pbalign --nproc 13 --forQuiver Input.fofn raw.contigs.fasta reads_vs_raw_assembly.cmp.h5

  quiver -j 16 reads_vs_raw_assembly.cmp.h5 -r raw.contigs.fasta -o $f -o polished_assembly.quiver1.variants.gff
fi

f=polished_assembly.quiver2.fasta
if [[ ! -e $f ]]; then
  samtools faidx polished_assembly.quiver1.fasta

  pbalign --nproc 13 --forQuiver Input.fofn polished_assembly.quiver1.fasta reads_vs_raw_assembly.quiver1.cmp.h5

  quiver -j 16 reads_vs_raw_assembly.quiver1.cmp.h5 -r polished_assembly.quiver1.fasta -o $f -o polished_assembly.quiver2.variants.gff
fi

###############################################################################

echo "Quiver complete! Next step: correction with pilon."
echo "run_pilon_and_comparison_stats.sh"

###############################################################################

echo "Read through script."

###############################################################################
