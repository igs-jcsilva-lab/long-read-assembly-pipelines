#!/bin/bash

###############################################################################

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace

###############################################################################

export PATH=/usr/local/bin:$PATH

java=/etc/alternatives/jre_1.8.0/bin/java

canu=/usr/local/packages/canu-1.3/bin/canu

###############################################################################

#January 2017, Kara Moser

#This is the same script as run_canu.sh,
#   but with correct paths to run on the new rhel7 machines
#   qlogin -q interactive.q -P jcsilva-gcid-proj4a-malaria -l centos7,mem_free=30G
#   (can also use a specific host: -l hostname=grid-1-3-1

#Input: Sample and corresponding country
#       -Assumes directory heirachy as in p_falciparum/samples/$country/$sample/$PACBIO_DATA
#       -How Canu is supposed to run: processes will be submitted by Canu to the grid automatically
#       -Canu is able to auto-detect resource requirements and scale itself to fit
#       -How Canu currently actually runs: Takes out nodes on the grid

#       -Until this is fixed, open up a qlogin session with at least ~30 GB

#       -Run this in an area with LOTS of space (scratch is prob. the best bet)

#Output: Raw assembly, ready for downstream analysis

###############################################################################

sample=$1
country=$2

###############################################################################

canu_d=${sample}_canu_assembly
mkdir -p $canu_d

f=$sample
if [[ ! -e $f ]]; then
  ln -s /local/projects-t3/p_falciparum/samples/$country/$sample
fi

###############################################################################

$canu -p $sample -d $canu_d genomeSize=24m -pacbio-raw ${sample}/PACBIO_DATA/*.fastq \
  corMaxEvidenceErate=0.15 \
  corMhapSensitivity=high \
  java=$java \
  useGrid=0 \
  ovsMemory=8g-256g
#  maxThreads=8

##DO NOT RUN BELOW CODE UNTIL THE GRID ISSUE IS RESOLVED

#$canu -p $sample -d $canu_d genomeSize=24m \
#  -pacbio-raw ${sample}/PACBIO_DATA/*.fastq \
#  gridOptionsJobName=canu_${sample} \
#  gridOptions="-P jcsilva-gcid-proj4a-malaria" \
#  gridEngineThreadsOption="-pe make THREADS" \
#  gridEngineMemoryOption="-l mem_free=MEMORY" \
#  corMaxEvidenceErate=0.15 \
#  corMhapSensitivity=normal java=$java \
#  gridOptions="-tc 10" \
#  maxThreads=2

###############################################################################

echo "Canu complete! Next step: polishing with quiver"
echo "run_quiver.sh"

###############################################################################

echo "Read through script."

###############################################################################
