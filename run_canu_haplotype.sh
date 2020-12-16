#!/bin/bash

# INPUT DATA FILES etc.

IN_CHILD=F1.fasta
IN_PARENT1=K12.parental.fasta
IN_PARENT2=O157.parental.fasta
PREFIX_PARENT1=K12
PREFIX_PARENT2=O157

CANU_WORK_DIR=ecoliTrio
CANU_COMMAND="canu -p asm -d ${CANU_WORK_DIR} useGrid=false genomeSize=5m -haplotype${PREFIX_PARENT1} ${IN_PARENT1} -haplotype${PREFIX_PARENT2} ${IN_PARENT1} -pacbio ${IN_CHILD}"

GB_MEMORY=100
N_THREADS=20

## Parameters for the modified parts (see README.md for details)
MIN_KMER_FREQ_IN_CHILD=20
MIN_DIFF_MATCH_KMERS_BETWEEN_PARENT_HAPLOTYPES=3
MIN_RATIO_MATCH_KMERS_BETWEEN_PARENT_HAPLOTYPES=2


# Set up executables: canu, splitHaplotype, meryl

## Case 1: Environment Modules
#ml canu/haplotype

## Case 2: Directly set PATH
PATH="/path/to/build/bin:$PATH"


# Run Canu as usual, which stops with an error at splitHaplotype (but we do not exit)
${CANU_COMMAND} || true


# Run meryl for the child dataset (NOTE: meryl should be that of Canu)
K=$(meryl print ${CANU_WORK_DIR}/haplotype/0-kmers/haplotype-${PREFIX_PARENT1}.meryl | head -1 | awk '{print length($1)}')
meryl count k=${K} memory=${GB_MEMORY} threads=${N_THREADS} ${IN_CHILD} output ${IN_CHILD}.meryl


# Run modified splitHaplotype
cd ${CANU_WORK_DIR}/haplotype

splitHaplotype \
  -cl 1000 \
  -cr ${MIN_RATIO_MATCH_KMERS_BETWEEN_PARENT_HAPLOTYPES} \
  -hf ${MIN_KMER_FREQ_IN_CHILD} \
  -hd ${MIN_DIFF_MATCH_KMERS_BETWEEN_PARENT_HAPLOTYPES} \
  -threads ${N_THREADS} \
  -memory ${GB_MEMORY} \
  -R ../../${IN_CHILD} \
  -M ../../${IN_CHILD}.meryl \
  -H ./0-kmers/haplotype-${PREFIX_PARENT1}.meryl ./0-kmers/reads-${PREFIX_PARENT1}.statistics ./haplotype-${PREFIX_PARENT1}.fasta.gz \
  -H ./0-kmers/haplotype-${PREFIX_PARENT2}.meryl ./0-kmers/reads-${PREFIX_PARENT2}.statistics ./haplotype-${PREFIX_PARENT2}.fasta.gz \
  -A ./haplotype-unknown.fasta.gz \
> haplotype.log.WORKING \
&& \
mv -f ./haplotype.log.WORKING ./haplotype.log

touch haplotyping.success

cd ../..


# Resume Canu
${CANU_COMMAND}
