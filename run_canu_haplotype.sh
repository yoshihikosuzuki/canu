#!/bin/bash

IN_CHILD=ONT.fastq.gz   # TODO: meryl accepts .gz?
IN_F1=F1.fastq.gz
IN_M2=M2.fastq.gz

CANU_WDIR=oikTrio

MIN_KMER_FREQ_IN_CHILD=20
MIN_DIFF_MATCH_KMERS_BETWEEN_PARENT_HAPLOTYPES=5

#ml canu/haplotype

meryl ${IN_CHILD} > ${IN_CHILD}.meryl   # NOTE: meryl should be that of Canu

# Run Canu as usual, which throws an error at splitHaplotype
canu ... || true

# Run modified splitHaplotype
cd ${CANU_WDIR}/haplotype

splitHaplotype \
  -cl 1000 \
  -hf ${MIN_KMER_FREQ_IN_CHILD} \
  -hd ${MIN_DIFF_MATCH_KMERS_BETWEEN_PARENT_HAPLOTYPES} \
  -threads 128 \
  -memory  500 \
  -R ../../${IN_CHILD} \
  -M ../../${IN_CHILD}.meryl \
  -H ./0-kmers/haplotype-F1.meryl ./0-kmers/reads-F1.statistics ./haplotype-F1.fasta.gz \
  -H ./0-kmers/haplotype-M2.meryl ./0-kmers/reads-M2.statistics ./haplotype-M2.fasta.gz \
  -A ./haplotype-unknown.fasta.gz \
> haplotype.log.WORKING \
&& \
mv -f ./haplotype.log.WORKING ./haplotype.log

cd ../..

# Resume Canu
canu ...   
