# Canu + modified haplotype separation for better trio binning

## Modifications and rationale

This branch (`haplotype-separation`; default branch of this fork) contains two modifications in the `splitHapplotype` program, which performs haplotype separation of child reads using parent datasets, as follows:

1. Ignore k-mers that appear in the child dataset with a frequency smaller than some threshold (`-hf` option; 20 by default). This removes potential false positive k-mer matches between child and parent datasets due to sequencing errors (Recall the spike around x=1 in the GenomeScope plot; we discard that part). This value should be determined by looking at the k-mer frequency histogram of the child dataset (using meryl or similar tools).
2. For each child read, Canu by default assigns one of the haplotypes (=parents) that has a larger number of k-mers shared with the child read. Our modification requires the difference in the number of shared k-mers between the two haplotypes to be greater than some threshold (`-hd` option; 5 by default). A larger threshold makes the decision of haplotype separation more conservative (for example, a read having only 2 mother-mers and 1 father-mer will be classified as ambiguous with the default threshold whereas original Canu regards it as a mother-specific read).

Currently these modifications are not incorporated into the whole pipeline of Canu but into only `splitHaplotype`, and one needs to use a script `run_canu_haplotype.sh` (put in the root directory of this branch) to run it.

## how to install/use

1. Download and compile:

```bash
git clone https://github.com/yoshihikosuzuki/canu
cd canu/src
make -j<number of threads>
```

This will generate a directory named `build/` in the root directory of this repository.

2. Modify the script `run_canu_haplotype.sh` (put in the root directory) according to your input files and environments:

```bash
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
```

---

# Canu 

Canu is a fork of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page), designed for high-noise single-molecule sequencing (such as the [PacBio](http://www.pacb.com) [RS II](http://www.pacb.com/products-and-services/pacbio-systems/rsii/)/[Sequel](http://www.pacb.com/products-and-services/pacbio-systems/sequel/) or [Oxford Nanopore](https://www.nanoporetech.com/) [MinION](https://nanoporetech.com/products)).

Canu is a hierarchical assembly pipeline which runs in four steps:

* Detect overlaps in high-noise sequences using [MHAP](https://github.com/marbl/MHAP)
* Generate corrected sequence consensus
* Trim corrected sequences
* Assemble trimmed corrected sequences

## Install:

* Do NOT download the .zip source code.  It is missing files and will not compile.  This is a [known flaw](https://github.com/dear-github/dear-github/issues/214) with git itself.

* The easiest way to get started is to download a binary [release](http://github.com/marbl/canu/releases).

* Installing with a 'package manager' is not encouraged, but if you have no other choice:
  * Conda: `conda install -c conda-forge -c bioconda -c defaults canu`
  * Homebrew: `brew install brewsci/bio/canu`

* Alternatively, you can use the latest unreleased version from the source code.  This version has not undergone the same testing as a release and so may have unknown bugs or issues generating sub-optimal assemblies. We recommend the release version for most users.

        git clone https://github.com/marbl/canu.git
        cd canu/src
        make -j <number of threads>

 * An *unsupported* Docker image made by Frank Förster is at https://hub.docker.com/r/greatfireball/canu/.

## Learn:

The [quick start](http://canu.readthedocs.io/en/latest/quick-start.html) will get you assembling quickly, while the [tutorial](http://canu.readthedocs.io/en/latest/tutorial.html) explains things in more detail.

## Run:

Brief command line help:

    ../<architecture>/bin/canu

Full list of parameters:

    ../<architecture>/bin/canu -options

## Citation:
 - Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM. [Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation](https://doi.org/10.1101/gr.215087.116). Genome Research. (2017). `doi:10.1101/gr.215087.116`
 - Koren S, Rhie A, Walenz BP, Dilthey AT, Bickhart DM, Kingan SB, Hiendleder S, Williams JL, Smith TPL, Phillippy AM. [De novo assembly of haplotype-resolved genomes with trio binning](http://doi.org/10.1038/nbt.4277).  Nature Biotechnology.  (2018). (If you use trio-binning)
 - Nurk S, Walenz BP, Rhiea A, Vollger MR, Logsdon GA, Grothe R, Miga KH, Eichler EE, Phillippy AM, Koren S. [HiCanu: accurate assembly of segmental duplications, satellites, and allelic variants from high-fidelity long reads](https://doi.org/10.1101/2020.03.14.992248).  biorXiv.  (2020). (If you use -pacbio-hifi)
 

