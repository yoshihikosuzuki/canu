# Canu v2.1.1 + custom haplotype separation

## Modifications and rationale

This fork contains basically two modifications in `splitHapplotype.C`, which performs classification of each child read into one of the two parents (or "unknwon"). There are three additional command-line arguments (run `$ splitHaplotype` for details):

1. `min-kmer-freq` in the `-H` argument
  - By default, Canu automatically determines the minimum k-mer frequency used as a thoreshold for haplotype-specific k-mers for each parent read set based on the shape of the k-mer count histogram. However, it sometimes failes to detect appropriate thresholds especially when the input datasets have some problems (e.g. contamination). We added this argument so that a user can specify arbitrary thresholds determined by manual inspection of k-mer count histograms.
2. `-M child-kmers.meryl min-kmer-freq`
  - This is used for ignoring k-mers that appear in the child dataset with a frequency smaller than some threshold, `min-kmer-freq`. K-mer counts of the child dataset, i.e. `child-kmers.meryl`, must be precomputed with Meryl provided in this repository.
  - This removes potential false positive k-mer matches between child and parent datasets due to sequencing errors in child reads (Recall the spike around x=1 in the GenomeScope plot; we discard that part). This value should be determined by looking at the k-mer count histogram of the child dataset.
4. `-cn number`
  - It requires each child read to have a difference in the number of k-mers shared with each of the two haplotypes greater than `-cn` in order to be classified. A larger threshold makes the decision of haplotype separation more conservative.

Currently these modifications are not incorporated into the whole pipeline of Canu but into only `splitHaplotype`, and one needs to use a custom script `run_canu_haplotype.sh` (put in the root directory of this branch) to run it.

## how to install/use

1. Download and compile:

```bash
git clone https://github.com/yoshihikosuzuki/canu
cd canu/src
make -j<number of threads>
```

This will generate a directory named `build/` in the root directory of this repository.

2. Modify the script `run_canu_haplotype.sh` (put in the root directory) according to your input files and environments and run it.


---


# Canu

Canu is a fork of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page), designed for high-noise single-molecule sequencing (such as the [PacBio](http://www.pacb.com) [RS II](http://www.pacb.com/products-and-services/pacbio-systems/rsii/)/[Sequel](http://www.pacb.com/products-and-services/pacbio-systems/sequel/) or [Oxford Nanopore](https://www.nanoporetech.com/) [MinION](https://nanoporetech.com/products)).

Canu is a hierarchical assembly pipeline which runs in four steps:

* Detect overlaps in high-noise sequences using [MHAP](https://github.com/marbl/MHAP)
* Generate corrected sequence consensus
* Trim corrected sequences
* Assemble trimmed corrected sequences

## Install:

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
 

