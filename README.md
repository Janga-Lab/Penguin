# Penguin
Penguin: A Tool for Predicting Pseudouridine Sites in Direct RNA Nanopore Sequencing Data

# Getting Started and pre-requisites
The following softwares and modules should be installed before using Penguin

python 3.6.10

minimpa2 (https://github.com/lh3/minimap2)

Nanopolish (https://github.com/jts/nanopolish)

samtools (http://www.htslib.org/)

numpy 1.18.1

pandas 1.0.1

sklearn 0.22.2.post1

tensorflow 2.0.0

keras 2.3.1 (using Tensorflow backend)


# Running Penguin:

In order to run Penguin, the user has do the following:

1- Ensure that the bedfile in the same path where penguin main.py file exists:
2- Run the following python command:

python main.py -r ref.fa -f reads.fastq

Where the penguin tool needs the following two inputs files when running it:

- A reference Genome file (ref.fa)
- The fastq reads file (reads.fastq)

# Note:
The user should enter the bed file name with the absolute path and extension 
