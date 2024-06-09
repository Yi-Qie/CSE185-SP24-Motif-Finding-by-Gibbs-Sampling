# Gibbs Sampling Motif Finding (gsmf)
Starting with a random selection of k-mers from the sequences, the Gibbs sampler we implemented iteratively refine this selection to converge on the most probable motifs. In each iteration, a profile matrix is created from the current set of motifs, excluding one sequence. This matrix represents the probabilities of each nucleotide at each position. Using this profile, a new k-mer is probabilistically selected for the excluded sequence based on its likelihood given the profile. This process is repeated, gradually reconfines the motifs by the probabilistic framework. By running the algorithm multiple times and using a background model to calculate p-values, the implementation aims to identify statistically significant motifs that are most representative of the biological signals.

# Usage
The gsmf (Gibbs Sampler Motif Finder) tool allows you to search for motifs in a set of DNA sequences using the Gibbs sampling algorithm. Below are the options and an example command to run the tool.

## Installation
1. Clone the repository: `git clone https://github.com/Yi-Qie/CSE185-SP24-Motif-Finding-by-Gibbs-Sampling.git`
2. Navigate to the project directory: `cd gsmf`
3. Install the package: `pip install .`

Once installed, you can use the gsmf command as indicated in the options session.

Requirements
Make sure you have the following dependencies installed:

* Python 3.x
* argparse module (usually included with Python)
* random module (usually included with Python)
* bedtools (if using bed file extraction functionality)

### Installing dependencies:
* Installing Python 3.x
Check Python Version: First, check if Python is already installed on your system and its version. Open a terminal or command prompt and type: `python --version`. If Python is installed, the version number will be displayed. Ensure it is Python 3.x (e.g., 3.7, 3.8).

If Python 3.x is not installed or if you need to update, download the latest version from the official Python website and follow the installation instructions for your operating system.

* Installing argparse Module: The argparse module is included with Python, so no separate installation is required.

* Installing random Module: The random module is also included with Python, so no separate installation is required.

* Installing bedtools (Optional): If you intend to use the bed file extraction functionality, you will need to install bedtools. Follow these steps:
1. Download bedtools: Visit the bedtools website and download the appropriate version for your operating system.
2. Install bedtools: Follow the installation instructions provided for your operating system.


Once you have installed Python 3.x and, optionally, bedtools, you are ready to use the gsmf tool. Follow the instructions provided in the repository's README or documentation to clone the repository, install the package, and use the tool.

Note: If the `gsmf` command was not found, you may need to run export `PATH=$PATH:/home/$USER/.local/bin` to include the script installation path in your $PATH.

## Options
The tool takes 5 mandatory input options:

* `-f`: specifies the path to the input FASTA file from which motifs will be searched, please make sure that you indicate a full path;

Note that this should be the extracted peak regions of complete reference genome, you would find instruction for extraction with `bed` file and `fasta` reference genome file in the end of options section.
* `-k`: length of the motif you're interested in;
* `-t`: number of sequences you would like the tool to output, the tool will print `t` sequences that strongly agree with the motif consensus when it finishes motif searching;
* `-n`: number of iterations you would like the tool to run, please note that for each iteration the algorithm would generate 1000 profile matrices by randomly replacing one of the motifs selected. We generally recommend `n` value below 1000, a larger value should also work but it is likely to take a longer time.

**Example Command:** `gsmf -f peak_sequences.fasta -t 5 -k 8 -n 1000`
This command will run the gsmf tool using:

The input FASTA file peak_sequences.fasta
5 sequences to be used in the motif search process
A motif length of 8
1000 iterations for the Gibbs sampler algorithm

\
To extract the peak regions from `fasta` reference genome using `bed` file information, you could use `bedtools`:
1. Install `bedtools` following [instruction](https://bedtools.readthedocs.io/en/latest/content/quick-start.html)
2. Use `bedtools` with this command: `bedtools getfasta -fi <reference_genome_file.fa> -bed <bed_file.bed> -fo <output_peak_region.fa>`

# Benchmark
Please find our benchmark documentation [here](/benchmark/Benchmark.md). 

# Contributors
This repository was created as part of a coursework project for CSE 185 Advanced Bioinformatics Lab at University of California San Diego. We would like to acknowledge the efforts and contributions of the following team members:

- **Qie Yi** - Contributed to data analysis, main algorithm implementation, and documentation.
- **Zhuoling Huang** - Contributed to testing, algorithm optimization, benchmarking and helped with code debugging.

We are grateful for the guidance and support provided by our course instructor, Professor Melissa Gymrek, and teaching assistants, Divya Prabhu and Hao Xu.

Special thanks to [ENCODE](https://www.encodeproject.org/) for providing test datasets, [HOMER](http://homer.ucsd.edu/homer/) and [RSAT](http://rsat.sb-roscoff.fr/) for benchmarking.

For any questions or suggestions, please feel free to contact us.
