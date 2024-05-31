# Gibbs Sampling Motif Finding
## Gibbs Sampling
Starting with a random selection of k-mers from the sequences, the Gibbs sampler we implemented iteratively refine this selection to converge on the most probable motifs. In each iteration, a profile matrix is created from the current set of motifs, excluding one sequence. This matrix represents the probabilities of each nucleotide at each position. Using this profile, a new k-mer is probabilistically selected for the excluded sequence based on its likelihood given the profile. This process is repeated, gradually reconfines the motifs by the probabilistic framework. By running the algorithm multiple times and using a background model to calculate p-values, the implementation aims to identify statistically significant motifs that are most representative of the biological signals.

## Installation


## Usage
This program is not yet developed to be a command-line executable utility. Therefore, to run this program, the myutils.py script has to be executed conventionally by using python myutils.py. To run the program as intended, the file location needs to be updated within the myutils.py script at line 105: fasta_file = "location/of/peakdatafile.fa". Parameters such as k, t, and n, corresponding to motif length, number of sequences, and number of iterations, can be varied to explore the resulting motifs.

## Options
The tool takes 5 mandatory input options:

* `-bed`: the bed file which contains coordinates information of the motif you're interested in. It serves as the guide to extract sequences at specific locations from the reference genome;
* `-fa`: the reference genome file which you would like to search motif from, please makesure it's in `fasta` format;
* `-k`: length of the motif you're interested in;
* `-t`: number of sequences you would like the tool to output, the tool will print `t` sequences that strongly agree with the motif consensus when it finishes motif searching;
* `-n`: number of iterations you would like the tool to run, please note that for each iteration the algorithm would generate 1000 profile matrices by randomly replacing one of the motifs selected. We generally recommend `n` value below 1000, a larger value should also work but it is likely to take a longer time.


## Contributors
This repository was created as part of a coursework project for CSE 185 Advanced Bioinformatics Lab at University of California San Diego. We would like to acknowledge the efforts and contributions of the following team members:

- **Qie Yi** - Contributed to data analysis, main algorithm implementation, and documentation.
- **Zhuoling Huang** - Contributed to testing, algorithm optimization, benchmarking and helped with code debugging.

We are grateful for the guidance and support provided by our course instructor, Professor Melissa Gymrek, and teaching assistants, Divya Prabhu and Hao Xu.

Special thanks to [ENCODE](https://www.encodeproject.org/) for providing test datasets, [HOMER](http://homer.ucsd.edu/homer/) and [RSAT](http://rsat.sb-roscoff.fr/) for benchmarking.

For any questions or suggestions, please feel free to contact us.
