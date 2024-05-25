# CSE185-SP24-Motif-Finding-by-Gibbs-Sampling
## Implementation
Starting with a random selection of k-mers from the sequences, the Gibbs sampler we implemented iteratively refine this selection to converge on the most probable motifs. In each iteration, a profile matrix is created from the current set of motifs, excluding one sequence. This matrix represents the probabilities of each nucleotide at each position. Using this profile, a new k-mer is probabilistically selected for the excluded sequence based on its likelihood given the profile. This process is repeated, gradually reconfines the motifs by the probabilistic framework. By running the algorithm multiple times and using a background model to calculate p-values, the implementation aims to identify statistically significant motifs that are most representative of the biological signals.
## Usage
This program is not yet developed to be a command-line executable utility. Therefore, to run this program, the myutils.py script has to be executed conventionally by using python myutils.py. To run the program as intended, the file location needs to be updated within the myutils.py script at line 105: fasta_file = "location/of/peakdatafile.fa". Parameters such as k, t, and n, corresponding to motif length, number of sequences, and number of iterations, can be varied to explore the resulting motifs.
## Coarse Benchmarking
A BED file of ChIP sequencing, `ENCFF221BLR.bed`, was downloaded from the ENCODE platform. The biosample for obtaining the peak data is *Homo sapiens* GM23338, and the sequencing was targeted to the STAT1 assay. The command `bedtools getfasta -fi hg38.fa -bed ENCFF221BLR.bed -fo peaks_sequences.fa` was used to extract the peak sequence segments from the reference *Homo sapiens* genome. The `peaks_sequences.fa` file was then used as the input for the Gibbs sampler program. When running for a k-mer length of 8, the output is:

```bash
Best motifs found:
Motif: CCCCGCCC, P-value: 0.0000000000
Motif: CCCCGCCC, P-value: 0.0000000000
Motif: CTCCGCCC, P-value: 0.0000000000
Motif: CCCCGCCC, P-value: 0.0000000000
Motif: CCGCGGCC, P-value: 0.0000000000
```

When running for a k-mer length of 14, the output is:

```bash
Best motifs found:
Motif: CGGCAGGCCTGCCT, P-value: 0.0000000000
Motif: CGCCCGCCCTGCCC, P-value: 0.0000000000
Motif: CCCCTGGCCTACCC, P-value: 0.0000000000
Motif: CCCCAGGCCTGCCC, P-value: 0.0000000000
Motif: GGCCAGGAATGCCT, P-value: 0.0000000000
```

However, variations were observed when these results were compared to the top-5 motif finding results by the HOMER tool:

```bash
Motif: TTGGCCCCGCCCCC, P-value: 1e-22
Motif: CTTCCCAG, P-value: 1e-10
Motif: CAGCTCAG, P-value: 1e-8
Motif: TTCAATTT, P-value: 1e-7
Motif: TTCTTTTT, P-value: 1e-7
```

There are several hypothesized reasons why the Gibbs sampler program might return different results compared to the HOMER tool. First, the underlying algorithms between the two tools differ; the Gibbs sampler relies on a probabilistic approach, while HOMER uses a deterministic approach based on position weight matrices (PWMs) and optimized scoring functions. This fundamental difference can lead to variations in the motifs identified. Second, the handling of background models and pseudocounts might differ between the tools, affecting the scoring and selection of motifs. Additionally, the initial conditions and random seed settings in the Gibbs sampler can introduce variability in the results, while HOMER might have more standardized or optimized initializations.

To improve and confirm the validity of our program, we can run more tests on different samples with more iterations. This will help us eliminate errors in implementation and understand the real differences between motif-finding mechanisms. Finally, there are many samples available for calling methods, and benchmarking can be improved in the following weeks.
