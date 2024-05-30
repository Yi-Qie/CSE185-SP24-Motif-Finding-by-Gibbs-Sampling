## Data information
A BED file of ChIP sequencing, [`ENCFF221BLR.bed`](https://www.encodeproject.org/files/ENCFF221BLR/), was downloaded from the ENCODE platform. The biosample for obtaining the peak data is [*Homo sapiens* GM23338](https://www.encodeproject.org/biosamples/ENCBS368AAA/), and the sequencing was targeted to the STAT1 assay. 

`bedtools` is incorporated to extract peak sequence segments from reference genome using coordinate information provided by bed file. For instance, with the two example files we mentioned above, the `bedtools` command is: `bedtools getfasta -fi hg38.fa -bed ENCFF221BLR.bed -fo peaks_sequences.fa`. 

Where the ouput `peaks_sequences.fa` file was then used as the input for the Gibbs sampler program. When running for a k-mer length of 8, the output is:

```bash
Best motifs found:
Motif: CCCCGCCC, P-value: 0.0000000000
Motif: CCCCGCCC, P-value: 0.0000000000
Motif: CTCCGCCC, P-value: 0.0000000000
Motif: CCCCGCCC, P-value: 0.0000000000
Motif: CCGCGGCC, P-value: 0.0000000000
```

## Benchmark Against HOMER
Below are the top 5 motifs found with lowest p-values among all using HOMER tool:

```bash
Motif: TTGGCCCCGCCCCC, P-value: 1e-22
Motif: CTTCCCAG, P-value: 1e-10
Motif: CAGCTCAG, P-value: 1e-8
Motif: TTCAATTT, P-value: 1e-7
Motif: TTCTTTTT, P-value: 1e-7
```
There is a rather noticeable difference between the motifs found by HOMER and our tool. 
There are several hypothesized reasons why the Gibbs sampler program might return different results compared to the HOMER tool. First, the underlying algorithms between the two tools differ; the Gibbs sampler relies on a probabilistic approach, while HOMER uses a deterministic approach based on position weight matrices (PWMs) and optimized scoring functions. This fundamental difference can lead to variations in the motifs identified. Second, the handling of background models and pseudocounts might differ between the tools, affecting the scoring and selection of motifs. Additionally, the initial conditions and random seed settings in the Gibbs sampler can introduce variability in the results, while HOMER might have more standardized or optimized initializations.

## Benchmark Against Regulatory Sequence Analysis Tools (RSAT)
[RSAT](http://rsat.sb-roscoff.fr/info-gibbs_form.cgi) is an online tool that provides a series of modular computer programs specifically designed for the detection of regulatory signals. One of its features is motif finding with Gibbs sampling approach, therefore we used it to benchmark the accuarcy of our tool and investigate the difference between HOMER result and result produced by our tool. 
It is worth noting that, to our surprise, we didn't find a commonly used established package/tool for Gibbs sampling-based motif finding implemented in python. 

We use the `peaks_sequences.fa` file as sequence input to run RSAT, which outputs the following logos:

![RSAT_logo_1](benchmark/RSAT_logos/RSAT_logo_1.png)

![RSAT_logo_2](benchmark/RSAT_logos/RSAT_logo_2.png)

The rest of parameters can be found at the end of this benchmark session. 
