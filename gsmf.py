#!/usr/bin/env python3
import random
import argparse

def profile_randomly_generated_kmer(sequence, k, profile):
    sequence = sequence.upper()
    probs = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        prob = 1
        for j in range(k):
            prob *= profile[j][kmer[j]]
        probs.append(prob)
    total_prob = sum(probs)
    probs = [p / total_prob for p in probs]
    return random.choices([sequence[i:i+k] for i in range(len(sequence) - k + 1)], weights=probs, k=1)[0]

def create_profile_with_pseudocounts(motifs, k):
    profile = [{'A': 1, 'C': 1, 'G': 1, 'T': 1} for _ in range(k)]
    for motif in motifs:
        motif = motif.upper()
        for i in range(k):
            profile[i][motif[i]] += 1
    for pos in profile:
        total = sum(pos.values())
        for nucleotide in pos:
            pos[nucleotide] /= total
    return profile

def score(motifs, k):
    consensus = ''
    for i in range(k):
        freq = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            motif = motif.upper()
            freq[motif[i]] += 1
        consensus += max(freq, key=freq.get)
    score = 0
    scores = [0 for _ in range(len(motifs))]
    for motif in motifs:
        motif = motif.upper()
        for i in range(k):
            if motif[i] != consensus[i]:
                score += 1
        scores[motifs.index(motif)] = score
    return score

def score_single_motif(motifs, k):
    consensus = ''
    for i in range(k):
        freq = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            motif = motif.upper()
            freq[motif[i]] += 1
        consensus += max(freq, key=freq.get)
    score = 0
    scores = [0 for _ in range(len(motifs))]
    for motif in motifs:
        motif is motif.upper()
        for i in range(k):
            if motif[i] != consensus[i]:
                score += 1
        scores[motifs.index(motif)] = score
    return scores

def gibbs_sampler(dna, k, t, n):
    best_motifs = []
    best_scores = []
    best_score = float('inf')

    for _ in range(1000):  # Run Gibbs Sampler 1000 times
        motifs = [random.choice([dna[i][j:j+k] for j in range(len(dna[i]) - k + 1)]) for i in range(t)]
        current_best_motifs = motifs[:]
        current_best_score = score(motifs, k)

        for _ in range(n):
            i = random.randint(0, t-1)
            temp_motifs = motifs[:i] + motifs[i+1:]
            profile = create_profile_with_pseudocounts(temp_motifs, k)
            motifs[i] = profile_randomly_generated_kmer(dna[i], k, profile)
            current_score = score(motifs, k)

            if current_score < current_best_score:
                current_best_motifs = motifs[:]
                current_best_score = current_score

        if current_best_score < best_score:
            best_motifs = current_best_motifs[:]
            best_scores = score_single_motif(best_motifs, k)
            best_score = current_best_score

    return best_motifs, best_scores

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = []
        sequence = ""
        for line in file:
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence.upper())
                    sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence.upper())
    return sequences

def generate_random_sequences(dna, num_sequences):
    random_sequences = []
    for _ in range(num_sequences):
        seq = ''.join(random.choice('ACGT') for _ in range(len(dna[0])))
        random_sequences.append(seq)
    return random_sequences

def calculate_p_value(observed_score, background_scores):
    count = sum(1 for score in background_scores if score <= observed_score)
    p_value = count / len(background_scores)
    return p_value

def main():
    parser = argparse.ArgumentParser(description="Gibbs Sampler Motif Finder")
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    parser.add_argument('-t', type=int, required=True, help='Number of sequences to be printed')
    parser.add_argument('-k', type=int, required=True, help='Length of motif to be found')
    parser.add_argument('-n', type=int, required=True, help='Number of iterations to be run')
    args = parser.parse_args()

    # Read DNA sequences
    dna_sequences = read_fasta(args.fasta)
    
    if len(dna_sequences) < args.t:
        raise ValueError(f"Expected at least {args.t} sequences in the FASTA file, but found {len(dna_sequences)}")
    dna_sequences = dna_sequences[:args.t]

    # Run Gibbs sampler 
    best_motifs, best_scores = gibbs_sampler(dna_sequences, args.k, len(dna_sequences), args.n)
    #print(best_scores)

    # Background score distribution
    num_random_sequences = 1000
    random_sequences = generate_random_sequences(dna_sequences, num_random_sequences)
    background_scores = [score([random.choice([random_sequences[i][j:j+args.k] for j in range(len(random_sequences[i]) - args.k + 1)]) for i in range(args.t)], args.k) for _ in range(num_random_sequences)]

    # Calculate p-values for each motif
    p_values = [calculate_p_value(s, background_scores) for s in best_scores]

    # Results
    print("Top motifs found:")
    for motif, p_value in zip(best_motifs, p_values):
        print(f"Motif: {motif}")

if __name__ == "__main__":
    main()
