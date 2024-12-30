import numpy as np
import subprocess
import os
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,  # Set to DEBUG for more verbose output
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),  # Log to console
    ],
)

def run_jellyfish(sequence, k, output_file="kmer_counts.jf"):
    """
    Run Jellyfish to count k-mers for a given sequence.
    The sequence is written to a temporary FASTA file for processing.
    """
    logging.info(f"Running Jellyfish for k={k} on a sequence of length {len(sequence)}.")

    # Create a temporary FASTA file
    temp_fasta = "temp_sequence.fasta"
    with open(temp_fasta, "w") as f:
        f.write(">seq\n")
        f.write(sequence + "\n")
    logging.debug(f"Temporary FASTA file created: {temp_fasta}")

    # Run Jellyfish to count k-mers
    jellyfish_command = [
        "jellyfish", "count", "-m", str(k), "-s", "100M", "-t", "8",
        "-C", temp_fasta, "-o", output_file
    ]
    subprocess.run(jellyfish_command, check=True)
    logging.info(f"Jellyfish count completed. Output: {output_file}")

    # Dump k-mer counts to a text file
    kmer_counts_file = "kmer_counts.txt"
    jellyfish_dump_command = [
        "jellyfish", "dump", output_file, "-o", kmer_counts_file
    ]
    subprocess.run(jellyfish_dump_command, check=True)
    logging.info(f"Jellyfish dump completed. K-mer counts saved to: {kmer_counts_file}")

    # Parse k-mer counts
    kmer_counts = {}
    with open(kmer_counts_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # Line with the count
                try:
                    count = int(line[1:])  # Parse the count after ">"
                except ValueError:
                    logging.error(f"Invalid count value: {line}")
                    raise ValueError(f"Invalid count value: {line}")
            else:  # Line with the k-mer
                kmer = line
                kmer_counts[kmer] = count
                logging.debug(f"K-mer: {kmer}, Count: {count}")

    # Clean up temporary files
    os.remove(temp_fasta)
    os.remove(output_file)
    os.remove(kmer_counts_file)
    logging.info("Temporary files cleaned up.")

    return kmer_counts


def generate_kmer_matrix(sequence, k, window_size, step):
    """
    Generate k-mer frequency vectors for sliding windows using Jellyfish.
    """
    logging.info(f"Generating k-mer matrix for sequence of length {len(sequence)}.")
    seq_length = len(sequence)
    vectors = []

    for start in range(0, seq_length - window_size + 1, step):
        window_seq = sequence[start:start + window_size]
        logging.debug(f"Processing window: Start={start}, Window size={window_size}")
        kmer_counts = run_jellyfish(window_seq, k)

        # Convert k-mer counts into a fixed-size vector
        from itertools import product
        nucleotides = ['A', 'C', 'G', 'T']
        all_kmers = [''.join(x) for x in product(nucleotides, repeat=k)]
        vec = np.zeros(len(all_kmers), dtype=float)
        for i, kmer in enumerate(all_kmers):
            vec[i] = kmer_counts.get(kmer, 0)
        vectors.append(vec)

    logging.info(f"K-mer matrix generation completed. Number of windows: {len(vectors)}")
    return np.array(vectors)


def normalized_cross_correlogram(freq_matrix, max_lag=20):
    """
    Compute normalized cross-correlogram features for a frequency matrix.
    """
    logging.info("Computing cross-correlogram features.")
    num_windows = freq_matrix.shape[0]
    correlogram = []

    for lag in range(1, max_lag + 1):
        corrs = []
        for i in range(num_windows - lag):
            vec1 = freq_matrix[i]
            vec2 = freq_matrix[i + lag]

            # Normalize the vectors
            norm_vec1 = vec1 / (np.linalg.norm(vec1) + 1e-6)
            norm_vec2 = vec2 / (np.linalg.norm(vec2) + 1e-6)

            # Pearson correlation
            if np.std(norm_vec1) == 0 or np.std(norm_vec2) == 0:
                corr = 0
            else:
                corr = np.corrcoef(norm_vec1, norm_vec2)[0, 1]
            corrs.append(corr)

        # Mean correlation for this lag
        correlogram.append(np.mean(corrs) if corrs else 0.0)

    logging.info("Cross-correlogram computation completed.")
    return np.array(correlogram)  # shape: (max_lag,)


def compute_4cac_features(sequence, k=4, window_size=200, step=200, max_lag=20):
    """
    Compute the cross-correlogram features for a single sequence
    based on the 4CAC methodology.
    """
    logging.info("Starting computation of 4CAC features.")
    freq_matrix = generate_kmer_matrix(sequence, k, window_size, step)
    if freq_matrix.shape[0] < 2:
        logging.warning("Fewer than 2 windows. Returning zero vector.")
        return np.zeros(max_lag)

    correlogram_feat = normalized_cross_correlogram(freq_matrix, max_lag)
    logging.info("4CAC feature computation completed.")
    return correlogram_feat
