from feature_extraction import compute_4cac_features
from data_ingestion import stream_fasta_files

def main():
    data_dir = '../../data/raw/genomes'

    import csv
    output_csv = './4cac_features_with_jellyfish.csv'
    k = 4
    window_size = 200
    step = 200
    max_lag = 20

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['sequence_id'] + [f"lag_{lag}" for lag in range(1, max_lag + 1)]
        writer.writerow(header)

        seq_count = 0

        for file_path, seq_generator in stream_fasta_files(data_dir):
            for rec in seq_generator:
                seq_str = str(rec.seq)
                seq_id = f"{file_path}::{rec.id}"

                # Compute cross-correlogram features using Jellyfish
                features = compute_4cac_features(seq_str, k, window_size, step, max_lag)
                row = [seq_id] + features.tolist()
                writer.writerow(row)

                seq_count += 1

    print(f"Done. Processed {seq_count} sequences. Features written to {output_csv}")


if __name__ == '__main__':
    main()
