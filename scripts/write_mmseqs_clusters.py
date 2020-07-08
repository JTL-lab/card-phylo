#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO

def run():
    parser = argparse.ArgumentParser(description='Write mmseq sequence cluster files and dump singeltons.')
    parser.add_argument('-c', '--cluster_fasta', type=str, required=True,
                    help="Path to MMSEQS2 easy-cluster fasta")
    parser.add_argument('-o', '--output_folder', type=str, required=True,
                    help="Folder to output family sequences")
    args = parser.parse_args()

    if not os.path.exists(args.output_folder):
        os.mkdir(args.output_folder)

    cluster_seqs = []
    out_fh = None
    cluster = 0
    num_seqs = 0

    singleton_record = open(os.path.join(args.output_folder, 'phylo_singletons.txt'), 'w')
    for record in SeqIO.parse(args.cluster_fasta, 'fasta'):
        if len(record.seq) == 0:
            if out_fh:
                if num_seqs >= 3:
                    SeqIO.write(cluster_seqs, out_fh, 'fasta')
                    cluster += 1
                else:
                    SeqIO.write(cluster_seqs, singleton_record, "fasta")
                    print(f"MMSEQS2 phylo singleton (<3 seqs): {cluster_seqs[0].id}")
                cluster_seqs = []
                num_seqs = 0
                out_fh.close()
            out_fh = open(os.path.join(args.output_folder, str(cluster) + ".faa"), 'w')
        else:
            cluster_seqs.append(record)
            num_seqs += 1
    singleton_record.close()

if __name__ == '__main__':
    run()
