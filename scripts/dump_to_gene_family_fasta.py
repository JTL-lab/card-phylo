#!/usr/bin/env python
import json
import os
import argparse
from Bio import SeqIO

class CARD():
    """
    Parser for the CARD database
    """
    def __init__(self, card_json_fp, rrna=False):
        with open(card_json_fp) as fh:
            with open(card_json_fp) as fh:
                self.card = json.load(fh)
            self.version = self.card['_version']

            # to avoid having to except them later when parsing the other
            # entries
            del self.card['_version']
            del self.card['_timestamp']
            del self.card['_comment']

            self.supported_models = ['protein homolog model',
                                     'protein variant model',
                                     'protein overexpression model']
            self.proteins, self.nucleotides = self.get_sequences()
            self.aro_to_gene_family = self.build_aro_to_gene_family()
            self.gene_family_to_aro = self.build_gene_family_to_aro()

    def build_aro_to_gene_family(self):
        aro_to_gene_family = {}
        for card_item in self.card.values():
            ARO_acc = card_item['ARO_accession']
            if card_item['model_type'] not in self.supported_models:
                continue

            # as multiple gene families are possible per ARO
            # although they are relatively rare so maybe I should talk with
            # Andrew about making them unique
            gene_families = []
            for category in card_item['ARO_category'].values():
                if category['category_aro_class_name'] == 'AMR Gene Family':
                    gene_families.append(category['category_aro_name'])

            # fix the situations where there are multiple gene families manually
            # all glycopeptide resistant gene clusters have 2 gene families one
            # indicating that it is a grgc and the other with its class
            # for now we are just using the cluster level and will deal with
            # the specifics at ARO level.
            if "glycopeptide resistance gene cluster" in gene_families:
                gene_families.remove('glycopeptide resistance gene cluster')

            # these all seem to be fusions so assigning a fusion class
            if 'AAC(3)' in gene_families and "AAC(6')" in gene_families:
                gene_families = ["AAC(3)_AAC(6')"]

            if "APH(2'')" in gene_families and "AAC(6')" in gene_families:
                gene_families = ["APH(2'')_AAC(6')"]

            if "ANT(3'')" in gene_families and "AAC(6')" in gene_families:
                gene_families = ["ANT(3'')_AAC(6')"]

            # also a fusion so assigned a new fusion class
            if 'class C LRA beta-lactamase' in gene_families and \
                    'class D LRA beta-lactamase' in gene_families:
                gene_families = ['class D/class C beta-lactamase fusion']

            ## efflux components
            if ARO_acc in ['300815','3000832', '3000676', '3000833', '3000263',
                            '3000817', '3003585', '3003381','3003511', '3003383',
                            '3003382', '3003896', '3003895', '3004107', '3000815',
                            '3000823']:
                gene_families = ['efflux regulator']
            if ARO_acc == '3000237':
                gene_families = ['efflux component']

            aro_to_gene_family.update({ARO_acc: gene_families})

        # tidy up and make unique aro:amr family relationship
        mapping_failure = False
        aros_without_protein = []
        for aro, gene_families in aro_to_gene_family.items():
            if len(gene_families) != 1:
                mapping_failure = True
                print(aro, gene_families)
            else:
                aro_to_gene_family[aro] = gene_families[0]

            if aro not in self.proteins:
                aros_without_protein.append(aro)

        if mapping_failure:
            raise ValueError("AROs and gene families don't map 1:1")

        # remove any aro without a protein sequence
        for aro in aros_without_protein:
            del aro_to_gene_family[aro]

        return aro_to_gene_family

    def build_gene_family_to_aro(self):
        gene_family_to_aro = {}
        for aro, gene_family in self.aro_to_gene_family.items():
            if gene_family not in gene_family_to_aro:
                gene_family_to_aro.update({gene_family: [aro]})
            else:
                gene_family_to_aro[gene_family].append(aro)

        return gene_family_to_aro

    def get_sequences(self):
        """
        Gather list of accession, prot and nucleotide sequence tuples from
        the card.json
        """
        data = {'protein_sequence': {}, 'dna_sequence': {}}
        for key, card_item in self.card.items():
            if card_item['model_type'] in self.supported_models:
                aro = card_item['ARO_accession']
                aro_name = card_item['ARO_name']
                sequences = card_item['model_sequences']['sequence']
                for seq_ix in sequences:
                    for sequence_type in sequences[seq_ix]:
                        if sequence_type in ['protein_sequence', 'dna_sequence']:
                            sequence = sequences[seq_ix][sequence_type]
                            if aro not in data[sequence_type]:
                                acc = ">gb|{}|{}|{}|".format(sequence['accession'],
                                                             aro,
                                                             aro_name.replace(' ', '_'))
                                data[sequence_type].update({aro: (acc,
                                                                  sequence['sequence'])})
        # what the fuck is happening why are we getting identical sequences
        # temporarily let's just add the first sequence
                            #else:
                            #    data[sequence_type][aro].append((acc,
                            #                                     sequence['sequence']))

        return data['protein_sequence'], data['dna_sequence']

    def write_seqs(self, seq_dict, seq_file_fp):
        if not os.path.exists(seq_file_fp):
            with open(seq_file_fp, 'w') as fh:
                for seq in seq_dict.values():
                    fh.write("{}\n{}\n".format(seq[0], seq[1]))

    def write_proteins(self, seq_file_fp):
        self.write_seqs(self.proteins, seq_file_fp)

    def write_nucleoties(self, seq_file_fp):
        self.write_seqs(self.nucleotides, seq_file_fp)

    def get_protein_per_family(self):
        """
        Write nucleotides sequences to per family fasta files
        """
        if not os.path.exists('family_fasta'):
            os.mkdir('family_fasta')

        for key, card_item in self.card.items():
            if card_item['model_type'] in self.supported_models:
                aro = card_item['ARO_accession']
                aro_name = card_item['ARO_name']
                sequences = card_item['model_sequences']['sequence']
                for seq_ix in sequences:
                    for sequence_type in sequences[seq_ix]:
                        if sequence_type == 'protein_sequence':
                            sequence = sequences[seq_ix][sequence_type]
                            with open(self.convert_amr_family_to_filename(self.aro_to_gene_family[aro]), 'a') as out_fh:
                                acc = ">gb|{}|{}|{}|".format(sequence['accession'],
                                                             aro,
                                                             aro_name.replace(' ', '_'))
                                seq = sequence['sequence']
                                out_fh.write(acc + '\n' + seq + '\n')

    def convert_amr_family_to_filename(self, family):
        fp = os.path.join('family_fasta', family.replace(' ', '_').replace('/', '_'))
        return fp

    def add_prevalence_to_family(self, prevalence_fasta):
        for record in SeqIO.parse(prevalence_fasta, "fasta"):
            try:
                aro = record.description.split('|')[2].replace('ARO:', '')
            except:
                print(record)
                assert False
            with open(self.convert_amr_family_to_filename(self.aro_to_gene_family[aro]), 'a') as out_fh:
                SeqIO.write(record, out_fh, "fasta")

def run():
    parser = argparse.ArgumentParser(description='Organise CARD sequences by annotated family.')
    parser.add_argument('-c', '--card_json', type=str, required=True,
                    help="Path to CARD canonical CARD.json file")
    parser.add_argument('-p', '--prevalence_fasta', type=str, required=True,
                    help="Path to CARD prevalence folder")
    parser.add_argument('-o', '--output_folder', type=str, required=True,
                    help="Folder to output family sequences")
    args = parser.parse_args()

    card = CARD(args.card_json)
    card.get_protein_per_family()
    card.add_prevalence_to_family(args.prevalence_fasta)

if __name__ == '__main__':
    run()
