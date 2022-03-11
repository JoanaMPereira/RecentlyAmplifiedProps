import subprocess as sp

import json
import os
import sys

# define inputs
indomain_annotation = sys.argv[1]
ecod_json = sys.argv[2]
selected_pdbs_folder = sys.argv[3]

# define working directory
working_dir = os.getcwd()

# HELPING ROUTINES

def get_ecod_domain_annotations(indomain_annotation):

    ecod_mapping = {}
    indomain_annotation = json.load(open(indomain_annotation, 'r'))

    for subject in indomain_annotation.keys():
        f_name = 'UNKNOWN'
        t_name = 'UNKNOWN'
        domain_name = 'UNKNOWN'
        domain_annotated = indomain_annotation[subject]['annotated_domains']
        prob = indomain_annotation[subject]['annotated_probability']

        if len(domain_annotated) > 0:
            # print(domain_annotated)
            for i, domain in enumerate(domain_annotated):
                if prob[i] > 50:
                    if 'propeller' in domain:
                        f_name = domain.split()[-1]
                        t_name = domain.split('T: ')[1].split(' F: ')[0]
                        domain_name = domain
                        if 'Partial' in domain_name:
                            domain_name = domain_name.split(') ')[-1]
                    else:
                        f_name = 'NO_PROPELLER'
                        t_name = 'NO_PROPELLER'
                        domain_name = 'NO_PROPELLER'
                else:
                    prob[i] = []


        if domain_name not in ecod_mapping:
            ecod_mapping[domain_name] = {'f': f_name, 'matched_propellers': [], 'match_probability': [], 't': t_name, 'number_blades': []}

        if f_name == 'NO_PROPELLER':
            print(subject, domain_annotated, indomain_annotation[subject]['annotated_probability'])

        ecod_mapping[domain_name]['matched_propellers'].append(subject)
        ecod_mapping[domain_name]['match_probability'].append(prob)

    return ecod_mapping

def add_t_names_and_blades(ecod_mapping, ecod_json, selected_pdbs_folder):

    ecod_data = json.load(open(ecod_json, 'r'))

    for pdb in ecod_data.keys():
        t_name = ecod_data[pdb]['t']
        f_name = ecod_data[pdb]['f']

        for domain_name in ecod_mapping:
            if ecod_mapping[domain_name]['f'] == f_name and ecod_mapping[domain_name]['t'] == t_name:

                blades_seq_file = '{}/{}_blades.fasta'.format(selected_pdbs_folder, pdb)

                if os.path.isfile(blades_seq_file):
                    blades = get_blades_from_seq_file(blades_seq_file)
                    ecod_mapping[domain_name]['number_blades'].append(len(blades))


    return ecod_mapping

def get_blades_from_seq_file(blades_seq_file):

    blades = {}

    with open(blades_seq_file, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                header = line.strip().replace('>', '')
            else:
                seq = line.strip()
                blades[header] = seq

    return blades

# MAIN CODE

ecod_mapping = get_ecod_domain_annotations(indomain_annotation)
ecod_mapping = add_t_names_and_blades(ecod_mapping, ecod_json, selected_pdbs_folder)

print(len(ecod_mapping))
json.dump(ecod_mapping, open('{}_ecod_mapping.json'.format(indomain_annotation.split('/')[-1][:-6]), 'w'), indent = 4)
