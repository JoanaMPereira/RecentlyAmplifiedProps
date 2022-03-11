import multiprocessing as mp
import subprocess as sp
import numpy as np

import argparse
import time
import os
import json
import ast

in_folder = 'selected_cluster_representatives_after_curation'

## HELPING ROUTINES

def parse_fasta(in_fasta):

    sequences = {}

    with open(in_fasta, 'r') as infasta:
        for line in infasta:
            if line.startswith('>'):
                ncbi_code = line[1:].split('|')[0]
            else:
                sequences[ncbi_code] = line.strip()
                
    return sequences

def update_sequences(sequences, curr_sequences):

    for ncbi_code in curr_sequences:
        if ncbi_code not in sequences:
            sequences[ncbi_code] = curr_sequences[ncbi_code]

    return sequences

def update_taxonomy(taxonomy, target_taxonomy):

    curr_taxonomy = json.load(open(target_taxonomy, 'r'))

    for superkingdom in curr_taxonomy:
        if superkingdom not in taxonomy:
            taxonomy[superkingdom] = {}

        for phylum in curr_taxonomy[superkingdom]:
            if phylum not in taxonomy[superkingdom]:
                taxonomy[superkingdom][phylum] = {}

            for taxclass in curr_taxonomy[superkingdom][phylum]:
                if taxclass not in taxonomy[superkingdom][phylum]:
                    taxonomy[superkingdom][phylum][taxclass] = {}

                for order in curr_taxonomy[superkingdom][phylum][taxclass]:
                    if order not in taxonomy[superkingdom][phylum][taxclass]:
                        taxonomy[superkingdom][phylum][taxclass][order] = {}

                    for genus in curr_taxonomy[superkingdom][phylum][taxclass][order]:
                        for species in curr_taxonomy[superkingdom][phylum][taxclass][order][genus]:
                            correct_genus = species.split()[0]
                            
                            if correct_genus not in taxonomy[superkingdom][phylum][taxclass][order]:
                                taxonomy[superkingdom][phylum][taxclass][order][correct_genus] = {}
                                
                            if species not in taxonomy[superkingdom][phylum][taxclass][order][correct_genus]:
                                taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species] = {'Medians': [], 'ncbi_codes': [], 'Propellers': [], 'Intervals': []}

                            for i, ncbi_code in enumerate(curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes']):
                                curr_propellers = ast.literal_eval(curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Propellers'][i])
                                curr_medians = ast.literal_eval(curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Medians'][i])
                                curr_intervals = ast.literal_eval(curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Intervals'][i])
                                
                                if ncbi_code not in taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['ncbi_codes']:
                                    taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['ncbi_codes'].append(ncbi_code)
                                    taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['Propellers'].append(curr_propellers)
                                    taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['Medians'].append(curr_medians)
                                    taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['Intervals'].append(curr_intervals)
                                else:
                                    ncbi_code_index = taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['ncbi_codes'].index(ncbi_code)
                                    old_intervals = taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['Intervals'][ncbi_code_index]
                                    
                                    curr_coverage = int(curr_intervals[-1][-1]) - int(curr_intervals[0][0])
                                    old_coverage = int(old_intervals[-1][-1]) - int(old_intervals[0][0])

                                    if curr_coverage > old_coverage:
                                        taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['Propellers'][ncbi_code_index] = curr_propellers
                                        taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['Intervals'][ncbi_code_index] = curr_intervals
                                        taxonomy[superkingdom][phylum][taxclass][order][correct_genus][species]['Medians'][ncbi_code_index] = curr_medians
  
    return taxonomy

def collect_all_sequences(in_folder):

    searches = ['{}/{}'.format(in_folder, folder) for folder in os.listdir(in_folder) if os.path.isdir('{}/{}'.format(in_folder, folder))]

    sequences = {}
    taxonomy = {}
    
    for search in searches:
        matches_found = 0
        databases = ['{}/{}'.format(search, db) for db in os.listdir(search) if os.path.isdir('{}/{}'.format(search, db))]

        for database in databases:
            rounds = ['{}/{}'.format(database, rnds) for rnds in os.listdir(database) if os.path.isdir('{}/{}'.format(database, rnds))]

            for rnd in rounds:
                target_fasta = '{}/selection_of_relevant_sequences_{}/selected_relevant_sequences_full_length.fasta_taxonomy'.format(rnd, database.split('/')[-1])
                target_taxonomy = '{}/selection_of_relevant_sequences_{}/selected_relevant_sequences_full_length.fasta_taxonomy.json'.format(rnd, database.split('/')[-1])
                if os.path.isfile(target_fasta):
                    matches_found += 1
                    curr_sequences = parse_fasta(target_fasta)
                    sequences = update_sequences(sequences, curr_sequences)

                    taxonomy = update_taxonomy(taxonomy, target_taxonomy)
##
##        if matches_found > 0:
##            print(search)

    return sequences, taxonomy

def write_to_file(all_sequences, curr_taxonomy, in_folder):

    with open('{}/all_searches_concatenated_taxonomy.json'.format(in_folder), 'w') as fp:
        json.dump(all_taxonomy, fp, indent = 4)
    with open('{}/all_searches_concatenated_sequences.json'.format(in_folder), 'w') as fp:
        json.dump(all_sequences, fp, indent = 4)

    out_fasta = '{}/all_searches_concatenated_full_length_sequences.fasta'.format(in_folder)
    with open(out_fasta, 'w') as fp:
        for superkingdom in curr_taxonomy:    
            for phylum in curr_taxonomy[superkingdom]:
                for taxclass in curr_taxonomy[superkingdom][phylum]:
                    for order in curr_taxonomy[superkingdom][phylum][taxclass]:
                        for genus in curr_taxonomy[superkingdom][phylum][taxclass][order]:
                            for species in curr_taxonomy[superkingdom][phylum][taxclass][order][genus]:
                                for i, ncbi_code in enumerate(curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes']):
                                    propellers = curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Propellers'][i]
                                    medians = curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Medians'][i]
                                    intervals = curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Intervals'][i]

                                    protID = '>{}|{}|{}|{}|{}|{}|{}| Intervals: {} | Medians: {} | Propellers: {}'.format(ncbi_code, superkingdom, phylum, taxclass, order, genus, species, intervals, medians, propellers)
                                    fp.write('{}\n{}\n'.format(protID, all_sequences[ncbi_code]))

    out_fasta = '{}/all_searches_concatenated_propellers_sequences.fasta'.format(in_folder)
    with open(out_fasta, 'w') as fp:
        for superkingdom in curr_taxonomy:    
            for phylum in curr_taxonomy[superkingdom]:
                for taxclass in curr_taxonomy[superkingdom][phylum]:
                    for order in curr_taxonomy[superkingdom][phylum][taxclass]:
                        for genus in curr_taxonomy[superkingdom][phylum][taxclass][order]:
                            for species in curr_taxonomy[superkingdom][phylum][taxclass][order][genus]:
                                for i, ncbi_code in enumerate(curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes']):
                                    propellers = curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Propellers'][i]
                                    medians = curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Medians'][i]
                                    intervals = curr_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Intervals'][i]
                                    full_sequence = all_sequences[ncbi_code]
                                    for j, propeller in enumerate(propellers):
                                        curr_intervals = [intervals[a] for a in propeller]
                                        curr_interval = [curr_intervals[0][0], curr_intervals[-1][-1]]
                                        curr_sequence = full_sequence[curr_interval[0]-1 : curr_interval[-1]]

                                        protID = '>{}_#{}({} blades: {}%)|{}|{}|{}|{}|{}|{}'.format(ncbi_code, j+1, len(propeller), medians[j], superkingdom, phylum, taxclass, order, genus, species)
                                        fp.write('{}\n{}\n'.format(protID, curr_sequence))

# To get single blade propellers

def parse_full_length_fasta(target_full_length_fasta):

    full_sequences = {}

    with open(target_full_length_fasta, 'r') as infasta:
        for line in infasta:
            if line.startswith('>'):
                ncbi_code = line[1:].split('|')[0]
            else:
                full_sequences[ncbi_code] = line.strip()

    return full_sequences

def get_significant_single_propellers(curr_propellers, max_evalue, max_len, target_full_length_fasta):

    all_full_lengths = parse_full_length_fasta(target_full_length_fasta)
    single_bladed = {}

    for propeller in curr_propellers:
        if len(curr_propellers[propeller]['Evalues']) == 1 and curr_propellers[propeller]['Evalues'][0] <= max_evalue and 'partial' not in curr_propellers[propeller]['Title']:
            if curr_propellers[propeller]['Intervals'][0][-1] - curr_propellers[propeller]['Intervals'][0][0] <= max_len:
                single_bladed[propeller] = curr_propellers[propeller]
                single_bladed[propeller]['full_length'] = all_full_lengths[propeller]

    return single_bladed

def find_single_blade_propellers(in_folder, max_evalue = 1e-10, max_len = 40):

    searches = ['{}/{}'.format(in_folder, folder) for folder in os.listdir(in_folder) if os.path.isdir('{}/{}'.format(in_folder, folder))]

    sequences = {}
    
    for search in searches:
        matches_found = 0
        search_label = '_'.join(search.split('/')[-1].split('_')[:2])
        databases = ['{}/{}'.format(search, db) for db in os.listdir(search) if os.path.isdir('{}/{}'.format(search, db))]

        for database in databases:
            database_name = database.split('/')[-1]
            rounds = ['{}/{}'.format(database, rnds) for rnds in os.listdir(database) if os.path.isdir('{}/{}'.format(database, rnds))]

            for rnd in rounds:
                target_fasta = '{}/{}_representative_blades_psiblast_15collected_rounds_{}_enriched_set_full_sequence_max_eval_1.fasta'.format(rnd, search_label, database_name)
                target_blades = '{}/missing_blade_searching_{}/sequences_with_missed_terminal_blades.json'.format(rnd, database_name)
                if os.path.isfile(target_blades):
                    matches_found += 1
                    curr_propellers = json.load(open(target_blades, 'r'))
                    single_propellers = get_significant_single_propellers(curr_propellers, max_evalue, max_len, target_fasta)
                    
                    if len(single_propellers) > 0:
                        for propeller in single_propellers:
                            single_propellers[propeller]['family'] = search.split('/')[-1]
                        sequences = update_sequences(sequences, single_propellers)
                        
    return sequences

def write_single_propellers_to_file(single_propellers, in_folder):
    
    with open('{}/all_searches_concatenated_single_blade_propellers.json'.format(in_folder), 'w') as fp:
        json.dump(single_propellers, fp, indent = 4)

    out_fasta = '{}/all_searches_concatenated_single_blade_propellers_annotated.fasta'.format(in_folder)
    with open(out_fasta, 'w') as fp:
        for propeller in single_propellers:
            propeller_interval = single_propellers[propeller]['Intervals'][0]
            identifier = '>{}| {} | Evalue: {} | Interval: {} | ECOD family: {}'.format(propeller, single_propellers[propeller]['Title'], single_propellers[propeller]['Evalues'][0], propeller_interval, single_propellers[propeller]['family'])
            full_length_seq = single_propellers[propeller]['full_length']
            fp.write('{}\n{}\n{} ### The blade\n{}\n\n'.format(identifier, single_propellers[propeller]['full_length'][:propeller_interval[0]-1], single_propellers[propeller]['full_length'][propeller_interval[0]-1:propeller_interval[-1]], single_propellers[propeller]['full_length'][propeller_interval[-1]:]))

    out_fasta = '{}/all_searches_concatenated_single_blade_propellers.fasta'.format(in_folder)
    with open(out_fasta, 'w') as fp:
        for propeller in single_propellers:
            propeller_interval = single_propellers[propeller]['Intervals'][0]
            identifier = '>{}| {} | Propellers: [[0]] | Evalue: {} | Intervals: [{}] | ECOD family: {}'.format(propeller, single_propellers[propeller]['Title'], single_propellers[propeller]['Evalues'][0], propeller_interval, single_propellers[propeller]['family'])
            full_length_seq = single_propellers[propeller]['full_length']
            fp.write('{}\n{}\n'.format(identifier, single_propellers[propeller]['full_length']))


## MAIN CODE

all_sequences, all_taxonomy = collect_all_sequences(in_folder)
write_to_file(all_sequences, all_taxonomy, in_folder)

single_propellers = find_single_blade_propellers(in_folder, max_evalue = 1e-10, max_len = 45)
write_single_propellers_to_file(single_propellers, in_folder)
