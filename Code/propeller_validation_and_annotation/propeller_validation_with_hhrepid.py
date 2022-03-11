
import multiprocessing as mp
import subprocess as sp
import numpy as np
import subprocess as sp
import pandas as pd

import argparse
import time
import os
import ast
import json
import string

from scipy import stats

# Define tools locations
muscle = '/Applications/muscle/muscle'
hhrepid = 'hhrepid'
hhrepid_cal_hmm = '/ebio/abt1_share/toolkit_sync/bioprogs/tools/hhrepid/cal_small.hhm'

# Get inputs
parser = argparse.ArgumentParser(description='')
requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments with defaults')

# required inputs
requiredNamed.add_argument('-i', dest='in_fasta', type=str, required=True, help='input fasta of full length sequences with the preliminary propellers')
requiredNamed.add_argument('-taxonomy', dest='taxonomy', type=str, required=True, help='input json file with taxonomy')
# optional inputs
optionalNamed.add_argument('-min_seqID', dest='min_seqID', default = 60,type=int, help='Minimum median sequence identity of blades to consider (default: 60)')
optionalNamed.add_argument('-min_pval', dest='min_pval', default = 1e-1,type=int, help='Minimum pvalue of a repeat unit to be considered (default: 1e-1)')
optionalNamed.add_argument('-cpu', dest='n_cpu', default = 2,type=int, help='Number of cpus to use (default: 2)')
optionalNamed.add_argument('-threads', dest='n_threads', default=1, type=int, help='Number of parallel jobs (default: 1)')
optionalNamed.add_argument('-force_search', dest='force_search', default='False', type=str, help='Boolean statement to force search even if it was already preformed (default: False)')
optionalNamed.add_argument('-tmp_folder', dest='tmp_folder', default = '/tmp/',type=str, help='Temporary folder (default: /tmp/)')

args = parser.parse_args()
curr_directory = os.getcwd()

# PRINT INPUTS SUMMARY
hhsearch_cpu = int(args.n_cpu/args.n_threads)

print('\nINPUT PARAMETERS:')
print(' ... Min. seq. identity: {}'.format(args.min_seqID))
print(' ...   Min. rep. pvalue: {}'.format(args.min_pval))
print(' ...  Number of threads: {}'.format(args.n_threads))
print(' ...  CPUs for HHsearch: {}\n'.format(hhsearch_cpu))

## HELPING ROUTINES

# 0. Helping math

def median_absolute_deviation(x):
    return np.median(np.abs(x-np.median(x)))

# 1. Job preparation

def parse_fasta(in_fasta):

    sequences = {}
    annotated_propellers = 0

    with open(in_fasta, 'r') as infst:
        for line in infst:
            if len(line) > 0:
                line = line.strip()
                if line.startswith('>') or '|' in line:
                    ncbi_code = line.split('|')[0].replace('>','')
                    sequences[ncbi_code] = {'seq': '', 'annotated_domains': [], 'annotated_intervals': [], 'annotated_blades': []}

                    if 'Propellers:' in line:
                        intervals = ast.literal_eval(line.split('Intervals:')[-1].split('|')[0].strip())
                        propellers = ast.literal_eval(line.split('Propellers:')[-1].split('|')[0].strip())
                        
                        for i, propeller in enumerate(propellers):

                            blades_intervals = [intervals[blade] for blade in propeller]

                            propeller_interval = [intervals[propeller[0]][0], intervals[propeller[-1]][-1]]
                            propeller_label = 'beta-propeller {}-blades'.format(len(propeller))

                            sequences[ncbi_code]['annotated_domains'].append(propeller_label)
                            sequences[ncbi_code]['annotated_intervals'].append(propeller_interval)
                            annotated_propellers += 1
                    
                            sequences[ncbi_code]['annotated_blades'].append(blades_intervals)

                else:
                    sequences[ncbi_code]['seq'] = line.strip()

    print("INPUT SEQUENCES: Found input {} sequences, totalling {} already annotated propellers".format(len(sequences), annotated_propellers))
    
    return sequences

# 2. Repeat annotation

def chunk_list(l, n):
    chunks = np.array_split(np.array(l), n)
    chunks = [list(chunk) for chunk in chunks]
    return chunks

def annotate_repeats_in_sequence(ncbi_code_data, tmp_label, add_types = True):

    #print(' ... ... ... Running HHrepID')
    # save sequence to file
    tmp_fasta = '{}{}_sequence.fasta'.format(args.tmp_folder, tmp_label)
    with open(tmp_fasta, 'w') as fp:
        fp.write('>{}\n{}\n'.format(tmp_label, ncbi_code_data['seq']))

    # run hhrepid for this file
    out_hhr = '{}.hhr'.format(tmp_fasta)
    run_hhrepid = sp.Popen([hhrepid, '-i', tmp_fasta, '-d', hhrepid_cal_hmm, '-plot', '0', '-o', out_hhr], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_hhrepid.communicate()

    if not os.path.isfile(out_hhr):
        repeats_found = []
        repeats_types = []

    else:
        # get repeats found
        repeats_found, repeats_types = parse_hhrepid_out(out_hhr, add_types = add_types)
        os.remove(out_hhr)

    ncbi_code_data['repeats_found'] = repeats_found
    if add_types:
        ncbi_code_data['repeats_types'] = repeats_types

    # remove temporary files
    os.remove(tmp_fasta)

    return ncbi_code_data

def parse_hhrepid_out(hhrepid_out, add_types = True, min_pval = args.min_pval):

    #print(' ... ... ... Parsing HHrepID')

    repeats_found = []
    repeats_types = []
    repeats_pvals = []

    found_table = False
    found_results = False
    with open(hhrepid_out, 'r') as hhrep:
        for line in hhrep:
            if line.startswith('ID') and 'P-value' in line:
                found_table = True
            elif 'Results for repeats' in line:
                found_results = True
            elif found_table and found_results and len(line.split()) == 8:
                interval = line.split()[6]
                interval = interval.split('-')
                interval = [int(i) for i in interval]
                pval = ast.literal_eval(line.split()[2])

                if interval not in repeats_found and pval <= min_pval:
                    repeats_found.append(interval)
                    repeats_pvals.append(pval)

                    if add_types:
                        repeats_types.append(line[0])

    if add_types:

        # if there are more than 1 type of repeats, sort and rename so that continuous repeated chunks of the same type 
        # have the same name. this alows to separate repeats of the same type that are separated by repeats of another type
        repeat_types_df = {'types': repeats_types, 'intervals': repeats_found, 'pval': repeats_pvals}
        repeat_types_df = pd.DataFrame(repeat_types_df)
        repeat_types_df = repeat_types_df.sort_values(by='intervals')
        repeat_types_df.index = range(len(repeat_types_df))

        new_letter_count = 0
        new_repeats_types = []
        new_repeats_found = []
        letters = list(string.ascii_uppercase)

        for index, row in repeat_types_df.iterrows():
            new_repeats_found.append(row.intervals)
            if len(new_repeats_types) == 0: 
                curr_letter = letters[0]
            elif row.types != repeat_types_df.types[index-1]:
                new_letter_count += 1
                curr_letter = letters[new_letter_count]
            new_repeats_types.append(curr_letter)

        repeat_types_df['new'] = new_repeats_types

        # now for each type, check whether there are clear long linkers and break the type into 2 if so
        # done by calculating the Z score of the linker length in the overall distribution of linker lengths
        # a linker is too large if it is 2 standard deviations larger than the average linker length 
        # AND if the block composed by the linker and the flanking intervals is also 2 standard deviations
        # larger than the average. This allows for large linkers that may contain pieces of a degenerated blade to 
        # not be considered too large and to not prematurately break a propeller

        linkers_lengths_to_previous = [0 for interval in repeat_types_df.intervals]
        pairs_lengths = [0 for interval in repeat_types_df.intervals]

        for index, row in repeat_types_df.iterrows():
            if index > 0:
                linker_length = row.intervals[0] - repeat_types_df.intervals[index-1][-1]
                linkers_lengths_to_previous[index] = linker_length

                pair_length = row.intervals[-1] - repeat_types_df.intervals[index-1][0]
                pairs_lengths[index] = pair_length

            else:
                pairs_lengths[index] = 2*(row.intervals[-1]-row.intervals[0])

        repeat_types_df['linker_lengths'] = linkers_lengths_to_previous
        repeat_types_df['linker_lengths_dev'] = stats.zscore(linkers_lengths_to_previous)

        repeat_types_df['pairs_lengths'] = pairs_lengths
        pairs_lengths_dev = []
        for repeat_type in list(set(repeat_types_df.new)):
            curr_pairs_lengths = repeat_types_df.loc[repeat_types_df.new == repeat_type]['pairs_lengths']
            pairs_lengths_dev += list(stats.zscore(curr_pairs_lengths))
        repeat_types_df['pairs_lengths_dev'] = pairs_lengths_dev
        
        new_letter_count = 0
        new_repeats_types = []

        for index, row in repeat_types_df.iterrows():
            if len(new_repeats_types) == 0: 
                curr_letter = letters[0]
            elif (row.pairs_lengths_dev > 2.0 and row.linker_lengths_dev > 2.0) or row.new != repeat_types_df.new[index-1]:
            # if the linker length is too large or we are now on another type, change the type
                new_letter_count += 1
                curr_letter = letters[new_letter_count]
            new_repeats_types.append(curr_letter)

        repeat_types_df['new_new'] = new_repeats_types
        
        # remove 'singleton' types (those that correspond to only one repeat)
        repeat_types_df = repeat_types_df[repeat_types_df.groupby('new_new').new.transform(len) > 1]

        # update lists
        repeats_found = list(repeat_types_df.intervals)
        repeats_types = list(repeat_types_df.new_new)

    return repeats_found, repeats_types

def update_propellers_intervals(ncbi_code_data):

    #print(' ... ... Updating propellers intervals')
    previous_propellers = ncbi_code_data['annotated_intervals']
    previous_propellers_labels = ncbi_code_data['annotated_domains']

    all_repeats_found = ncbi_code_data['repeats_found']
    repeats_types = ncbi_code_data['repeats_types']

    previous_coverage = sum([interval[-1]-interval[0] for interval in previous_propellers])
    current_coverage = sum([interval[-1]-interval[0] for interval in all_repeats_found])

    # if the repeats found cover more of the sequence than the previously found ones
    if len(all_repeats_found) > 0 and current_coverage > previous_coverage:

        # find to which repeat types our previous propellers correspond to
        curr_repeat_types_intervals = {repeat_type: [] for repeat_type in set(repeats_types)}
        propelles_types = {str(propeller): {'repeat_type': '', 'repeat_type_overlap': 0, 'corresponding_repeats': []} for propeller in previous_propellers}

        for repeat_type in set(repeats_types):
            curr_intervals = [interval for i, interval in enumerate(all_repeats_found) if repeats_types[i] == repeat_type]
            curr_repeat_types_intervals[repeat_type] = [curr_intervals[0][0], curr_intervals[-1][-1]]

            for index, propeller in enumerate(previous_propellers):
                # check how much the current propeller overlaps with the current repeat type and assign it the type with the largest overlap
                overlap_to_repeat = set(range(propeller[0], propeller[1])).intersection(set(range(curr_repeat_types_intervals[repeat_type][0], curr_repeat_types_intervals[repeat_type][1])))
                
                if len(overlap_to_repeat) > propelles_types[str(propeller)]['repeat_type_overlap']:

                    propelles_types[str(propeller)]['repeat_type'] = repeat_type
                    propelles_types[str(propeller)]['repeat_type_overlap'] = len(overlap_to_repeat)
                    propelles_types[str(propeller)]['corresponding_repeats'] = curr_intervals
                    propelles_types[str(propeller)]['name'] = previous_propellers_labels[index]

        # select only those repeat types that overlap with propellers
        types_with_propellers = list(set([propelles_types[propeller]['repeat_type'] for propeller in propelles_types if propelles_types[propeller]['repeat_type'] != '']))

        if len(types_with_propellers) > 0:

            new_propeller_intervals = []

            for repeat_type in types_with_propellers:

                repeats_found = [interval for i, interval in enumerate(all_repeats_found) if repeats_types[i] == repeat_type]
                
                # get the propellers in this repeat type. It could be that it either contains multiple propellers or
                # the propellers within it are just fragments of a bigger propeller
                # to figure that out, calculate the average blade length in the previously identified propellers and
                # check whether the lengths of the repeats found comprise roughly 1 blade or more
                previous_propellers_in_type = [propeller for propeller in previous_propellers if propelles_types[str(propeller)]['repeat_type'] == repeat_type]

                # get the general blade lengths for each propeller previously identified
                general_blade_length = np.median([(propeller[-1]-propeller[-0])/int(propelles_types[str(propeller)]['name'].split()[-1].split('-')[0]) for propeller in previous_propellers_in_type])

                # check if the detected repeats correspond to a full propeller or to individual blades 
                # (in cases where there are multiple propellers, it is possible that it identifies repeats as the propeller and not the blades)
                median_blades_per_repeat = np.median([round((repeat[-1]-repeat[0])/general_blade_length) for i, repeat in enumerate(repeats_found)])

                if int(median_blades_per_repeat) <= 1.0: # each repeat could be a blade. 
                    # define the chunck as a propeller
                    new_interval = [repeats_found[0][0], repeats_found[-1][-1]]
                    new_propeller_intervals.append(new_interval)
                           
                else: 
                    # we are dealing with a repeat unit corresponding to multiple blades. this is likely a repeat of multiple propellers. 
                    # each repeat will be an individual propeller
                    new_propeller_intervals += repeats_found  

                # add to the intervals any propeller that did not match any repeat type
                for propeller in propelles_types:
                    if propelles_types[propeller]['repeat_type'] == '':
                        new_propeller_intervals.append(ast.literal_eval(propeller))
      
        else:
            new_propeller_intervals = previous_propellers

    else:
        new_propeller_intervals = previous_propellers

    ncbi_code_data['updated_intervals'] = new_propeller_intervals
    #print(new_propeller_intervals)

    return ncbi_code_data

def update_propellers_annotations(ncbi_code_data, tmp_label):

    #print(' ... ... Updating propellers annotations')

    previous_annotations = ncbi_code_data['annotated_domains']
    new_domain_annotations = []

    previous_blades_intervals = ncbi_code_data['annotated_blades']
    new_blades_intervals = []

    previous_propellers = ncbi_code_data['annotated_intervals']

    protein_sequence = ncbi_code_data['seq']

    for i, updated_interval in enumerate(ncbi_code_data['updated_intervals']):

        try:
            chunk_sequence = protein_sequence[updated_interval[0]-1 : updated_interval[1]]

            chunk_data = ncbi_code_data.copy()
            chunk_data['seq'] = chunk_sequence

            # find repeats in this sequence chunk
            chunk_data = annotate_repeats_in_sequence(chunk_data, tmp_label = '{}_propeller{}'.format(tmp_label, i), add_types=True)

            # update the numberings
            chunk_data['repeats_found'] = [[x+updated_interval[0]-1 for x in interval] for interval in chunk_data['repeats_found']]
            all_repeats_found = chunk_data['repeats_found']

            # check if there are multiple types of repeats and which ones are in regions with previously detected propellers
            repeats_types = chunk_data['repeats_types']
            curr_repeat_types_intervals = {repeat_type: [] for repeat_type in set(repeats_types)}
            propelles_types = {str(propeller): {'repeat_type': '', 'repeat_type_overlap': 0, 'corresponding_repeats': []} for propeller in previous_propellers}

            for repeat_type in set(repeats_types):

                curr_intervals = [interval for x, interval in enumerate(all_repeats_found) if repeats_types[x] == repeat_type]
                curr_repeat_types_intervals[repeat_type] = [curr_intervals[0][0], curr_intervals[-1][-1]]

                for index, propeller in enumerate(previous_propellers):

                    # check how much the current propeller overlaps with the current repeat type and assign it the type with the largest overlap
                    overlap_to_repeat = set(range(propeller[0], propeller[1])).intersection(set(range(curr_repeat_types_intervals[repeat_type][0], curr_repeat_types_intervals[repeat_type][1])))
                    
                    if len(overlap_to_repeat) > propelles_types[str(propeller)]['repeat_type_overlap']:

                        propelles_types[str(propeller)]['repeat_type'] = repeat_type
                        propelles_types[str(propeller)]['repeat_type_overlap'] = len(overlap_to_repeat)
                        propelles_types[str(propeller)]['corresponding_repeats'] = curr_intervals
                        propelles_types[str(propeller)]['corresponding_blades'] = previous_blades_intervals[index]
                        propelles_types[str(propeller)]['name'] = previous_annotations[index]

            # select only those repeat types that overlap with propellers
            types_with_propellers = list(set([propelles_types[propeller]['repeat_type'] for propeller in propelles_types if propelles_types[propeller]['repeat_type'] != '']))
            
            if len(types_with_propellers) > 0:

                for repeat_type in types_with_propellers:

                    repeats_found = [interval for x, interval in enumerate(all_repeats_found) if repeats_types[x] == repeat_type]
                    previous_propellers_in_type = [ast.literal_eval(propeller) for propeller in propelles_types if propelles_types[propeller]['repeat_type'] == repeat_type]

                    # get the general blade lengths for each propeller previously identified
                    general_blade_length = np.median([(propeller[-1]-propeller[0])/int(propelles_types[str(propeller)]['name'].split()[-1].split('-')[0]) for propeller in previous_propellers_in_type])
                    
                    print([(propeller[-1]-propeller[0]) for propeller in previous_propellers_in_type])
                    print([(propeller[-1]-propeller[0])/int(propelles_types[str(propeller)]['name'].split()[-1].split('-')[0]) for propeller in previous_propellers_in_type])
                    print(general_blade_length)

                    # check if the detected repeats correspond to a full propeller or to individual blades 
                    # (in cases where there are multiple propellers, it is possible that it identifies repeats as the propeller and not the blades)
                    median_blades_per_repeat = np.median([round((repeat[-1]-repeat[0])/general_blade_length) for i, repeat in enumerate(repeats_found)])
                    
                    print([round((repeat[-1]-repeat[0])) for i, repeat in enumerate(repeats_found)])
                    print([round((repeat[-1]-repeat[0])/general_blade_length) for i, repeat in enumerate(repeats_found)])
                    print(median_blades_per_repeat)

                    if int(median_blades_per_repeat) <= 1.0: # each repeat could be a blade. 

                        if len(previous_propellers_in_type) == 1: 
                            # if only one previous propeller is inside this chunk
                            # check the number of repeats in the new and the old propeller. if it is bigger now, update it
                            propeller_previous_blades = propelles_types[str(previous_propellers_in_type[0])]['corresponding_blades']

                            if len(repeats_found) > len(propeller_previous_blades): 
                                blades_intervals = repeats_found
                            else: # else, keep the old one
                                blades_intervals = propeller_previous_blades

                        else: 
                            # this chunk covers more than one propeller previously identified.
                            blades_intervals = repeats_found

                    else: 
                    # we are dealing with a repeat unit corresponding to multiple blades. 
                    # we should split it, but this should never happen 
                        blades_intervals = repeats_found


                    # check if there is velcro closure
                    if len(blades_intervals) > 2: # if there is at least one middle blade
                        average_middle_blade_length = np.median([blade[-1]-blade[0] for blade in blades_intervals[1:-1]])

                        # if the sum of the lengths of the 2 last blades is as large as the median blade length, then there
                        # is velcro and the number of blades shall be reduced by 1
                        if round(((blades_intervals[-1][-1]-blades_intervals[-1][0])+(blades_intervals[0][-1]-blades_intervals[0][0]))/average_middle_blade_length) == 1:
                            new_annotation = 'beta-propeller {}-blades'.format(len(blades_intervals)-1)
                        else:    
                            new_annotation = 'beta-propeller {}-blades'.format(len(blades_intervals))
                    else:
                        new_annotation = 'beta-propeller {}-blades'.format(len(blades_intervals))

                    if blades_intervals not in new_blades_intervals:
                        new_domain_annotations.append(new_annotation)         
                        new_blades_intervals.append(blades_intervals)
        except:
            print('ERROR')
            print(tmp_label, updated_interval, protein_sequence)
            print(ncbi_code_data)


    new_propellers_intervals = [[propeller_block[0][0], propeller_block[-1][-1]] for propeller_block in new_blades_intervals]

    ncbi_code_data['updated_propeller_annotation'] = new_domain_annotations
    ncbi_code_data['updated_intervals'] = new_propellers_intervals
    ncbi_code_data['updated_blades_intervals'] = new_blades_intervals

    return ncbi_code_data

def update_propellers_blades_similarity(ncbi_code_data, tmp_label, min_seqID = args.min_seqID):

    #print('... ... Calculating blades similarities for {} propellers'.format(len(ncbi_code_data['updated_propeller_annotation'])))

    protein_sequence = ncbi_code_data['seq']
    ncbi_code_data['propellers_median_blade_similarity'] = []
    indeces_to_remove = []

    for i, propeller in enumerate(ncbi_code_data['updated_propeller_annotation']):
        propeller_blades_intervals = ncbi_code_data['updated_blades_intervals'][i]

        # save sequences to file
        tmp_fasta = '{}{}_propeller{}.fasta'.format(args.tmp_folder, tmp_label, i)
        with open(tmp_fasta, 'w') as fp:
            for j, blade in enumerate(propeller_blades_intervals):
                blade_seq = protein_sequence[blade[0]-1:blade[1]]
                fp.write('>blade_{}\n{}\n'.format(j, blade_seq))

        # align blades with muscle
        run_muscle = sp.Popen([muscle, '-in', tmp_fasta], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_muscle.communicate()
    
        # parse clustal output
        alignment = parse_muscle_out(stdout.decode('ascii'))

        # get blade identities matrix
        seq_IDmatrix = compute_identity_matrix(alignment)

        # get the median sequence identity
        curr_values = []
        for a, blade_a in enumerate(seq_IDmatrix):
            for b, blade_b in enumerate(seq_IDmatrix):
                if a > b:
                    curr_values.append(seq_IDmatrix[a][b])
        sequence_identitity = np.median(seq_IDmatrix)

        if sequence_identitity >= min_seqID:
            ncbi_code_data['propellers_median_blade_similarity'].append(sequence_identitity)
        else:
            indeces_to_remove.append(i)

        # remove temporary files
        os.remove(tmp_fasta)

    # update the dicitionary by removing those 'propellers' with a blade similarity lower than the requested (< 60%)
    ncbi_code_data['updated_propeller_annotation'] = [ncbi_code_data['updated_propeller_annotation'][i] for i in range(len(ncbi_code_data['updated_propeller_annotation'])) if i not in indeces_to_remove]
    ncbi_code_data['updated_intervals'] = [ncbi_code_data['updated_intervals'][i] for i in range(len(ncbi_code_data['updated_intervals'])) if i not in indeces_to_remove] 
    ncbi_code_data['updated_blades_intervals'] = [ncbi_code_data['updated_blades_intervals'][i] for i in range(len(ncbi_code_data['updated_blades_intervals'])) if i not in indeces_to_remove]

    return ncbi_code_data

def parse_muscle_out(clustal_out):
    
    alignment = []

    seq = ''
    for line in clustal_out.split('\n'):
        if line.startswith('>'):
            if seq != '':
                alignment.append(seq)
            seq = ''
        else:
            seq += line.strip()
    alignment.append(seq)
    
    return alignment

def compute_identity_matrix(alignment):
    
    identity_matrix = [[0 for seq in alignment] for seq in alignment]
    
    for i, seqi in enumerate(alignment):
        for j, seqj in enumerate(alignment):
            identity = 0
            alignment_length = 0
            for aa in range(len(seqi)):
                if seqi[aa] != '-' and seqj[aa] != '-':
                    alignment_length += 1
                    if seqi[aa] == seqj[aa]:
                        identity += 1

            if alignment_length > 0:
                identity_matrix[i][j] = identity*100/alignment_length
        
    return identity_matrix

## MAIN ROUTINES

def update_propellers_composition(arguments):

    jobID = arguments[0]
    ncbi_codes = arguments[1]
    all_sequences = arguments[2]
    
    print('\nANNOTATING SEQUENCES (thread {}):\n'.format(jobID+1))
    
    new_sequences = {}
    for i, ncbi_code in enumerate(ncbi_codes):      
        ncbi_code_data = all_sequences[ncbi_code]   

        print(' ... Thread {}: Working on {} ({}/{})'.format(jobID+1, ncbi_code, i+1, len(ncbi_codes)))
        
        thread_start = time.time()
        
        ncbi_code_data = annotate_repeats_in_sequence(ncbi_code_data, tmp_label = ncbi_code)
        ncbi_code_data = update_propellers_intervals(ncbi_code_data)
        ncbi_code_data = update_propellers_annotations(ncbi_code_data, tmp_label = ncbi_code)
        ncbi_code_data = update_propellers_blades_similarity(ncbi_code_data, tmp_label = ncbi_code)
                   
        new_sequences[ncbi_code] = ncbi_code_data

        thread_end = time.time()
        numb_seconds = thread_end - thread_start

        print(' ... Thread {}: Finished working on {} ({}/{}) after {} days {}'.format(jobID+1, ncbi_code, i+1, len(ncbi_codes), int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
            
    return new_sequences

def write_updated_files(all_sequences, taxonomy):

    print('\nSAVING EVERYTHING TO FILES\n')

    # create output folder
    out_folder = '{}/results_after_validation_with_hhrepID'.format(curr_directory)
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    # define output files
    out_full_length = '{}/all_searches_concatenated_full_length_sequences_after_validation_more_than_{}_seqID.fasta'.format(out_folder, args.min_seqID)
    out_propellers = '{}/all_searches_concatenated_propellers_sequences_after_validation_more_than_{}_seqID.fasta'.format(out_folder, args.min_seqID)
    
    out_tax = '{}/all_searches_concatenated_taxonomy_after_validation_more_than_{}_seqID.json'.format(out_folder, args.min_seqID)
    out_overall = '{}/all_searches_concatenated_validation_results_more_than_{}_seqID.json'.format(out_folder, args.min_seqID)

    taxonomy = json.load(open(taxonomy, 'r'))

    updated_taxonomy = {}

    with open(out_full_length, 'w') as outfull:
        with open(out_propellers, 'w') as outprop:
            for superkingdom in taxonomy:
                if superkingdom not in updated_taxonomy:
                    updated_taxonomy[superkingdom] = {}

                for phylum in taxonomy[superkingdom]:
                    if phylum not in updated_taxonomy[superkingdom]:
                        updated_taxonomy[superkingdom][phylum] = {}

                    for taxclass in taxonomy[superkingdom][phylum]:
                        if taxclass not in updated_taxonomy[superkingdom][phylum]:
                            updated_taxonomy[superkingdom][phylum][taxclass] = {}

                        for order in taxonomy[superkingdom][phylum][taxclass]:
                            if order not in updated_taxonomy[superkingdom][phylum][taxclass]:
                                updated_taxonomy[superkingdom][phylum][taxclass][order] = {}

                            for genus in taxonomy[superkingdom][phylum][taxclass][order]:
                                if genus not in updated_taxonomy[superkingdom][phylum][taxclass][order]:
                                    updated_taxonomy[superkingdom][phylum][taxclass][order][genus] = {}

                                for species in taxonomy[superkingdom][phylum][taxclass][order][genus]:
                                    if species not in updated_taxonomy[superkingdom][phylum][taxclass][order][genus]:
                                        updated_taxonomy[superkingdom][phylum][taxclass][order][genus][species] = {key: [] for key in taxonomy[superkingdom][phylum][taxclass][order][genus][species]}
                                        updated_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Propellers_labels'] = []

                                    for ncbi_code in taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes']:
                                        if ncbi_code in all_sequences:
                                            medians = all_sequences[ncbi_code]['propellers_median_blade_similarity']

                                            if len(medians) > 0:

                                                intervals = all_sequences[ncbi_code]['updated_blades_intervals']
                                                propellers_labels = all_sequences[ncbi_code]['updated_propeller_annotation']
                                                updated_intervals = all_sequences[ncbi_code]['updated_intervals']

                                                df = {'intervals': intervals, 'medians': medians, 'propellers_labels': propellers_labels, 'updated_intervals': updated_intervals}
                                                df = pd.DataFrame(df)
                                                df = df.sort_values(by='intervals')

                                                all_blades = [blade for propeller in df['intervals'] for blade in propeller]
                                                propellers = [[all_blades.index(blade) for blade in propeller] for propeller in df['intervals']]

                                                intervals = list(df.intervals)
                                                propellers_labels = list(df.propellers_labels)
                                                medians = list(df.medians)
                                                updated_intervals = list(df.updated_intervals)

                                                # update the taxonomy dictionary
                                                updated_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes'].append(ncbi_code)
                                                updated_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Medians'].append(medians)
                                                updated_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Propellers'].append(propellers)
                                                updated_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Intervals'].append(intervals)
                                                updated_taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Propellers_labels'].append(propellers_labels)

                                                # write full length sequence to fasta
                                                header = '>{}|{}|{}|{}|{}|{}|{}| Intervals: {} | Medians: {} | Propellers:{} | Propellers_labels:{}'.format(ncbi_code, superkingdom, phylum, taxclass, genus, order, species, intervals, medians, propellers, propellers_labels)
                                                sequence = all_sequences[ncbi_code]['seq']
                                                outfull.write('{}\n{}\n'.format(header, sequence))

                                                # write propellers sequences
                                                for i, propeller_interval in enumerate(updated_intervals):
                                                    propeller_sequence = sequence[propeller_interval[0]-1:propeller_interval[1]]
                                                    propeller_similarity = medians[i]
                                                    propeller_label = propellers_labels[i]
                                                    propeller_length = int(propeller_label.split()[-1].split('-')[0])

                                                    header = '>{}_#{}({} blades: {}%)|{}|{}|{}|{}|{}|{}'.format(ncbi_code, i+1, propeller_length, round(propeller_similarity,2), superkingdom, phylum, taxclass, genus, order, species)
                                                    outprop.write('{}\n{}\n'.format(header, propeller_sequence))

    # save updated json files
    json.dump(all_sequences, open(out_overall, 'w'), indent = 4)
    json.dump(updated_taxonomy, open(out_tax, 'w'), indent = 4)


## MAIN CODE

start = time.time()

out_folder = '{}/results_after_validation_with_hhrepID'.format(curr_directory)
if os.path.isfile('{}/all_searches_concatenated_validation_results_more_than_{}_seqID.json'.format(out_folder, args.min_seqID)) and args.force_search == 'False':
    
    print('Searches were already preformed, will only write files')

    all_sequences = json.load(open('{}/all_searches_concatenated_validation_results_more_than_{}_seqID.json'.format(out_folder, args.min_seqID)))

else:
    # 1. Get sequences
    all_sequences = parse_fasta(args.in_fasta)

    # 2. Find repeats in input sequences
    ncbi_codes = [ncbi_code for ncbi_code in all_sequences.keys()]
    separated_jobs = chunk_list(ncbi_codes, args.n_threads)

    list_arguments = [i for i in zip(range(args.n_threads), separated_jobs, [all_sequences for job in separated_jobs])]

    start = time.time()
    pool = mp.Pool(args.n_threads)
    results = pool.map(update_propellers_composition, list_arguments)
    all_sequences = {key: dic[key] for dic in results for key in dic.keys()}

    pool.close()
    pool.join()

# 3. Write sequences to fasta files and update taxonomy dictionary
write_updated_files(all_sequences, args.taxonomy)

end = time.time()
numb_seconds = end - start
print("END: Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
