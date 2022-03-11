import os
import ast
import random
import json
import argparse
import time

import numpy as np
import pandas as pd
import seaborn as sns
import subprocess as sp
import matplotlib.pyplot as plt
import multiprocessing as mp

from scipy import stats
from scipy.signal import find_peaks

# Define tools locations
muscle = '/Applications/muscle/muscle'
tmp_folder = '/tmp/'

# Get inputs
parser = argparse.ArgumentParser(description='')

### ... Required inputs
parser.add_argument('in_frags', metavar='in_frags', type=str, help='the json file with the hits matched (e.g., the blades)')
parser.add_argument('in_full', metavar='in_full', type=str, help='the fasta file with the full sequences')
# ... Optional inputs
parser.add_argument('-cpu', dest='cpu', default = 2, type=int, help='number of cpus to use (default: 2)')
parser.add_argument('-min_seqID', dest='min_seqID', default = 60, type=int, help='the minimum median sequence ID of the blades for a propeller to be considered (default: 60)')
parser.add_argument('-exclude_partial', dest='exclude_partial', default = True, type=bool, help='decide if you want to exclude sequences annotated with partial (default: True)')
parser.add_argument('-exclude_with_single_blade', dest='exclude_with_single_blade', default = True, type=bool, help='decide if you want to exclude sequences that have propellers that are a single blade (default: True)')
parser.add_argument('-label', dest='label', default = 'default',type=str, help='the label to append to the out folder')
parser.add_argument('-working_directory', dest='working_directory', default = '.',type=str, help='the working directory')

args = parser.parse_args()
if args.working_directory == '.':
    curr_directory = os.getcwd()
else:
    curr_directory = args.working_directory

# ROUTINES

# 0. Helping math

def median_absolute_deviation(x):
    return np.median(np.abs(x-np.median(x)))

# 1. For parsing blade and full length sequences (including ploting)

def convert_interval(interval_string):
    
    interval = interval_string.split()
    
    if len(interval) > 2:
        interval.remove(interval[0])

    if '[' in interval[0]:
        interval[0] = interval[0].replace('[','')
    if ']' in interval[1]:
        interval[1] = interval[1].replace(']','')

    interval[0] = int(interval[0])
    interval[1] = int(interval[1])
    
    return interval

def parse_blade_sequences_file(infasta):

    print('Parsing blades from {}'.format(infasta))
    
    sequences = {}
    ncbi_entrez = 'nan'

    with open(infasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>') and 'pdb' not in line and '|sp|' not in line:
                line_data = line[1:].split('|')
                
                ncbi_entrez = line_data[0].split('_#')[0]

                e_value = float(line_data[3].split(': ')[-1].strip())
                title = line_data[2]
                cur_interval = line_data[4].split(': ')[-1]

                interval = convert_interval(cur_interval)

                if ncbi_entrez not in sequences:
                    sequences[ncbi_entrez] = {'Intervals': [], 'Evalues': [], 'Sequences': [], 'Title': title}
                
                sequences[ncbi_entrez]['Intervals'].append(interval)
                sequences[ncbi_entrez]['Evalues'].append(e_value)
            
            elif ncbi_entrez != 'nan':
                sequence = line.strip()
                sequences[ncbi_entrez]['Sequences'].append(sequence)
                ncbi_entrez = 'nan'
    
    return sequences

def parse_full_sequences_file(infasta):

    print('Parsing full sequences from {}'.format(infasta))
    
    full_sequences = {}
    sequence_lengths = []
    ncbi_entrez = 'nan'
    
    with open(infasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>') and 'pdb' not in line and '|sp|' not in line:
                line_data = line[1:].split('|')
                
                ncbi_entrez = line_data[0].split('_#')[0]
                species = line_data[1]

                full_sequences[ncbi_entrez] = {}
                
            elif ncbi_entrez != 'nan':
                sequence = line.strip()
                full_sequences[ncbi_entrez]['Sequence'] = sequence
                full_sequences[ncbi_entrez]['Species'] = species
                sequence_lengths.append(len(sequence))
                ncbi_entrez = 'nan'
    
    return full_sequences

def get_blade_median_length(sequences):
    
    print('Getting median blade length')

    blade_lengths = []

    for ncbi_code in sequences:
        seqs = sequences[ncbi_code]['Sequences']
        
        for seq in seqs:
            blade_lengths.append(len(seq))
    
    median_length = np.median(blade_lengths)
    mad_length = round(median_absolute_deviation(blade_lengths))

    print(' ... The median blade length is {} +/- {} aa (error is median absolute deviation)'.format(median_length, mad_length))

    return median_length, mad_length

def get_linker_median_length(sequences, find_peaks = True):
    
    print('Getting median linker length')

    blade_distances = []

    for ncbi_code in sequences:
        intervals = sequences[ncbi_code]['Intervals']
        
        if len(intervals) > 1:
            for i, interval in enumerate(intervals[:-1]):
                curr_end = interval[1]
                next_start = intervals[i+1][0]

                dist = next_start - curr_end
                if dist > 0:
                    blade_distances.append(dist)

    if find_peaks:
        median_distance, mad_distance = find_major_peak(blade_distances)
        message = ' ... The most common linker length is {} +/- {} aa (error is the width of the peak)'.format(median_distance, mad_distance)
    else:
        median_distance = np.median(blade_distances)
        mad_distance = round(stats.median_absolute_deviation(blade_distances))
        message = ' ... The median linker length is {} +/- {} aa (error is median absolute deviation)'.format(median_distance, mad_distance)

    print(message)

    return median_distance, mad_distance

def find_major_peak(blade_distances):

    lengths = range(min(blade_distances), max(blade_distances)+1)
    counts = np.array([blade_distances.count(value) for value in lengths])

    height = int(len(blade_distances)*0.01)
    peaks, properties = find_peaks(counts, height=height, distance=1, width=1)

    highest_peak = 0
    highest_count = 0
    peak_count = 0
    for i, peak in enumerate(peaks):
        if counts[peak] > highest_count:
            highest_peak = peak
            highest_count = counts[peak]
            peak_count = i
        
    median_distance = lengths[highest_peak]
    mad_distance = round(properties['widths'][peak_count])

    return median_distance, mad_distance

# 2. To find propellers

def add_consecutive_blades(sequences, distance_theshold):

    print('Finding propellers using a distance threshold of up to {} to consider 2 blades as consecutive'.format(distance_theshold))
    
    for ncbi_code in sequences:
        intervals = sequences[ncbi_code]['Intervals']
        
        sequences[ncbi_code]['Consecutive_blades'] = []

        if len(intervals) > 1:
            for i, interval in enumerate(intervals[:-1]):
                
                if i == 0:
                    curr_propeller = [i]
                    
                curr_end = interval[1]
                next_start = intervals[i+1][0]

                dist = next_start - curr_end
                
                if dist <= distance_theshold:
                    curr_propeller.append(i+1)
                else:
                    sequences[ncbi_code]['Consecutive_blades'].append(curr_propeller)
                    curr_propeller = [i+1]     
            sequences[ncbi_code]['Consecutive_blades'].append(curr_propeller)
        
        else:
            sequences[ncbi_code]['Consecutive_blades'] = ['nan']
    
    return sequences

def exclude_propellers_with_single_blade(sequences):

    print('Excluding sequences that contain at least one propeller with only one blade')
    
    new_sequences = {}

    for ncbi_code in sequences:
        propellers = sequences[ncbi_code]['Consecutive_blades']

        has_single = False

        for propeller in propellers:
            if len(propeller) == 1 or propeller == 'nan':
                has_single = True

        if not has_single:
            new_sequences[ncbi_code] = sequences[ncbi_code]

    return new_sequences

def fix_blade_sequences(sequences, full_length_sequences):

    print('Fixing blade sequences so that they are complete')

    new_sequences = {}

    for ncbi_code in sequences:
        propellers = sequences[ncbi_code]['Consecutive_blades']
        intervals = sequences[ncbi_code]['Intervals']

        sequence = full_length_sequences[ncbi_code]['Sequence']

        new_sequences[ncbi_code] = sequences[ncbi_code]

        for propeller in propellers:
            for blade in propeller[1:-1]:
                curr_interval = [intervals[blade-1][-1]+1, intervals[blade][1]]
                curr_sequence = sequence[curr_interval[0]-1:curr_interval[1]]
                
                new_sequences[ncbi_code]['Intervals'][blade] = curr_interval
                new_sequences[ncbi_code]['Sequences'][blade] = curr_sequence            

    return new_sequences


# 3. To add sequence identity matrix

def chunk_list(l, n):
    chunks = np.array_split(np.array(l), n)
    chunks = [list(chunk) for chunk in chunks]
    return chunks

def align_blades(blade_sequences, jobid):

    random_number = random.randrange(100000)
        
    # write sequences to temporary fasta
    tmp_fasta = '{}tmp_blades_{}_{}_{}.fasta'.format(tmp_folder,jobid, random_number, args.label)
    with open(tmp_fasta, 'w') as tmpfasta:
        for i, seq in enumerate(blade_sequences):
            tmpfasta.write('>blade_{}\n{}\n'.format(i, seq))
            
    # run muscle for sequences
    run_muscle = sp.Popen([muscle, '-in', tmp_fasta], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_muscle.communicate()
    
    # parse clustal output
    alignment = parse_muscle_out(stdout.decode('ascii'))

    os.system('rm {}'.format(tmp_fasta))
    
    return alignment

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

def add_propellers_IDmatrix(arguments):

    jobid = arguments[0]
    ncbi_codes = arguments[1]
    sequences = arguments[2]

    new_sequences = {}

    print('Job {}. Calculating blades median seqID'.format(jobid))
    
    for count, ncbi_code in enumerate(ncbi_codes):

        new_sequences[ncbi_code] = sequences[ncbi_code]
        
        propellers = sequences[ncbi_code]['Consecutive_blades']
        blade_sequences = sequences[ncbi_code]['Sequences']
        
        new_sequences[ncbi_code]['BladesID_median'] = []
        
        if len(propellers) > 0:
            for i, propeller in enumerate(propellers):
                if len(propeller) > 1 and propeller != 'nan':
                                        
                    curr_blades = blade_sequences[propeller[0]:propeller[-1]+1]
                    curr_alignment = align_blades(curr_blades, jobid)
                    curr_IDmatrix = compute_identity_matrix(curr_alignment)

                    curr_values = []
                    for a, blade_a in enumerate(curr_IDmatrix):
                        for b, blade_b in enumerate(curr_IDmatrix):
                            if a > b:
                                curr_values.append(curr_IDmatrix[a][b])
                    sequence_identitity = np.median(curr_values)
                    
                    new_sequences[ncbi_code]['BladesID_median'].append(sequence_identitity)
                                   
                else:
                    print(' ... ... {} propeller has {} blades {}'.format(ncbi_code, len(propeller), propeller))
                    new_sequences[ncbi_code]['BladesID_median'].append('nan')
        else:
            print(' ... ... {} has {} propellers'.format(ncbi_code, len(propellers)))
            new_sequences[ncbi_code]['BladesID_median'].append('nan')

        print("... add_propellers_IDmatrix: Job {}. {} ({}/{}) ({}%)".format(jobid, ncbi_code, count, len(ncbi_codes), round(count*100/len(ncbi_codes))))
        
    return new_sequences

# 4. To select relevant sequences

def select_relevant_ncbi_codes(sequences, full_length_sequences, out_dir, min_seqID = args.min_seqID, exclude_partial = args.exclude_partial):

    print("Selecting relevant sequences with at least one propeller with a median blade seqID of {}%".format(min_seqID))
    
    new_sequences = {}

    out_fasta = "{}/selected_relevant_sequences_full_length.fasta".format(out_dir)

    count = 0
    with open(out_fasta, 'w') as out:
        for ncbi_code in sequences:
            medianIDs = sequences[ncbi_code]['BladesID_median']

            if any(i >= min_seqID for i in medianIDs):
                full_length_seq = full_length_sequences[ncbi_code]['Sequence']
                species = full_length_sequences[ncbi_code]['Species']
                
                title = sequences[ncbi_code]['Title']
                intervals = sequences[ncbi_code]['Intervals']
                propellers = sequences[ncbi_code]['Consecutive_blades']

                if not exclude_partial or 'partial' not in title.lower():
                    prot_identifier = '>{}|{}|{}| Medians: {}| Propellers: {}| Intervals: {}'.format(ncbi_code, species, title, [round(i, 2) for i in medianIDs], propellers, intervals)
                    out.write('{}\n'.format(prot_identifier))
                    out.write('{}\n'.format(full_length_seq))
                    count += 1

                    new_sequences[ncbi_code] = sequences[ncbi_code]

    if len(new_sequences) == 0:
        os.system('rm {}'.format(out_fasta))
        
    print(" ... Wrote {} sequences".format(count))
    
    return new_sequences

# 5. To make plots

def plot_relevant_propeller_histogram(sequences, out_dir, min_seqID = args.min_seqID):

    propellers_sizes = []

    for ncbi_code in sequences:
        medianIDs = sequences[ncbi_code]['BladesID_median']
        propellers = sequences[ncbi_code]['Consecutive_blades']

        for i, median in enumerate(medianIDs):
            if median >= min_seqID:
                propellers_sizes.append(len(propellers[i]))
    
    x = sorted(set(propellers_sizes))
    y = [propellers_sizes.count(i) for i in x]

    plt.figure(figsize=(int(len(x)/2),5))
    plt.bar(x,y)
    plt.yscale('log', nonposy='clip')
    plt.xticks(x)
    plt.title('Histogram of propellers with median blade ID >= {}%'.format(min_seqID))
    plt.savefig('{}/relevant_propellers_histogram.pdf'.format(out_dir), format = 'pdf')

    with open('{}/relevant_propellers_histogram.txt'.format(out_dir), 'w') as hist_txt:
        for i, value in enumerate(x):
            hist_txt.write('{}\t{}\n'.format(value, y[i]))
    
# 6. To write size representatives to doc file

def write_size_representatives(sequences, full_length_sequences, out_dir, min_seqID = args.min_seqID, topN = 1):

    propeller_sizes_topN = {}
    propeller_summary = {'ncbi_code': [], 'num_propellers': [], 'size': [], 'medianID': [], 'interval': [], 'full_seq_left': []}

    for ncbi_code in sequences:
        
        propellers = sequences[ncbi_code]['Consecutive_blades']
        seqIDs = sequences[ncbi_code]['BladesID_median']
        intervals = sequences[ncbi_code]['Intervals']
        
        if len(propellers) > 0:
            for i, medianID in enumerate(seqIDs):
                size = len(propellers[i])
                interval = [intervals[propellers[i][0]][0], intervals[propellers[i][-1]][-1]]
                    
                if medianID >= min_seqID:
                    propeller_summary['ncbi_code'].append(ncbi_code)
                    propeller_summary['size'].append(size)
                    propeller_summary['medianID'].append(medianID)
                    propeller_summary['interval'].append(interval)
                    propeller_summary['num_propellers'].append(len(propellers))
                    propeller_summary['full_seq_left'].append((len(full_length_sequences[ncbi_code]['Sequence'])-(interval[-1] - interval[0]))*100/len(full_length_sequences[ncbi_code]['Sequence']))
    
    propeller_summary = pd.DataFrame(propeller_summary)
    propeller_summary['score'] = propeller_summary['full_seq_left']*propeller_summary['medianID']
    propeller_summary = propeller_summary.sort_values(by = ['score'], ascending = False)

    tfile = open('{}/representatives_summary.txt'.format(out_dir), 'w')
    repfile = open('{}/representative_sequences.doc'.format(out_dir), 'w')
    
    for uniq_value in sorted(set(propeller_summary['size']), reverse = True):
        curr_propellers = propeller_summary.loc[propeller_summary['size'] == uniq_value].head(topN)
        tfile.write('{}\n'.format(curr_propellers))
        
        repfile.write('\n\n### {} BLADES\n\n'.format(uniq_value))
        
        for index, row in curr_propellers.iterrows():
            ncbi_code = row.ncbi_code
            interval = row.interval
            
            prot_identifier = '{}|{}|{}| Medians: {}| Propellers: {}| Intervals: {}'.format(ncbi_code, full_length_sequences[ncbi_code]['Species'], sequences[ncbi_code]['Title'], sequences[ncbi_code]['BladesID_median'], sequences[ncbi_code]['Consecutive_blades'], sequences[ncbi_code]['Intervals'])
            n_term = full_length_sequences[ncbi_code]['Sequence'][:interval[0]-1]
            propeller = full_length_sequences[ncbi_code]['Sequence'][interval[0]-1:interval[-1]-1]
            c_term = full_length_sequences[ncbi_code]['Sequence'][interval[-1]-1:]
            repfile.write('\n> {}\n\n{}\n\n{}\n\n{}\n\n'.format(prot_identifier, n_term, propeller, c_term))            

    tfile.close()
    repfile.close()
                   

# 7. To write N and C termini to fasta file

def write_NC_termini_fasta(sequences, full_length_sequences, out_dir, min_seqID = args.min_seqID, min_length = 10):

    outN = open('{}/Ntermini_from_selected_sequences.fasta'.format(out_dir), 'w')
    outC = open('{}/Ctermini_from_selected_sequences.fasta'.format(out_dir), 'w')

    for ncbi_code in sequences:
        
        propellers = sequences[ncbi_code]['Consecutive_blades']
        seqIDs = sequences[ncbi_code]['BladesID_median']
        intervals = sequences[ncbi_code]['Intervals']
        
        if len(propellers) > 0:
            for i, medianID in enumerate(seqIDs):
                size = len(propellers[i])
                interval = [intervals[propellers[i][0]][0], intervals[propellers[i][-1]][-1]]
                    
                if medianID >= min_seqID:

                    prot_identifier = '{}|{}|{}| Medians: {}| Propellers: {}| Intervals: {}'.format(ncbi_code, full_length_sequences[ncbi_code]['Species'], sequences[ncbi_code]['Title'], sequences[ncbi_code]['BladesID_median'], sequences[ncbi_code]['Consecutive_blades'], sequences[ncbi_code]['Intervals'])

                    n_term = full_length_sequences[ncbi_code]['Sequence'][:interval[0]-1]
                    c_term = full_length_sequences[ncbi_code]['Sequence'][interval[-1]-1:]

                    if len(n_term) >= min_length:
                        outN.write('\n>{}\n{}\n'.format(prot_identifier, n_term))
                    if len(c_term) >= min_length:
                        outC.write('\n>{}\n{}\n'.format(prot_identifier, c_term))

    outN.close()
    outC.close()
    
# START MAIN PIPELINE

# 1. Prepare for the run
out_dir = '{}/selection_of_relevant_sequences_{}'.format(curr_directory, args.label)
os.system('mkdir {}'.format(out_dir))

# 1.1. Parse the fragments json file, get their median length and the median linker length
fragments = json.load(open(args.in_frags, 'r'))
median_length, mad_length = get_blade_median_length(fragments)
median_linker, mad_linker = get_linker_median_length(fragments, find_peaks = False)

# 1.2. Parse full length sequences file
full_length_sequences = parse_full_sequences_file(args.in_full)

# 2. Find the propellers and add them to the dictionary
fragments = add_consecutive_blades(fragments, median_linker + 3*mad_linker)

# 2.1. Exclude those that have at least one propeller with just one blade (in case that option is turned on)
if args.exclude_with_single_blade:
    fragments = exclude_propellers_with_single_blade(fragments)

with open('{}/sequences_with_propellers.json'.format(out_dir), 'w') as fp:
    json.dump(fragments, fp, indent=4)

# 2.2. Fix the blade sequences so that they are complete
fragments = fix_blade_sequences(fragments, full_length_sequences)

with open('{}/sequences_with_fixed_sequence.json'.format(out_dir), 'w') as fp:
    json.dump(fragments, fp, indent=4)
    
# 3. Calculate the pairwise sequence identities between blades in the same propeller (do this in parallel)

ncbi_codes = [ncbi_code for ncbi_code in fragments.keys()]
separated_jobs = chunk_list(ncbi_codes, args.cpu)

list_arguments = [i for i in zip(range(args.cpu), separated_jobs, [fragments for job in separated_jobs])]

start = time.time()
pool = mp.Pool(args.cpu)
results = pool.map(add_propellers_IDmatrix, list_arguments)
fragments = {key: dic[key] for dic in results for key in dic.keys()}

pool.close()
pool.join()

end = time.time()
numb_seconds = end - start
print(" ... Finished after: {} ".format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

with open('{}/sequences_with_bladesID.json'.format(out_dir), 'w') as fp:
    json.dump(fragments, fp, indent=4)

# 4. Select only the ncbi_codes where at least one of the propellers has a median blade seqID larger than the threshold
#    and write their full length sequence to a fasta file
selected_sequences = select_relevant_ncbi_codes(fragments, full_length_sequences, out_dir, min_seqID = args.min_seqID, exclude_partial = args.exclude_partial)

if len(selected_sequences.keys()) > 0:

    print("\n ... ### Found {} sequences containing highly repetitive propellers".format(len(selected_sequences.keys())))
    
    with open('{}/selected_sequences_with_bladesID.json'.format(out_dir), 'w') as fp:
        json.dump(selected_sequences, fp, indent=4)

    ##selected_sequences = json.load(open('{}/selected_sequences_with_bladesID.json'.format(out_dir), 'r'))

    # 5. Make plots
    plot_relevant_propeller_histogram(selected_sequences, out_dir, min_seqID = args.min_seqID)

    # 6. Select size representatives and write to a doc file
    write_size_representatives(selected_sequences, full_length_sequences, out_dir, min_seqID = args.min_seqID, topN = 2)

    # 7. Write the N- and C-termini of the propellers into a fasta file
    write_NC_termini_fasta(selected_sequences, full_length_sequences, out_dir, min_seqID = args.min_seqID, min_length = median_length - 1*mad_length)

else:
    print("\n ... ### Did not find sequences containing highly repetitive propellers")
