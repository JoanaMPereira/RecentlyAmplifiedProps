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
hhalign = 'hhalign'
hhmake = 'hhmake'
muscle = '/Applications/muscle/muscle'
trimal = '/Applications/trimal-trimAl/source/trimal'

tmp_folder = '/tmp/'

# Get inputs
parser = argparse.ArgumentParser(description='')

# ... Required inputs
parser.add_argument('in_frags', metavar='in_frags', type=str, help='the fasta file with the hits matched (e.g., the blades)')
parser.add_argument('in_full', metavar='in_full', type=str, help='the fasta file with the full sequences')
# ... Optional inputs
parser.add_argument('-cpu', dest='cpu', default = 2,type=int, help='number of cpus to use (default: 2)')
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

def parse_full_sequences_file(infasta, out_dir):
    
    full_sequences = {}
    sequence_lengths = []
    ncbi_entrez = 'nan'

    with open(infasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>') and 'pdb' not in line and '|sp|' not in line:
                line_data = line[1:].split('|')
                
                ncbi_entrez = line_data[0].split('_#')[0]
                
            elif ncbi_entrez != 'nan':
                sequence = line.strip()
                full_sequences[ncbi_entrez] = sequence
                sequence_lengths.append(len(sequence))
                ncbi_entrez = 'nan'

    # plot their lengths distribution and get median and mad
    
    print(' ... Making length distribution plots')

    os.system('mkdir {}'.format(out_dir))

    plt.clf()
    plt.hist(sequence_lengths, bins = int((max(sequence_lengths)-min(sequence_lengths))/20))
    plt.savefig('{}/full_sequences_length_histogram.pdf'.format(out_dir), format = 'pdf')

    plt.clf()
    plt.hist(sequence_lengths, bins = int((max(sequence_lengths)-min(sequence_lengths))/20))
    plt.yscale('log', nonposy='clip')
    plt.savefig('{}/full_sequences_length_histogram_ylog10.pdf'.format(out_dir), format = 'pdf')

    plt.clf()
    plt.hist(sequence_lengths, bins = int((max(sequence_lengths)-min(sequence_lengths))/20))
    plt.yscale('log', nonposy='clip')
    plt.xlim(0, 2500)
    plt.savefig('{}/full_sequences_length_histogram_ylog10_crop.pdf'.format(out_dir), format = 'pdf')
    
    return full_sequences

def make_blade_length_hist(sequences, out_dir, label = 'initial'):
    
    print('Making blade length distribution plots')

    blade_lengths = []

    for ncbi_code in sequences:
        seqs = sequences[ncbi_code]['Sequences']
        
        for seq in seqs:
            blade_lengths.append(len(seq))

    median_length = np.median(blade_lengths)
    mad_length = round(median_absolute_deviation(blade_lengths))

    os.system('mkdir {}'.format(out_dir))

    plt.clf()
    plt.hist(blade_lengths, bins = max(blade_lengths)-min(blade_lengths))
    plt.yscale('log', nonposy='clip')
    plt.vlines(x = median_length, ymin = 0, ymax = 50000, color = 'black', linestyle = '--')
    plt.vlines(x = median_length + 3*mad_length, ymin = 0, ymax = 50000, color = 'red', linestyle = '--')
    plt.title('median: {} +/- {} aa (error is median absolute deviation)'.format(median_length, mad_length))
    plt.savefig('{}/{}_blade_length_histogram_ylog10.pdf'.format(out_dir, label), format = 'pdf')

    print(' ... The median blade length is {} +/- {} aa (error is median absolute deviation)'.format(median_length, mad_length))

    return median_length, mad_length

def make_linker_length_hist(sequences, out_dir, label = 'initial', find_peaks = True):
    
    print('Making linker length distribution plots')

    blade_distances = []

    for ncbi_code in sequences:
        intervals = sequences[ncbi_code]['Intervals']
        
        if len(intervals) > 1:
            for i, interval in enumerate(intervals[:-1]):
                curr_end = interval[1]
                next_start = intervals[i+1][0]

                dist = next_start - curr_end
                blade_distances.append(dist)

    if len(blade_distances) > 0:
        if find_peaks:
            median_distance, mad_distance = find_major_peak(blade_distances)
                        
            title = 'highest peak: {} +/- {} aa (error is the width of the peak)'.format(median_distance, mad_distance)
            message = ' ... The most common linker length is {} +/- {} aa (error is the width of the peak)'.format(median_distance, mad_distance)
        else:
            median_distance = np.median(blade_distances)
            mad_distance = round(stats.median_absolute_deviation(blade_distances))
                
            title = 'median: {} +/- {} aa (error is median absolute deviation)'.format(median_distance, mad_distance)
            message = ' ... The median linker length is {} +/- {} aa (error is median absolute deviation)'.format(median_distance, mad_distance)

        os.system('mkdir {}'.format(out_dir))

        plt.clf()
        plt.hist(blade_distances, bins = max(blade_distances)-min(blade_distances))
        plt.yscale('log', nonposy='clip')
        plt.vlines(x = median_distance, ymin = 0, ymax = 50000, color = 'black', linestyle = '--', linewidth = 0.75)
        plt.xlim(-10, 400)
        plt.title(title)
        plt.savefig('{}/{}_linker_length_histogram_crop.pdf'.format(out_dir, label), format = 'pdf')

        print(message)
    else:
        median_distance = 'nan'
        mad_distance = 'nan'
        print(" ... There are not enough linkers")

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
    
# 2. For making hhm profile

def make_blades_hhm(sequences, out_dir, label = '', blade_length = 38, confidence = 0, nmax = 200):
    
    print('Making HHM for {} random sequences'.format(nmax))
    
    # write n random blades that follow the median_length to a fasta file
    print(' ... Collecting all sequences')
    
    random_number_1 = random.randrange(1,100000)
    out_blades = '{}median_blades_{}_{}.fasta'.format(tmp_folder, args.label, random_number_1)
    sequences_selected = []
    
    for ncbi_code in sequences.keys():
        for blade in set(sequences[ncbi_code]['Sequences']):
            if len(blade) >= blade_length - confidence and len(blade) <= blade_length + confidence:
                sequences_selected.append(blade)

    if len(sequences_selected) == 0:
        for ncbi_code in sequences.keys():
            for blade in set(sequences[ncbi_code]['Sequences']):
                sequences_selected.append(blade)
   
                
    with open(out_blades, 'w') as outfasta:
        if nmax < len(sequences_selected): 
            random_sequences = random.sample(sequences_selected, nmax)
        else:
            random_sequences = sequences_selected
            
        for i, seq in enumerate(random_sequences):
            outfasta.write('>blade{}\n{}\n'.format(i, seq))
    
    # align blades with muscle
    print(' ... Aligning selected')

    run_muscle = sp.Popen([muscle, '-in', out_blades], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_muscle.communicate()
    
    # parse muscle output and write it to a file
    alignment = parse_muscle_out(stdout.decode('ascii'))
    out_alignment = '{}/median_blades_alignment{}.a2m'.format(out_dir, label)
    
    with open(out_alignment, 'w') as outal:
        for i, seq in enumerate(alignment):
            outal.write('>blade{}\n{}\n'.format(i, seq))

    # trim alignment with trimal
    print(' ... Triming alignment')
    
    out_trimal = "{}_trimal.a2m".format(out_alignment.split('.')[0])
    run_trimal1 = sp.Popen([trimal, '-in', out_alignment, '-out', out_trimal, '-gt', str(0.30)], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_trimal1.communicate()
    run_trimal2 = sp.Popen([trimal, '-in', out_trimal, '-out', out_trimal, '-resoverlap', str(0.75), '-seqoverlap', str(0.80)], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_trimal2.communicate()
    
    # make hhm with hhmake for this fasta file
    print(' ... Making HHM')

    out_hhm = '{}.hhm'.format(out_trimal.split('.')[0])
    make_hhm = sp.Popen([hhmake, '-i', out_trimal, '-cons', '-o', out_hhm], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = make_hhm.communicate()
    
    os.system('rm {}'.format(out_blades))

    return out_hhm
    

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

# 3. For unusual and missing blade searching

def chunk_list(l, n):
    chunks = np.array_split(np.array(l), n)
    chunks = [list(chunk) for chunk in chunks]
    return chunks

def align_blade_to_hhm(blade_seq, blades_hhm, jobid, label = 'tmp'):

    random_number_1 = random.randrange(1,100000)
    random_number_2 = random.randrange(1,100000)
    random_number_3 = random.randrange(1,100000)
    
    # align linker to hhm
    out_fasta = '{}{}_blade_{}_{}_{}.fasta'.format(tmp_folder, label, jobid, int(random_number_1*random_number_3/random_number_2), args.label)
    
    with open(out_fasta,'w') as tmpblade:
        tmpblade.write('>blade\n{}\n'.format(blade_seq))

    run_hhalign = sp.Popen([hhalign, '-i', out_fasta, '-t', blades_hhm, '-cov', '80'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_hhalign.communicate()
            
    if not os.path.isfile('{}.hhr'.format(out_fasta.split('.')[0])):
        print(out_fasta, stderr)

    results = parse_hhr('{}.hhr'.format(out_fasta.split('.')[0]))

    os.system('rm {}.*'.format(out_fasta.split('.')[0]))
        
    return results
                        
def parse_hhr(hhr_file):
        
    found_table = False
    results = {'Prob': [],
               'Evalue': [],
               'Query interval': [],
               'Template coverage': []}
    
    with open(hhr_file, 'r') as outhhalign:
        for line in outhhalign:
            if 'No Hit                             Prob E-value' in line:
                found_table = True
                
            elif found_table and len(line) < 4:
                found_table = False
            
            elif found_table and len(line) > 4:
                line_results = line.split()
                results['Prob'].append(float(line_results[2]))
                results['Evalue'].append(float(line_results[3]))
                results['Query interval'].append([int(line_results[8].split('-')[0]), int(line_results[8].split('-')[1])])
                
                tmpl_cov = (int(line_results[9].split('-')[1]) - int(line_results[9].split('-')[0]))/int(line_results[10].replace('(','').replace(')',''))
                
                results['Template coverage'].append(round(tmpl_cov*100,1))
                
    return results

def order_selected_intervals(selected_intervals):
    
    intervals = []
    
    for selected_interval in selected_intervals:
        intervals.append(selected_interval[0])
    
    ordered_selected_intervals = []
    
    for interval in np.sort(np.array(intervals),axis=0):
        for selected_interval in selected_intervals:
            if selected_interval[0] == list(interval):
                ordered_selected_intervals.append(selected_interval)
    
    return ordered_selected_intervals
        
def search_profile_iteratively(sequence, hhm, blade_length, confidence, jobid, min_prob, min_cov, selected_intervals = [], correction_factor = 0, label = 'tmp'):
    
    hhalign_results = align_blade_to_hhm(sequence, hhm, jobid, label = label)
    collected = False
        
    for j,prob in enumerate(hhalign_results['Prob']):
        if prob >= min_prob or hhalign_results['Template coverage'][j] >= min_cov:
            
            curr_interval = hhalign_results['Query interval'][j]
            found_interval = [curr_interval[0]+correction_factor, curr_interval[1]+correction_factor]
            evalue = hhalign_results['Evalue'][j]
            
            new_blade_seq = sequence[curr_interval[0]-1: curr_interval[1]]
            selected_intervals.append([found_interval, new_blade_seq, evalue])
            
            collected = True

    if collected:
        if len(new_blade_seq)/len(sequence) < 0.5 and (len(sequence) - len(new_blade_seq) > (blade_length - confidence)):
            
            if curr_interval[0] > blade_length - confidence:
                remaining_sequence = sequence[:curr_interval[0]]
                selected_intervals = search_profile_iteratively(remaining_sequence, hhm, blade_length, confidence, jobid, min_prob, min_cov, selected_intervals = selected_intervals, correction_factor = correction_factor)
           
            if (len(sequence) - curr_interval[1]) > blade_length - confidence:
                remaining_sequence = sequence[curr_interval[1]:]
                correction_factor = found_interval[1]
                selected_intervals = search_profile_iteratively(remaining_sequence, hhm, blade_length, confidence, jobid, min_prob, min_cov, selected_intervals = selected_intervals, correction_factor = correction_factor)
                        
    return selected_intervals


def fix_long_blades_with_HHalign(arguments):

    job_id = arguments[0]
    ncbi_codes = arguments[1]
    sequences = arguments[2]
    blades_hhm = arguments[3]
    blade_length = arguments[4]
    confidence = arguments[5]
    min_prob = arguments[6]
    min_cov = arguments[7]
    
    new_sequences = {}
    
    print('Job {}. Analysing blades'.format(job_id))
    
    for count, ncbi_code in enumerate(ncbi_codes):

        print(' ... fix_long_blades_with_HHalign: Job {}.{} ({}%)'.format(job_id, count, round(count*100/len(ncbi_codes))))
        
        curr_data = sequences[ncbi_code]
        new_sequences[ncbi_code] = {key:[] for key in curr_data.keys()}
        new_sequences[ncbi_code]['Title'] = curr_data['Title']
                
        for i, blade_seq in enumerate(curr_data['Sequences']):

            curr_interval = curr_data['Intervals'][i]
            curr_evals = curr_data['Evalues'][i]
            
            if len(blade_seq) >= blade_length+confidence*2: # analyse only the blades that are larger thant the reference length by 2 confidence intervals
                
                selected_intervals = search_profile_iteratively(blade_seq, blades_hhm, blade_length, confidence, jobid = job_id, min_prob = min_prob,min_cov = min_cov, selected_intervals = [], correction_factor = curr_interval[0]-1, label = ncbi_code.split('.')[0])
                    
                #add the new intervals to the dictionary       
                if len(selected_intervals) > 0:
                    ordered_intervals = order_selected_intervals(selected_intervals)
                    
                    for interval in ordered_intervals:
                        new_sequences[ncbi_code]['Sequences'].append(interval[1])
                        new_sequences[ncbi_code]['Evalues'].append(interval[2])
                        new_sequences[ncbi_code]['Intervals'].append(interval[0])
                  
                
                # if we continue having only 1 match, keep the default from the blast search
                else:
                    new_sequences[ncbi_code]['Sequences'].append(blade_seq)
                    new_sequences[ncbi_code]['Evalues'].append(curr_evals)
                    new_sequences[ncbi_code]['Intervals'].append(curr_interval)
                 
            else:
                new_sequences[ncbi_code]['Sequences'].append(blade_seq)
                new_sequences[ncbi_code]['Evalues'].append(curr_evals)
                new_sequences[ncbi_code]['Intervals'].append(curr_interval)
                
                    
    return new_sequences

def add_internal_missing_blades_with_HHalign(arguments):

    job_id = arguments[0]
    ncbi_codes = arguments[1]
    sequences = arguments[2]
    blades_hhm = arguments[3]
    blade_length = arguments[4]
    confidence = arguments[5]
    min_prob = arguments[6]
    min_cov = arguments[7]
    full_length_sequences = arguments[8]
    
    new_sequences = {}
    
    print('Job {}. Analysing linkers'.format(job_id))
    
    for count, ncbi_code in enumerate(ncbi_codes):

        print(' ... add_internal_missing_blades_with_HHalign: Job {}.{} ({}%)'.format(job_id, count, round(count*100/len(ncbi_codes))))
        
        curr_data = sequences[ncbi_code]
        new_sequences[ncbi_code] = {key:[] for key in curr_data.keys()}
        new_sequences[ncbi_code]['Title'] = curr_data['Title']
        
        intervals = curr_data['Intervals']
    
        if len(intervals) > 1:
            for i, interval in enumerate(intervals[:-1]):
                
                new_sequences[ncbi_code]['Sequences'].append(curr_data['Sequences'][i])
                new_sequences[ncbi_code]['Evalues'].append(curr_data['Evalues'][i])
                new_sequences[ncbi_code]['Intervals'].append(curr_data['Intervals'][i])
                        
                curr_end = interval[1]
                next_start = intervals[i+1][0]

                dist = next_start - curr_end
                
                if dist > median_length - confidence*2: # check the linkers that are at least as large as the median blade length minus 2 times de confidence interval
                    
                    linker_seq = full_length_sequences[ncbi_code][curr_end:next_start]
                    
                    selected_intervals = search_profile_iteratively(linker_seq, blades_hhm, blade_length, confidence, jobid = job_id, min_prob = min_prob,min_cov = min_cov, selected_intervals = [], correction_factor = curr_end, label = ncbi_code.split('.')[0])
                    
                    if len(selected_intervals) > 0:
                        ordered_intervals = order_selected_intervals(selected_intervals)
                    
                        for interval in ordered_intervals:
                            new_sequences[ncbi_code]['Sequences'].append(interval[1])
                            new_sequences[ncbi_code]['Evalues'].append(interval[2])
                            new_sequences[ncbi_code]['Intervals'].append(interval[0])
        
        new_sequences[ncbi_code]['Sequences'].append(curr_data['Sequences'][-1])
        new_sequences[ncbi_code]['Evalues'].append(curr_data['Evalues'][-1])
        new_sequences[ncbi_code]['Intervals'].append(curr_data['Intervals'][-1])
    
    return new_sequences

def add_NC_blades_with_HHalign(arguments):

    job_id = arguments[0]
    ncbi_codes = arguments[1]
    sequences = arguments[2]
    blades_hhm = arguments[3]
    blade_length = arguments[4]
    confidence = arguments[5]
    min_prob = arguments[6]
    min_cov = arguments[7]
    full_length_sequences = arguments[8]
    
    new_sequences = {}
    
    print('Job {}. Analysing N and C termini'.format(job_id))
    
    for count, ncbi_code in enumerate(ncbi_codes):

        print(' ... add_NC_blades_with_HHalign: Job {}.{} ({}%)'.format(job_id, count, round(count*100/len(ncbi_codes))))
        
        curr_data = sequences[ncbi_code]
        new_sequences[ncbi_code] = {key:[] for key in curr_data.keys()}
        new_sequences[ncbi_code]['Title'] = curr_data['Title']
        
        intervals = curr_data['Intervals']
    
        # Check N-terminus
        n_terminal_seq = full_length_sequences[ncbi_code][:intervals[0][0]]
        
        if len(n_terminal_seq) > median_length - confidence*2:
            selected_intervals = search_profile_iteratively(n_terminal_seq, blades_hhm, blade_length, confidence, jobid = job_id, min_prob = min_prob,min_cov = min_cov, selected_intervals = [], correction_factor = 0, label = ncbi_code.split('.')[0])

            if len(selected_intervals) > 0:
                ordered_intervals = order_selected_intervals(selected_intervals)

                for interval in ordered_intervals:
                    new_sequences[ncbi_code]['Sequences'].append(interval[1])
                    new_sequences[ncbi_code]['Evalues'].append(interval[2])
                    new_sequences[ncbi_code]['Intervals'].append(interval[0])
        
        # Add all blades found before
        for i, interval in enumerate(intervals):

            new_sequences[ncbi_code]['Sequences'].append(curr_data['Sequences'][i])
            new_sequences[ncbi_code]['Evalues'].append(curr_data['Evalues'][i])
            new_sequences[ncbi_code]['Intervals'].append(curr_data['Intervals'][i])
   
        # Check C-terminus
        c_terminal_seq = full_length_sequences[ncbi_code][intervals[-1][-1]:]
        
        if len(c_terminal_seq) > median_length - confidence*2:
            selected_intervals = search_profile_iteratively(c_terminal_seq, blades_hhm, blade_length, confidence, jobid = job_id, min_prob = min_prob,min_cov = min_cov, selected_intervals = [], correction_factor = intervals[-1][-1], label = ncbi_code.split('.')[0])

            if len(selected_intervals) > 0:
                ordered_intervals = order_selected_intervals(selected_intervals)

                for interval in ordered_intervals:
                    new_sequences[ncbi_code]['Sequences'].append(interval[1])
                    new_sequences[ncbi_code]['Evalues'].append(interval[2])
                    new_sequences[ncbi_code]['Intervals'].append(interval[0])
        
    return new_sequences

# START MAIN PIPELINE
min_cov = 50
min_prob = 70
nmax = 1000

# 1. Prepare for the searches
out_dir = '{}/missing_blade_searching_{}'.format(curr_directory, args.label)

# 1.1. Parse the fragments file, plot their length distribution and get their median size and associated median absolute deviation
fragments = parse_blade_sequences_file(args.in_frags)
median_length, mad_length = make_blade_length_hist(fragments, out_dir = out_dir, label = 'initial')

# 1.2. Parse full length sequences file
full_length_sequences = parse_full_sequences_file(args.in_full, out_dir = out_dir)

# 1.3. Make hhm profile for the fragments
fragments_hmm = make_blades_hhm(fragments, out_dir = out_dir, label = '_for_long_blades', blade_length = median_length, confidence = 1*mad_length, nmax = nmax)

# 2. Do the searches (and define the running parameters) by constantly updating the hmm
confidence = 3*mad_length

# 2.1. First search for unusual blades that may correspond to more than one blade, separating the jobs by the different cpus
ncbi_codes = [ncbi_code for ncbi_code in fragments.keys()]
separated_jobs = chunk_list(ncbi_codes, args.cpu)

list_arguments = [i for i in zip(range(args.cpu), separated_jobs, [fragments for job in separated_jobs], [fragments_hmm for job in separated_jobs], [median_length for job in separated_jobs], [confidence for job in separated_jobs], [min_prob for job in separated_jobs], [min_cov for job in separated_jobs])]

start = time.time()
pool = mp.Pool(args.cpu)
results = pool.map(fix_long_blades_with_HHalign, list_arguments)
fragments = {key: dic[key] for dic in results for key in dic.keys()}

pool.close()
pool.join()

end = time.time()
numb_seconds = end - start
print(" ... Finished after: {} ".format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

with open('{}/sequences_with_fixed_long_blades.json'.format(out_dir), 'w') as fp:
    json.dump(fragments, fp, indent=4)

# 2.2. Now analyse the linkers
fragments_hmm = make_blades_hhm(fragments, out_dir = out_dir, label = '_for_linkers', blade_length = median_length, confidence = 1*mad_length, nmax = nmax)

median_linker, mad_linker = make_linker_length_hist(fragments, out_dir = out_dir, label = 'initial', find_peaks = False)

if median_linker != 'nan':
    list_arguments = [i for i in zip(range(args.cpu), separated_jobs, [fragments for job in separated_jobs], [fragments_hmm for job in separated_jobs], [median_length for job in separated_jobs], [confidence for job in separated_jobs], [min_prob for job in separated_jobs], [min_cov for job in separated_jobs], [full_length_sequences for job in separated_jobs])]

    start = time.time()
    pool = mp.Pool(args.cpu)
    results = pool.map(add_internal_missing_blades_with_HHalign, list_arguments)
    fragments = {key: dic[key] for dic in results for key in dic.keys()}

    pool.close()
    pool.join()

    end = time.time()
    numb_seconds = end - start
    print(" ... Finished after: {} ".format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

    with open('{}/sequences_with_fixed_long_linkers.json'.format(out_dir), 'w') as fp:
        json.dump(fragments, fp, indent=4)

# 2.3. Finally, analyse the N and C termini for missed blades
fragments_hmm = make_blades_hhm(fragments, out_dir = out_dir, label = '_for_termini', blade_length = median_length, confidence = 1*mad_length, nmax = nmax)

list_arguments = [i for i in zip(range(args.cpu), separated_jobs, [fragments for job in separated_jobs], [fragments_hmm for job in separated_jobs], [median_length for job in separated_jobs], [confidence for job in separated_jobs], [min_prob for job in separated_jobs], [min_cov for job in separated_jobs], [full_length_sequences for job in separated_jobs])]

start = time.time()
pool = mp.Pool(args.cpu)
results = pool.map(add_NC_blades_with_HHalign, list_arguments)
fragments = {key: dic[key] for dic in results for key in dic.keys()}

pool.close()
pool.join()

end = time.time()
numb_seconds = end - start
print(" ... Finished after: {} ".format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

#2.4. Plot the histograms again and save json results
median_length, mad_length = make_blade_length_hist(fragments, out_dir = out_dir, label = 'final')
median_linker, mad_linker = make_linker_length_hist(fragments, out_dir = out_dir, label = 'final', find_peaks = False)

with open('{}/sequences_with_missed_terminal_blades.json'.format(out_dir), 'w') as fp:
    json.dump(fragments, fp, indent=4)

print('Saved results to {}/sequences_with_missed_terminal_blades.json'.format(out_dir))
