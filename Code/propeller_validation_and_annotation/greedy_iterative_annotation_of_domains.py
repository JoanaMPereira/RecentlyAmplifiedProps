''' Follows an iterative greedy algorithm to find domains in a set of input sequences'''

import multiprocessing as mp
import subprocess as sp
import numpy as np

import argparse
import time
import os
import ast
import json
import sys

hhsearch = 'hhsearch' 
hhdatabases_location = # add here

# Get inputs
parser = argparse.ArgumentParser(description='')
requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments with defaults')

# required inputs
requiredNamed.add_argument('-i', dest='in_fasta', type=str, required=True, help='input fasta to annotate')
# optional inputs
optionalNamed.add_argument('-db', dest='db_label', default='ECOD', type=str, help='HH-suite database to use (default: ECOD)')
optionalNamed.add_argument('-rounds', dest='rounds', default=1, type=int, help='Number of rounds to search for domains (default: 1)')
optionalNamed.add_argument('-min_len', dest='min_seq_length', default=40, type=int, help='Minimum length of a sequence to be annotated (default: 40 aminoacids)')
optionalNamed.add_argument('-min_prob', dest='min_prob', default=50, type=float, help='Minimum probability of matches to be considerd (default: 50%)')
optionalNamed.add_argument('-cpu', dest='n_cpu', default = 2,type=int, help='Number of cpus to use (default: 2)')
optionalNamed.add_argument('-threads', dest='n_threads', default=1, type=int, help='Number of parallel jobs (default: 1)')
optionalNamed.add_argument('-tmp_folder', dest='tmp_folder', default = '/tmp/',type=str, help='Temporary folder (default: /tmp/)')
optionalNamed.add_argument('-out_dir', dest='out_dir', default = '.',type=str, help='Output folder (default: .)')
optionalNamed.add_argument('-print_tname', dest='print_tname', default = 'False',type=str, help='Boolean to print t-name also in domain description (default: False)')

args = parser.parse_args()
curr_directory = os.getcwd()

if args.out_dir == '.':
    out_dir = curr_directory
else:
    out_dir = args.out_dir

print('INPUTS:')
for arg in vars(args):
    print(arg, getattr(args, arg))

# Define automatic logger

class Logger(object):
    def __init__(self, name):
        self.terminal = sys.stdout
        self.log = open(name, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    

# PRINT INPUTS SUMMARY
hhsearch_cpu = int(args.n_cpu/args.n_threads)

print('\nINPUT PARAMETERS:')
print(' ... Database:                                 {}'.format(args.db_label))
print(' ... Number of threads:                        {}'.format(args.n_threads))
print(' ... CPUs for HHsearch:                        {}'.format(hhsearch_cpu))
print(' ... Maximum no. of rounds per input sequence: {}'.format(args.rounds))
print(' ... Minimum probability of a match:           {}'.format(args.min_prob))
print(' ... Minimum sequence length to annotate:      {}'.format(args.min_seq_length))

## HELPING ROUTINES

# 1. Job preparation

def find_database(database_label):

    db_folder = '{}{}/'.format(hhdatabases_location, database_label)
    db_file = ['{}{}'.format(db_folder, fl) for fl in os.listdir(db_folder) if fl.endswith('_hhm.ffdata')][0]
    db = db_file.replace('_hhm.ffdata', '')
    db_version = time.ctime(os.path.getmtime(db_file))
    
    print('\nDATABASE: Searches will be performed over database {} (version from {})\n'.format(db.split('/')[-1], db_version))
    
    return db

def parse_fasta(in_fasta):

    sequences = {}
    annotated_propellers = 0

    with open(in_fasta, 'r') as infst:
        for line in infst:
            if len(line) > 0:
                line = line.strip()
                if line.startswith('>') or '|' in line:
                    ncbi_code = line.split('|')[0].replace('>','').split('_#A')[0]
                    sequences[ncbi_code] = {'seq': '', 'annotated_domains': [], 'annotated_intervals': [], 'annotated_probability':[], 'annotated_probability_diff':[]}

                    if 'Propellers:' in line:
                        intervals = ast.literal_eval(line.split('Intervals:')[-1].split('|')[0].strip())
                        intervals = [interval for intervals_block in intervals for interval in intervals_block]
                        propellers = ast.literal_eval(line.split('Propellers:')[-1].split('|')[0].strip())
                        propellers_labels = ast.literal_eval(line.split('Propellers_labels:')[-1].split('|')[0].strip())
                        
                        for i, propeller in enumerate(propellers):

                            propeller_interval = [intervals[propeller[0]][0], intervals[propeller[-1]][-1]]
                            propeller_label = propellers_labels[i]

                            sequences[ncbi_code]['annotated_domains'].append(propeller_label)
                            sequences[ncbi_code]['annotated_intervals'].append(propeller_interval)
                            annotated_propellers += 1
                else:
                    sequences[ncbi_code]['seq'] = line.strip()

    print("INPUT SEQUENCES: Found input {} sequences, totalling {} already annotated propellers".format(len(sequences), annotated_propellers))
    
    return sequences

# 2. Domain annotation

def chunk_list(l, n):
    chunks = np.array_split(np.array(l), n)
    chunks = [list(chunk) for chunk in chunks]
    return chunks

def find_sequence_stretches_to_annotate(ncbi_code_data):

    stretches_to_annotate = {}

    if len(ncbi_code_data['annotated_domains']) == 0:
        stretches_to_annotate = {'full': {'seq': ncbi_code_data['seq'], 'correction_factor': 0}}
    else:
        intervals = ncbi_code_data['annotated_intervals']
        for i, interval in enumerate(intervals):
            if i == 0:
                stretches_to_annotate['N'] = {'seq': ncbi_code_data['seq'][:interval[0]-1], 'correction_factor': 0}

            if i < len(ncbi_code_data['annotated_intervals'])-1:
                next_interval = intervals[i+1]
                stretches_to_annotate['gap{}'.format(i+1)] = {'seq': ncbi_code_data['seq'][interval[-1]-1:next_interval[0]-1], 'correction_factor': interval[-1]}

            else:
                stretches_to_annotate['C'] = {'seq': ncbi_code_data['seq'][interval[-1]-1:], 'correction_factor': interval[-1]}
    
    return stretches_to_annotate

def parse_hhsearch_out(out_hhr, min_hit_length):

    domains_found = {}

    found_table = False
    found_end_of_table = False

    curr_index = 0
    with open(out_hhr, 'r') as hhr:
        for line in hhr:
            if 'No Hit' in line:
                found_table = True
            elif 'No 1' in line:
                found_end_of_table = True

            if found_table and not found_end_of_table:
                if len(line.strip()) > 0 and 'No Hit' not in line:
                    hit_index = int(line.split()[0])
                    hit_stats = line[35:].split()
                    
                    hit_prob = float(hit_stats[0])
                    hit_interval = [int(i) for i in hit_stats[6].split('-')]
            
                    template_interval = hit_stats[7].split('(')[0]
                    template_interval = [int(i) for i in template_interval.split('-')]
                    try:
                        template_total_len = int(hit_stats[8].replace('(','').replace(')',''))
                    except:
                        template_total_len = int(hit_stats[7].split('(')[-1].replace('(','').replace(')',''))
                    template_coverage = (template_interval[-1]-template_interval[0])*100/template_total_len

                    if hit_interval[-1] - hit_interval[0] >= min_hit_length:
                        domains_found[hit_index] = {'Prob': hit_prob, 'Interval': hit_interval, 'Template_coverage': round(template_coverage, 2)}

            if line.startswith('>'):
                curr_index += 1
                correct_label = line.strip()[1:]
                if args.db_label == 'ECOD':
                    x_level = correct_label.split(', X: ')[-1].split(', ')[0].strip().replace('"','')
                    t_level = correct_label.split(', T: ')[-1].split(', ')[0].strip().replace('"','')
                    f_level = correct_label.split(', F: ')[-1].split('|')[0].strip().replace('"','')

                    if args.print_tname == 'True':
                        corr_label = 'X: {} T: {} F: {}'.format(x_level, t_level, f_level)
                    else:
                        corr_label = 'X: {} F: {}'.format(x_level, f_level)
                    
                    if x_level == 'NO_X_NAME' or f_level == 'F_UNCLASSIFIED':
                        a_level = correct_label.split('A: ')[-1].split(', ')[0].strip().replace('"','')
                        if x_level == 'NO_X_NAME' and f_level == 'F_UNCLASSIFIED':
                            corr_label = '{}'.format(a_level)
                        else:
                            corr_label = '{} {}'.format(a_level, corr_label)

                    correct_label = corr_label
                        
                elif args.db_label == 'pfama':
                    correct_label = correct_label.split(';')[-1].strip()
                    
                if curr_index in domains_found:
                    if template_coverage < 90:
                        correct_label = 'Partial ({}%) {}'.format(round(template_coverage), correct_label)
                        
                    domains_found[curr_index]['Domain_label'] = correct_label

    return domains_found

def find_significant_domains(domains_found):

    domains_to_remove = []
    probs_diff = {}

    for domain_i in domains_found:
        prob_i = domains_found[domain_i]['Prob']
        interval_i = domains_found[domain_i]['Interval']
        curr_prob_diffs = np.nan

        for domain_j in domains_found:
            prob_j = domains_found[domain_j]['Prob']
            interval_j = domains_found[domain_j]['Interval']

            if domain_i < domain_j:
                range_i = set(range(interval_i[0], interval_i[1]))
                range_j = set(range(interval_j[0], interval_j[1]))

                if len(range_i.intersection(range_j)) > 0:

                    if np.isnan(curr_prob_diffs):
                        curr_prob_diffs = abs(prob_i - prob_j)

                    if prob_i < prob_j:
                        domains_to_remove.append(domain_i)
                    else:
                        domains_to_remove.append(domain_j)

        domains_found[domain_i]['Prob_diff'] = curr_prob_diffs

    domains_accepted = {domain: domains_found[domain] for domain in domains_found if domain not in domains_to_remove}
    return domains_accepted

def add_domains_to_ncbi_code_data(ncbi_code_data, domains_found, correction_factor):

    for domain in domains_found:
        domain_interval = [domains_found[domain]['Interval'][0]+correction_factor,domains_found[domain]['Interval'][-1]+correction_factor]
        domain_label = domains_found[domain]['Domain_label']
        domain_probability = domains_found[domain]['Prob']
        domain_probability_diff = domains_found[domain]['Prob_diff']

        ncbi_code_data['annotated_intervals'].append(domain_interval)
        ncbi_code_data['annotated_domains'].append(domain_label)
        ncbi_code_data['annotated_probability'].append(domain_probability)
        ncbi_code_data['annotated_probability_diff'].append(domain_probability_diff)

    return ncbi_code_data

def organize_domain_order(ncbi_code_data):

    new_ncbi_code_data = {'seq': ncbi_code_data['seq'], 'annotated_domains': [], 'annotated_intervals': [], 'annotated_probability': [], 'annotated_probability_diff': []}

    intervals = []

    for i, domain in enumerate(ncbi_code_data['annotated_domains']):
        curr_interval = ncbi_code_data['annotated_intervals'][i]
        intervals.append(curr_interval)

    intervals = np.array(intervals)
    intervals.view('i8,i8').sort(order=['f1'], axis=0)
    intervals = [[int(value) for value in interval] for interval in intervals]

    for interval in intervals:
        interval_index_in_data = ncbi_code_data['annotated_intervals'].index(list(interval))
        domain_label = ncbi_code_data['annotated_domains'][interval_index_in_data]
        domain_probability = ncbi_code_data['annotated_probability'][interval_index_in_data]
        domain_probability_diff = ncbi_code_data['annotated_probability_diff'][interval_index_in_data]

        new_ncbi_code_data['annotated_domains'].append(domain_label)
        new_ncbi_code_data['annotated_probability'].append(domain_probability)
        new_ncbi_code_data['annotated_probability_diff'].append(domain_probability_diff)
        new_ncbi_code_data['annotated_intervals'].append(list(interval))

    return new_ncbi_code_data

def annotate_domains_in_stretches(ncbi_code_data, stretches_to_annotate, db, tmp_label = 'tmp', maxres = 50001, min_prob = args.min_prob, min_length = args.min_seq_length):

    for stretch in stretches_to_annotate:
        stretch_seq = stretches_to_annotate[stretch]['seq']
        correction_factor = stretches_to_annotate[stretch]['correction_factor']

        if len(stretch_seq) > min_length:
            # save sequence to file
            tmp_fasta = '{}{}_stretch_{}.fasta'.format(args.tmp_folder, tmp_label, stretch)
            with open(tmp_fasta, 'w') as fp:
                fp.write('>stretch_{}\n{}\n'.format(stretch, stretch_seq))

            # run hhsearch for this file
            out_hhr = '{}.hhr'.format(tmp_fasta)
            if not os.path.isfile(out_hhr):
                run_hhsearch = sp.Popen([hhsearch, '-i', tmp_fasta, '-d', db, '-o', out_hhr, '-norealign', '-maxres', str(maxres), '-p', str(min_prob), '-z', str(1)], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
                stdout, stderr = run_hhsearch.communicate()

            # get domains found
            domains_found = parse_hhsearch_out(out_hhr, min_hit_length = min_length)

            # find significant domains (remove overlapping matches)
            domains_found = find_significant_domains(domains_found)
                
            # update domain annotation for ncbi_code
            ncbi_code_data = add_domains_to_ncbi_code_data(ncbi_code_data, domains_found, correction_factor)

            # remove temporary files
            os.remove(tmp_fasta)
            os.remove(out_hhr)

    # reorganize the domain order
    ncbi_code_data = organize_domain_order(ncbi_code_data)
            
    return ncbi_code_data

## MAIN ROUTINES

def update_domain_composition(arguments, cpus = hhsearch_cpu, number_of_rounds = args.rounds, min_prob = args.min_prob, min_length = args.min_seq_length):

    jobID = arguments[0]
    ncbi_codes = arguments[1]
    all_sequences = arguments[2]
    db = arguments[3]
    
    print('\nANNOTATING SEQUENCES (thread {}):\n'.format(jobID+1))
    
    new_sequences = {}

    for i, ncbi_code in enumerate(ncbi_codes):      
        ncbi_code_data = all_sequences[ncbi_code]
        converged = False
        search_round = 0

        out_json = '{}{}_domains_annotated_{}rounds.json'.format(args.tmp_folder, ncbi_code, number_of_rounds)

        if not os.path.isfile(out_json):            
            while not converged:
                search_round += 1
                if search_round <= number_of_rounds:
                    print(' ... Thread {}: Working on {} (round {}) ({}/{})'.format(jobID+1, ncbi_code, search_round, i+1, len(ncbi_codes)))
                    
                    previous_number_of_domains = len(ncbi_code_data['annotated_domains'])
                    
                    stretches_to_annotate = find_sequence_stretches_to_annotate(ncbi_code_data)
                    ncbi_code_data = annotate_domains_in_stretches(ncbi_code_data, stretches_to_annotate, db, tmp_label = ncbi_code, min_prob = min_prob, min_length = min_length)

                    # check if it converged
                    number_new_domains = len(ncbi_code_data['annotated_domains']) - previous_number_of_domains
                    if number_new_domains == 0:
                        converged = True
                        print(' ... {}: Converged after {} rounds'.format(ncbi_code, search_round))
                else:
                    converged = True
                    print(' ... {}: Reached maximum number of rounds'.format(ncbi_code))

            json.dump(ncbi_code_data, open(out_json, 'w'), indent = 4)

        else:
            print('Already searched domains for {}. Will just load the data'.format(ncbi_code))
            ncbi_code_data = json.load(open(out_json, 'r'))
            
        new_sequences[ncbi_code] = ncbi_code_data
                    
    return new_sequences

def save_out_table(all_sequences, out_file):

    out = open(out_file, 'w')
    out.write('\nINPUT PARAMETERS:\n')
    out.write(' ... Database:                                 {}\n'.format(args.db_label))
    out.write(' ... Number of threads:                        {}\n'.format(args.n_threads))
    out.write(' ... CPUs for HHsearch:                        {}\n'.format(hhsearch_cpu))
    out.write(' ... Maximum no. of rounds per input sequence: {}\n'.format(args.rounds))
    out.write(' ... Minimum probability of a match:           {}\n'.format(args.min_prob))
    out.write(' ... Minimum sequence length to annotate:      {}\n\n'.format(args.min_seq_length))

    out.write('Sequence identifier\tDomain\tInterval\tProb.\tProb. different to next match{} label\tSequence\n'.format(args.db_label))
    for ncbi_code in all_sequences:
        out.write('\n')
        for i, domain in enumerate(all_sequences[ncbi_code]['annotated_domains']):
            interval = all_sequences[ncbi_code]['annotated_intervals'][i]
            seq = all_sequences[ncbi_code]['seq'][interval[0]-1:interval[-1]]
            prob = all_sequences[ncbi_code]['annotated_probability'][i]
            prob_diff = all_sequences[ncbi_code]['annotated_probability_diff'][i]
            out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ncbi_code, i+1, interval, prob, prob_diff, domain, seq))

    print('\nSAVING: Results saved in {}\n'.format(out_file))
    out.close()
      
## MAIN CODE

logfile = '{}/{}_domain_annotation.log'.format(out_dir, args.in_fasta.split('/')[-1].split('.')[0])
sys.stdout = Logger(logfile)

start = time.time()

# 1. Find the database and print its version
db = find_database(args.db_label)

# 2. Get sequences
all_sequences = parse_fasta(args.in_fasta)

# 3. Annotate domains in non-annotated regions in parallel
ncbi_codes = [ncbi_code for ncbi_code in all_sequences.keys()]
separated_jobs = chunk_list(ncbi_codes, args.n_threads)

list_arguments = [i for i in zip(range(args.n_threads), separated_jobs, [all_sequences for job in separated_jobs], [db for job in separated_jobs])]

start = time.time()
pool = mp.Pool(args.n_threads)
results = pool.map(update_domain_composition, list_arguments)
all_sequences = {key: dic[key] for dic in results for key in dic.keys()}

pool.close()
pool.join()

json.dump(all_sequences, open('{}/{}_domains_annotated.json'.format(out_dir, args.in_fasta.split('/')[-1].split('.')[0]), 'w'), indent = 4)

# 4. Save data into a table
out_file = '{}/{}_domains_annotated.csv'.format(out_dir, args.in_fasta.split('/')[-1].split('.')[0])
save_out_table(all_sequences, out_file)

end = time.time()
numb_seconds = end - start
print("END: Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
