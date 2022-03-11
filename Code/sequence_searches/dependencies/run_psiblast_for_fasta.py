import subprocess as sp
import sys 
import os
import os.path
import json
import re

from Bio.Blast import NCBIXML
from Bio import Entrez
from numpy import *
from math import *
from time import time

# define databases location
seq_database = '/ebio/abt1_share/toolkit_sync/databases/standard/'
# define programs location
blast = '/ebio/abt1_share/toolkit_sync/bioprogs/tools/blastplus/bin/' 

# define tmp folder
#tmp_folder = '/tmp/'
tmp_folder = '/tmp/jpereira/'

Entrez.email = 'joana.pereira@tuebingen.mpg.de'

# define parameters

#psiblast parameters for blades search
seg_on = 'yes'
max_evalue = 1
inclusion_ethresh = 1e-4
num_alignments = 50000
hit_cov = 80.0
max_overlap_length = 10
exclude_highly_repetitive = False
get_only_strains = False 

seq_fasta = sys.argv[1]
mode = sys.argv[2]	# avaliable modes: 'aligned_only' and 'full_sequence', or 'both'
database = sys.argv[3]
num_iterations = sys.argv[4] # until convergence is 0
collect_rounds = int(sys.argv[5])
num_threads = sys.argv[6]
working_directory = sys.argv[7]

print "\nGoing to run Psiblast for file {} over {} database with following parameters:".format(seq_fasta, database)
print " ... Hit dispalying E-value: {}".format(max_evalue)
print " ... Using segmasker?      : {}".format(seg_on)
print " ... PSSM inclusion E-value: {}".format(inclusion_ethresh)
print " ... Number of rounds:       {}".format(num_iterations)
print " ... Rounds collected:       {}".format(collect_rounds)
print " ... Maximum number of alig: {}".format(num_alignments)
print " ... Minimum query coverag.: {}".format(hit_cov)
print "\nAnd will save the hit sequences based on:"
print " ... Mode:                   {}".format(mode)
print " ... Maximum overlap:        {}".format(max_overlap_length)
print " ... Exclude repetitive:     {}".format(exclude_highly_repetitive)

# subroutines

def process_sequence_title(title):

    try:
        ncbi_code = title.split('|')[1]
        species = title.split('[')[-1].split(']')[0]
        prot_name = title.split('|')[2].split('[')[0]
    except:
        try:
            ncbi_code = title.split(' ')[0]
            species = title.split('[')[-1].split(']')[0]
            prot_name = title.replace('{} '.format(ncbi_code),'')
            prot_name = prot_name.replace(' [{}]'.format(species),'')
        except:
            ncbi_code = 'not available'
            species   = 'not available'
            prot_name = 'not available'

            print "   ERROR: couldn't take care of the results for {}.".format(title)
            print "   .....: something is wrong with the title."

    return ncbi_code, species, prot_name

def find_repeats(target_sequence):

    r = re.compile(r"(.+?)\1+")
    for match in r.finditer(target_sequence):
        yield (match.group(1), len(match.group(0))/len(match.group(1)))
           
def compute_repeat_level(target_sequence):

    target_sequence = target_sequence.replace('-','')  
    repeats = list(find_repeats(target_sequence))
  
    repeat_level = 0
    for repeat in repeats:
        seq_repeat = repeat[0]
        count_repeat = repeat[1]

        repeat_level += len(seq_repeat)*count_repeat

    repeat_level = repeat_level*100.0/len(target_sequence)
    repeat_lengths = [len(repeat[0]) for repeat in repeats]
    if len(repeat_lengths) > 0:
        max_repeat_length = max(repeat_lengths)
    else:
        max_repeat_length = 0

    return repeat_level, max_repeat_length

# main routine

def run_psiblast_and_mmseqs2(seq_fasta, mode, get_only_strains = get_only_strains, exclude_highly_repetitive = exclude_highly_repetitive, hit_cov = hit_cov, num_threads = num_threads, database = database, seg_on = seg_on, max_evalue = max_evalue, inclusion_ethresh = inclusion_ethresh, num_alignments = num_alignments):

    # first psiblast all sequences against the target database
    
    blast_outfile = '{}{}_psiblast_{}rounds_{}_max_eval_{}.xml'.format(tmp_folder, seq_fasta.split('/')[-1], num_iterations, database, max_evalue)
    #blast_outfile = '{}_psiblast_{}rounds_{}_max_eval_{}.xml'.format(seq_fasta, num_iterations, database, max_evalue)

    print "1. Psiblasting {} over {} database".format(seq_fasta, database)
    ran_psiblast = False
    if os.path.isfile(blast_outfile):
        print " ... {} file already exists. Will only parse it".format(blast_outfile)
    else:
        run_blast = sp.Popen(['{}psiblast'.format(blast), '-query', seq_fasta, '-db', '{}{}'.format(seq_database, database), '-out', blast_outfile, '-seg', seg_on, '-num_threads', str(num_threads), '-num_iterations', str(num_iterations),'-num_alignments', str(num_alignments), '-evalue', str(max_evalue), '-inclusion_ethresh', str(inclusion_ethresh), '-outfmt', '5'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_blast.communicate()
        ran_psiblast = True

        print stdout
        print stderr

    # parse the blast outfile with biopython and save all sequences to a fasta file

    print "\n2. Parsing Blast output file"
    
    result_handle = open(blast_outfile)
    blast_records = NCBIXML.parse(result_handle)
    #blast_records = list(blast_records)

    hit_sequences = {}
    previous_query = ''
    
    print " 2.1. Going through each record"
    
    for record in blast_records:

        query_name = record.query
        query_total_length = record.query_length

        if query_name != previous_query:
            collected = 0
        
        if collected < collect_rounds:

            print " ... Processing blast results for {}".format(query_name)

            num_sequences = len(record.alignments)

            print " ... ... There are {} sequences found".format(num_sequences)

            for alignment in record.alignments:
                title = alignment.title
                try:
                    ncbi_code, species, prot_name = process_sequence_title(title)
                    
                    if get_only_strains:
                        if len(species.split()) > 2:
                            species_is_strain = True
                        else:
                            species_is_strain = False
                    
                    if not get_only_strains or (get_only_strains and species_is_strain):
                    
                        for hsp in alignment.hsps:
                        
                            query_coverage = (hsp.query_end - hsp.query_start)*100.0/float(query_total_length)
                            match_query_coverage = (hsp.sbjct_end - hsp.sbjct_start)*100.0/float(query_total_length)
                        
                            if query_coverage >= hit_cov: 
            
                                target_sequence = hsp.sbjct

                                # check if match is highly repetitive
                                if exclude_highly_repetitive:
                                    repetition_level, max_repeat_length = compute_repeat_level(target_sequence)
                                else:
                                    repetition_level = 0
                                    max_repeat_length = 0
                                   
                                if repetition_level <= 20.0 and max_repeat_length <= 4:

                                    if ncbi_code not in hit_sequences.keys():
                                        hit_sequences[ncbi_code] = {'Interval': [],
                                                                    'Query_name': [],
                                                                    'E-value': [],
                                                                    'Species': species,
                                                                    'Prot_name': prot_name}
                                
                                    query_sequence = hsp.query
                                    evalue = hsp.expect
                                    interval = [hsp.sbjct_start, hsp.sbjct_end]

                                    if interval not in hit_sequences[ncbi_code]['Interval']:

                                        hit_sequences[ncbi_code]['Interval'].append(interval)
                                        hit_sequences[ncbi_code]['Query_name'].append(query_name)
                                        hit_sequences[ncbi_code]['E-value'].append(evalue)
                                
                                    else:
                                        index_interval = hit_sequences[ncbi_code]['Interval'].index(interval)   

                                        if hit_sequences[ncbi_code]['E-value'][index_interval] < evalue:

                                            hit_sequences[ncbi_code]['Interval'][index_interval] = interval
                                            hit_sequences[ncbi_code]['Query_name'][index_interval] = query_name
                                            hit_sequences[ncbi_code]['E-value'][index_interval] = evalue         


                except:
                    print "   ERROR: couldn't take care of the results for {}:".format(title)

            print ' ... ... ... {} sequences collected so far'.format(len(hit_sequences.keys()))

            outputfile_hit_sequences_name = '{}/{}_psiblast_{}collected_rounds_{}_max_eval_{}_unique_hit_sequences.json'.format(working_directory, seq_fasta.split('/')[-1], collect_rounds, database, max_evalue)
            outputfile_hit_sequences = open(outputfile_hit_sequences_name, "w")
            json.dump(hit_sequences, outputfile_hit_sequences, indent = 4)

            previous_query = query_name
            collected += 1
            
    #if ran_psiblast:
        # try to copy the blast file to the working directory
        #copy_blast = sp.Popen(['cp', blast_outfile, working_directory], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        #stdout, stderr = copy_blast.communicate()
        #print(stderr)

    return hit_sequences, blast_outfile

# ROUTINES TO SELECT THE BEST INTERVALS IN EACH UNIQUE MATCHED NCBI ID

def get_non_overlapping_intervals(intervals, evalues, queries, max_evalue = max_evalue, max_overlap = max_overlap_length):

    indeces_to_remove = []
    
    for i in range(len(intervals)):
        for j in range(len(intervals)):
            if i < j and i != j:

                if evalues[i] <= max_evalue and evalues[j] <= max_evalue:
                    range_i = set(range(intervals[i][0], intervals[i][1]))
                    range_j = set(range(intervals[j][0], intervals[j][1]))

                    if len(range_i.intersection(range_j)) > max_overlap:
                        if evalues[j] < evalues[i]:
                            indeces_to_remove.append(i)
                        else:
                            indeces_to_remove.append(j)

                elif evalues[i] > max_collect_evalue:
                    indeces_to_remove.append(i)
                elif evalues[j] > max_collect_evalue:
                    indeces_to_remove.append(j)

    new_intervals = [intervals[i] for i in range(len(intervals)) if i not in indeces_to_remove]
    new_evalues = [evalues[i] for i in range(len(evalues)) if i not in indeces_to_remove]
    new_queries = [queries[i] for i in range(len(queries)) if i not in indeces_to_remove]
    
    return new_intervals, new_evalues, new_queries

def find_best_non_overlapping_matches(hit_sequences):

    print "Finding the best match intervals for each subject sequence"
    new_hit_sequences = {}

    count = 0
    for ncbi_code in hit_sequences.keys():

        count+=1
        print " ... {} ({}/{} = {}%)".format(ncbi_code, count, len(list(hit_sequences.keys())), round(float(count*100/len(list(hit_sequences.keys()))), 2))
        
        species = hit_sequences[ncbi_code]['Species']
        prot_name = hit_sequences[ncbi_code]['Prot_name']
        
        intervals = hit_sequences[ncbi_code]['Interval']
        evalues = hit_sequences[ncbi_code]['E-value']
        queries = hit_sequences[ncbi_code]['Query_name']
        new_intervals, new_evalues, new_queries = get_non_overlapping_intervals(intervals, evalues, queries)

        if len(new_intervals) > 1: 
            print " ... Subject {} has {} non-overlapping matches".format(ncbi_code, len(new_intervals))
            
        new_hit_sequences[ncbi_code] = {'Interval': new_intervals,
                                        'Query_name': new_queries,
                                        'E-value': new_evalues,
                                        'Species': species,
                                        'Prot_name': prot_name}

    outputfile_hit_sequences_name = '{}/{}_psiblast_{}collected_rounds_{}_max_eval_{}_unique_hit_sequences_processed.json'.format(working_directory, seq_fasta.split('/')[-1], collect_rounds, database, max_evalue)
    outputfile_hit_sequences = open(outputfile_hit_sequences_name, "w")
    json.dump(new_hit_sequences, outputfile_hit_sequences, indent = 4)
    
    return new_hit_sequences

def  write_hits_to_fasta(hit_sequences, mode):

    print "Writting hit sequences to file"
    
    enriched_file = '{}/{}_psiblast_{}collected_rounds_{}_enriched_set_{}_max_eval_{}.fasta'.format(working_directory, seq_fasta.split('/')[-1][:-6], collect_rounds, database, mode, max_evalue)
    enriched_out = open(enriched_file, 'w')
    
    written_sequences = 0
    for ncbi_code in hit_sequences.keys():

        #print " ... Searching for {} on nr database".format(ncbi_code)
                            
        # get protein from nr database
        run_blastdb = sp.Popen(['{}blastdbcmd'.format(blast), '-entry', ncbi_code, '-db', '{}nr'.format(seq_database)], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_blastdb.communicate()
        title, target_sequence = stdout.split('\n')[0], ''.join(stdout.split('\n')[1:])

        species = hit_sequences[ncbi_code]['Species']
        prot_name = hit_sequences[ncbi_code]['Prot_name']
	
	if len(target_sequence) > 0:
            if mode == 'aligned_only':
                intervals = array(hit_sequences[ncbi_code]['Interval'])
                intervals.sort(axis = 0)

                for i, interval in enumerate(intervals):
                    evalue = hit_sequences[ncbi_code]['E-value'][i]
	       	    coverage = float(interval[-1]-interval[0])/len(target_sequence)

                    prot_identifier = '>{}_#{}"{}"|{}|{}| E-val: {}| Interval: {}| Full_seq_coverage: {}%'.format(ncbi_code, i+1, len(intervals), species, prot_name, evalue, interval, round(coverage*100, 1))
                    match_sequence = target_sequence[interval[0]-1:interval[-1]]

                    enriched_out.write('{}\n'.format(prot_identifier))
                    enriched_out.write('{}\n'.format(match_sequence))
		    written_sequences += 1

            elif mode == 'full_sequence':

                prot_identifier = '>{}|{}|{}| E-vals: {}| Intervals: {}'.format(ncbi_code, species, prot_name, hit_sequences[ncbi_code]['E-value'], hit_sequences[ncbi_code]['Interval'])
                enriched_out.write('{}\n'.format(prot_identifier))
                enriched_out.write('{}\n'.format(target_sequence))  
		written_sequences += 1
	    
            else:
		print "Mode '{}' not recognised. Use 'aligned_only' or 'full_sequence'"

    enriched_out.close()
    
    print " ... Wrote {} sequences to file {}".format(written_sequences, enriched_file)
    
    return enriched_file

# run main routine
hit_sequences, blast_outfile = run_psiblast_and_mmseqs2(seq_fasta, mode, exclude_highly_repetitive = exclude_highly_repetitive, hit_cov = hit_cov, num_threads = num_threads, database = database, max_evalue = max_evalue, num_alignments = num_alignments)
hit_sequences = find_best_non_overlapping_matches(hit_sequences)

if mode == 'both':
    out_fasta = write_hits_to_fasta(hit_sequences, 'aligned_only')
    out_fasta = write_hits_to_fasta(hit_sequences, 'full_sequence')
else:
    out_fasta = write_hits_to_fasta(hit_sequences, mode)
