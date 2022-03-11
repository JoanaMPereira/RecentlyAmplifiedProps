import matplotlib
matplotlib.use('Agg')
import json
import os
import time
import random
import urllib.parse
import urllib.request
import requests
import sys
import statistics
import argparse
import ast
import matplotlib.backends.backend_pdf

import subprocess as sp
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns

from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

# Get inputs
parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments with defaults')

# required inputs
requiredNamed.add_argument('-i', dest='assemblies_csv', type=str, required=True, help='csv file with the assemblies to which targets belong to')
# optional inputs
optionalNamed.add_argument('-n', dest='n_cases', default = 1, type=int, help='Number of cases to take (default: 1)')
optionalNamed.add_argument('-level', dest='level', default = 'Phylum', type=str, help='Level at which the cases should be taken  (default: Superkingdom)')
optionalNamed.add_argument('-cpu', dest='n_cpu', default = 1,type=int, help='Number of cpus to use (default: 1)')
optionalNamed.add_argument('-threads', dest='n_threads', default=1, type=int, help='Number of parallel jobs (default: 1)')
optionalNamed.add_argument('-rounds', dest='rounds', default=1, type=int, help='Number of GIAD rounds (default: 1)')
optionalNamed.add_argument('-tmp_folder', dest='tmp_folder', default = '/tmp/',type=str, help='Temporary folder (default: /tmp/)')

# Define inputs
args = parser.parse_args()

GIAD = '../propeller_validation_and_annotation/greedy_iterative_annotation_of_domains.py'
main_start = time.time()

if args.n_cases == -1:
	n_cases = 'all'
else:
	n_cases = args.n_cases

rounds = args.rounds

# HELPING ROUTINES

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

# 0. General helpers

def chunk_list(l, n):
    chunks = np.array_split(np.array(l), n)
    chunks = [list(chunk) for chunk in chunks]
    return chunks

# 1. To download files

def download_and_parse_refseq_and_gb_databases(databases = ['genbank', 'refseq']):
	
	database_assembly_mapping = {}

	for db in databases:
		print(' ... Taking care of {} summary table'.format(db))

		summary_table = '{}/assembly_summary_{}.txt'.format(os.getcwd(), db)

		if not os.path.isfile(summary_table):
			link = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/{}/assembly_summary_{}.txt'.format(db, db)

			print(' ... ... Downloading')
			download_db_table = sp.Popen(['wget', link], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
			stdout, stderr = download_db_table.communicate()
			print(' ... ... Done Downloading')

		else:
			print(' ... ... Summary table already exists (version from: {})'.format(time.ctime(os.path.getmtime(summary_table))))
		
		with open(summary_table, 'r') as summary:
			print(' ... ... Parsing summary table')
			for line in summary:
				if not line.startswith('#'):
					data = line.strip().split('\t')
					database_assembly_mapping[data[0]] = data[19]
		print(' ... ... Done parsing')
		
	return database_assembly_mapping

def download_and_extract_assembly(assembly_id, database_assembly_mapping, tmp_folder = args.tmp_folder, label = ''):

	print(' ... ... Downloading and extracting assembly {}'.format(assembly_id))

	if assembly_id in database_assembly_mapping:

		assembly_link = database_assembly_mapping[assembly_id]

		assembly_label = assembly_link.split('/')[-1]
		
		assembly_genomic_file_link = '{}/{}_genomic.gff.gz'.format(assembly_link, assembly_label)
		out_assembly_file_gz = '{}/{}'.format(tmp_folder, assembly_genomic_file_link.split('/')[-1])
		out_assembly_file = out_assembly_file_gz[:-3]

		if not os.path.isfile(out_assembly_file) and not os.path.isfile(out_assembly_file_gz):
			download_assembly = sp.Popen(['wget', assembly_genomic_file_link, '-O', out_assembly_file_gz], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
			stdout, stderr = download_assembly.communicate()

		if os.path.isfile(out_assembly_file_gz) or os.path.isfile(out_assembly_file):

			print(' ... ... ... Downloaded assembly {} ....'.format(assembly_id))

			extract_assembly_file = sp.Popen(['gunzip', out_assembly_file_gz], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
			stdout, stderr = extract_assembly_file.communicate()

			if os.path.isfile(out_assembly_file):

				assembly = parse_assembly(out_assembly_file)
				print(' ... ... ... Finished parsing assembly {} and found {} CDS entries'.format(assembly_id, len(assembly['ncbi_codes'])))
				return assembly, assembly_link

			else:
				print(' ... ... ... It was not possible to extract assembly {}'.format(assembly_id))
				return 'nan', assembly_link
		else:
			print(' ... ... ... It was not possible to save assembly {}'.format(assembly_id))
			return 'nan', assembly_link
	else:
		print(' ... ... ... Assembly {} is not in the database. It may have been replaced.'.format(assembly_id))
		return 'nan', 'nan'

def download_and_extract_assembly_seq(assembly_id, assembly_link, tmp_folder = args.tmp_folder, label = ''):

	print(' ... ... Downloading and extracting sequence for assembly {}'.format(assembly_id))

	assembly_label = assembly_link.split('/')[-1]
	assembly_seq_file_link = '{}/{}_genomic.fna.gz'.format(assembly_link, assembly_label)

	out_assembly_file_gz = '{}/{}'.format(tmp_folder, assembly_seq_file_link.split('/')[-1])
	out_assembly_file = out_assembly_file_gz[:-3]

	if not os.path.isfile(out_assembly_file) and not os.path.isfile(out_assembly_file_gz):
		download_assembly = sp.Popen(['wget', assembly_seq_file_link, '-O', out_assembly_file_gz], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
		stdout, stderr = download_assembly.communicate()

	if os.path.isfile(out_assembly_file_gz) or os.path.isfile(out_assembly_file):

		print(' ... ... ... Downloaded sequence for assembly {} ....'.format(assembly_id))

		extract_assembly_file = sp.Popen(['gunzip', out_assembly_file_gz], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
		stdout, stderr = extract_assembly_file.communicate()

		if os.path.isfile(out_assembly_file):

			assembly_seq = parse_genome_fasta(out_assembly_file)
			print(' ... ... ... Finished parsing sequence for assembly {} and found {} scaffolds'.format(assembly_id, len(assembly_seq)))
			return assembly_seq

		else:
			print(' ... ... ... It was not possible to extract sequence for assembly {}'.format(assembly_id))
			return 'nan'
	else:
		print(' ... ... ... It was not possible to save sequence for assembly {}'.format(assembly_id))
		return 'nan'

	return assembly_seq

# 2. Parsers

def select_target_assemblies(assemblies_df, n, level):

	if n == -1:
		n = 'all'

	print('\nSelecting {} top targets at level {}:\n'.format(n, level))

	if n != 'all':
		targets = assemblies_df.groupby(level).head(1).head(n)
	else:
		targets = assemblies_df

	targets = targets.reset_index(drop=True)

	print(targets)

	return targets

def parse_assembly(assembly_file):

	assembly = {'ncbi_codes': [], 'starts': [], 'ends': [], 'directions': [], 'names': [], 'scaffolds': []}
	curr_scaffold = 0
	
	with open(assembly_file, 'r') as in_assembly:
		for line in in_assembly:
			if not line.startswith('#'):
				line_data = line.split('\t')

				if line_data[2] == 'CDS':

					start = int(line_data[3])
					end = int(line_data[4])
					direction = line_data[6]

					if 'cds-' in line:
						ncbi_code = line_data[8].split('ID=cds-')[1].split(';')[0]
					else:
						if 'Name=' in line:
							ncbi_code = line_data[8].split('Name=')[1].split(';')[0]
						else:
							ncbi_code = 'unk'

					if 'pseudo=' not in line and 'product=' in line and 'fragment' not in line:
						prot_name = line_data[8].split('product=')[1].split(';')[0]
					else:
						prot_name = 'pseudogene'

					if ncbi_code in assembly['ncbi_codes'] and assembly['ncbi_codes'][-1] == ncbi_code: # it means that this is some kind of fragmented gene (has introns?...) and so we have to collect the largest interval encompassing it
						if start < assembly['starts'][-1]:
							assembly['starts'][-1] = start
						if end > assembly['ends'][-1]:
							assembly['ends'][-1] = end
					else:
						if '|' in ncbi_code:
							ncbi_code = ncbi_code.replace('|','_')

						assembly['ncbi_codes'].append(ncbi_code)
						assembly['names'].append(prot_name)
						assembly['scaffolds'].append(curr_scaffold)
						assembly['starts'].append(start)
						assembly['ends'].append(end)
						assembly['directions'].append(line_data[6])
					
			elif line.startswith('##sequence-region'):
				curr_scaffold = '{}'.format(line.split()[1])

	return assembly	

def parse_genome_fasta(assembly_seq):

	genome = {}

	scaffold_id = None
	sequence = ''
	with open(assembly_seq, 'r') as assemblytxt:
		for line in assemblytxt:
			if line.startswith('>'):
				if scaffold_id != None:
					genome[scaffold_id] = sequence

				scaffold_id = line.split()[0].replace('>', '')
				sequence = ''
			else:
				sequence += line.strip().upper()
	genome[scaffold_id] = sequence.upper()

	return genome

# 3. To find, extract and process intergenic regions

def get_reverse_complemente(sequence, code = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}):

	complement = ''

	for nt in sequence:
		complement += code[nt]

	reverse_complement = ''.join(reversed(complement))

	return reverse_complement


def get_intergenic_regions_next_to_targets(target_proteins, assembly_gff, assembly_seq, assembly_id):

	intergenic_regions = []

	print(' ... ... Collecting intergenic region close to the targets in {}. Pseudogenes will be considered as intergenic spaces'.format(assembly_id))

	targets_not_found = []
	for protein in target_proteins:
		if protein in assembly_gff['ncbi_codes']:
			target_index = assembly_gff['ncbi_codes'].index(protein)

			protein_direction = assembly_gff['directions'][target_index]
			protein_scaffold = assembly_gff['scaffolds'][target_index]

			neighbours_indeces = [j for j, x in enumerate(assembly_gff['scaffolds']) if x == protein_scaffold and assembly_gff['directions'][j] == protein_direction and assembly_gff['names'][j] != 'pseudogene']
			
			if len(neighbours_indeces) > 1:
				if neighbours_indeces.index(target_index) > 0 and neighbours_indeces.index(target_index) < len(neighbours_indeces)-1:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)-1], neighbours_indeces[neighbours_indeces.index(target_index)], neighbours_indeces[neighbours_indeces.index(target_index)+1]]
				elif neighbours_indeces.index(target_index) == 0:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)], neighbours_indeces[neighbours_indeces.index(target_index)+1]] 
				else:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)-1], neighbours_indeces[neighbours_indeces.index(target_index)]]

				# first add the spaces between genes
				for j, curr_gene in enumerate(genomic_context[:-1]):
					next_gene = genomic_context[j+1]

					curr_end = assembly_gff['ends'][curr_gene]
					next_start = assembly_gff['starts'][next_gene]

					if next_start - curr_end > 0: # if there is no overlap/if they are not transcriptionally coupled
						intergenic_region = '{} : {}-{} ({})'.format(protein_scaffold, curr_end+1, next_start, protein_direction)
						intergenic_regions.append(intergenic_region)

				# now, if our target only has a neighbour, it means it is either at the start or the end of the assembly. check if there is some 
				# genomic sequence between it and the closest terminus of the genome
				if genomic_context[0] == target_index:
					curr_start = assembly_gff['starts'][target_index]
					if curr_start > 1:
						intergenic_region = '{} : {}-{} ({})'.format(protein_scaffold, 1, curr_start, protein_direction)
						intergenic_regions.append(intergenic_region)
				elif genomic_context[-1] == target_index:
					curr_end = assembly_gff['ends'][target_index]
					if curr_end < len(assembly_seq[protein_scaffold]):
						intergenic_region = '{} : {}-{} ({})'.format(protein_scaffold, curr_end+1, len(assembly_seq[protein_scaffold]), protein_direction)
						intergenic_regions.append(intergenic_region)
			else:				
				intergenic_region = '{} : {}-{} ({})'.format(protein_scaffold, 1, assembly_gff['starts'][target_index], protein_direction)
				intergenic_regions.append(intergenic_region)

				intergenic_region = '{} : {}-{} ({})'.format(protein_scaffold, assembly_gff['ends'][target_index]+1, len(assembly_seq[protein_scaffold]), protein_direction)
				intergenic_regions.append(intergenic_region)
			
		else:
			targets_not_found.append(protein)

	intergenic_regions = sorted(list(set(intergenic_regions)))

	# if len(targets_not_found) == 0:
	print(' ... ... ... Found {} intergenic_regions to analyse in assembly {}.'.format(len(intergenic_regions), assembly_id))
	# else:
	# 	print(' ... ... ... Found {} intergenic_regions to analyse in assembly {}. {} targets where not found, which means they may have been removed and may not represent genes anymore: {}'.format(len(intergenic_regions), assembly_id, len(targets_not_found), targets_not_found))

	return intergenic_regions, assembly_gff

def collect_orfs(intergenic_regions, assembly_seq, assembly_id, min_len = 20):

	print(' ... ... Collection all protein translations of intergenic region close to the targets in {}'.format(assembly_id))

	orfs = {}

	for intergenic_region in intergenic_regions:
		scaffold = intergenic_region.split()[0]
		start = int(intergenic_region.split()[2].split('-')[0])
		end = int(intergenic_region.split()[2].split('-')[1])
		direction = intergenic_region.split()[-1].replace('(','').replace(')','')

		sequence = assembly_seq[scaffold][start-1:end]

		if direction == '-':
			#sequence = get_reverse_complemente(sequence)
			sequence = str(Seq(sequence, IUPAC.ambiguous_dna).reverse_complement())

		for i in range(3):
			curr_sequence = sequence[i:]
			if len(curr_sequence) % 3 != 0:
				curr_sequence = curr_sequence[:-(len(curr_sequence) % 3)]

			orf_seq = str(Seq(curr_sequence, IUPAC.ambiguous_dna).translate())
			orf_seq = orf_seq.replace('*', 'X')

			if len(orf_seq) >= min_len:
				orf_name = '{}:{}-{}({})'.format(scaffold, start+i, start+i+((len(orf_seq)-1)*3), direction)
				orfs[orf_name] = orf_seq
	
	print(' ... ... ... Found {} protein-like sequences in assembly {}'.format(len(orfs), assembly_id))

	return orfs


# MAIN ROUTINES

def collect_assemblies_and_intergenic_orfs(arguments):

	jobID = arguments[0]
	target_assemblies = arguments[1]
	target_assemblies_df = arguments[2]
	database_assembly_mapping = arguments[3]

	target_assemblies_data = {assembly_id: {'intergenic_orfs': '', 'relevant_scaffolds_summary': {}} for assembly_id in target_assemblies} 

	print('\nANALYSING GENOMES (thread {}):\n'.format(jobID+1))


	for i, assembly_id in enumerate(target_assemblies):

		curr_tax_label = list(target_assemblies_df.loc[target_assemblies_df.Assemblies == assembly_id].Phylum_label)[0]
		curr_species = list(target_assemblies_df.loc[target_assemblies_df.Assemblies == assembly_id].Species)[0]

		print('\n ... Thread {}: Working on {} from {} {} ({}/{})'.format(jobID+1, assembly_id, curr_tax_label, curr_species, i+1, len(target_assemblies)))

		assembly_gff, assembly_link = download_and_extract_assembly(assembly_id, database_assembly_mapping, tmp_folder = args.tmp_folder, label = assembly_id)

		if assembly_link != 'nan' and assembly_gff != 'nan':

			assembly_seq = download_and_extract_assembly_seq(assembly_id, assembly_link, tmp_folder = args.tmp_folder, label = assembly_id)

			if  assembly_seq != 'nan':

				curr_target_proteins = ast.literal_eval(list(target_assemblies_df.loc[target_assemblies_df.Assemblies == assembly_id].Members)[0])
				
				# for each target protein, get the intergenic regions between it and the two flanking genes
				intergenic_regions, assembly_gff = get_intergenic_regions_next_to_targets(curr_target_proteins, assembly_gff, assembly_seq, assembly_id)

				# collect open reading frames in intergenic regions
				orfs = collect_orfs(intergenic_regions, assembly_seq, assembly_id)

				target_assemblies_data[assembly_id]['intergenic_orfs'] = orfs

				scaffolds_data = {}
				for i, scaffold_id in enumerate(assembly_gff['scaffolds']):
					if scaffold_id not in scaffolds_data:
						scaffolds_data[scaffold_id] = {'length': len(assembly_seq[scaffold_id]),'ncbi_codes': [], 'starts': [], 'ends': [], 'directions': [], 'names': []}

					for key in scaffolds_data[scaffold_id].keys():
						if key != 'length':
							scaffolds_data[scaffold_id][key].append(assembly_gff[key][i])

				# save only the scaffolds that have any of the target proteins
				for scaffold_id in scaffolds_data:
					if any(target in scaffolds_data[scaffold_id]['ncbi_codes'] for target in curr_target_proteins):
						target_assemblies_data[assembly_id]['relevant_scaffolds_summary'][scaffold_id] = scaffolds_data[scaffold_id]

	return target_assemblies_data

def annotate_domains_and_find_propellers(target_assemblies_data, rounds = 1, max_x_content = 15, select_min_prob = 50, GIAD_min_prob = 0, min_len = 20, cpu = args.n_cpu, threads = args.n_cpu):

	prot_count = 0

	tmp_fasta = '{}/all_intergenic_orfs_{}_{}_{}.fasta'.format(args.tmp_folder, args.assemblies_csv.split('/')[-1].split('.')[0], n_cases, args.level)
	with open(tmp_fasta, 'w') as out_fasta:
		for assembly_id in target_assemblies_data:
			for orf in target_assemblies_data[assembly_id]['intergenic_orfs']:
				out_fasta.write('>Assembly-{}:{}\n{}\n'.format(assembly_id, orf, target_assemblies_data[assembly_id]['intergenic_orfs'][orf]))
				prot_count += 1

	print('\nFINDING GRAVEYARDS in {} intergenic protein-like proteins around the targets in ALL assemblies (may take some time)'.format(prot_count))

	print('\n ... Parameters used:')
	print(' ... ... GIAD rounds:       {}'.format(rounds))
	print(' ... ... GIAD min. Prob.:   {}'.format(GIAD_min_prob))
	print(' ... ... Minimum length:    {}\n'.format(min_len))
	print(' ... ... Max_*_content:     {}'.format(max_x_content))
	print(' ... ... Select min. Prob.: {}'.format(select_min_prob))


	#run GIAD on ECOD, collecting all hits even at low probability

	graveyards = {}
	out_annotated = '{}/{}_domains_annotated.json'.format(args.tmp_folder, tmp_fasta.split('/')[-1].split('.')[0])

	if not os.path.isfile(out_annotated):
		print(' ... ... Annotating domains')
		annotate_domains = sp.Popen(['nice', 'python3', GIAD, '-i', tmp_fasta, '-db', 'ECOD', '-rounds', str(rounds), '-min_prob', str(GIAD_min_prob), '-min_len', str(min_len), '-cpu', str(cpu), '-threads', str(threads), '-tmp_folder', args.tmp_folder, '-out_dir', args.tmp_folder], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
		stdout, stderr = annotate_domains.communicate()

	if os.path.isfile(out_annotated):
		print(' ... ... Collecting all beta-propeller matches')

		graveyard_count = 0

		out_annotated = json.load(open(out_annotated, 'r'))

		for orf in out_annotated:
			for i, domain in enumerate(out_annotated[orf]['annotated_domains']):
				if 'beta-propeller' in domain:
					probability = out_annotated[orf]['annotated_probability'][i]
					interval = out_annotated[orf]['annotated_intervals'][i]
					seq = out_annotated[orf]['seq'][interval[0]-1:interval[-1]] 
					x_content = sum([seq.count('X')])*100/len(seq)

					if probability >= select_min_prob and x_content <= max_x_content:
						assembly_id = orf.split(':')[0].split('Assembly-')[-1]
						scaffold_id = orf.split(':')[1]
						orf_direction = orf.split('(')[-1].replace(')','')

						orf_position_in_scaffold = orf.split(':')[2].split('(')[0].split('-')
						orf_position_in_scaffold = [int(pos) for pos in orf_position_in_scaffold]

						grave_start = orf_position_in_scaffold[0]+(interval[0]*3)-3
						grave_end = orf_position_in_scaffold[0]+(interval[1]*3)

						if assembly_id not in graveyards:
							graveyards[assembly_id] = {}
						if scaffold_id not in graveyards[assembly_id]:
							graveyards[assembly_id][scaffold_id] = {'starts': [], 'ends': [], 'directions': [], 'probabilities': [], 'orf_names': [], 'graveyard_names': []}

						graveyards[assembly_id][scaffold_id]['starts'].append(grave_start)
						graveyards[assembly_id][scaffold_id]['ends'].append(grave_end)
						graveyards[assembly_id][scaffold_id]['directions'].append(orf_direction)
						graveyards[assembly_id][scaffold_id]['probabilities'].append(probability)
						graveyards[assembly_id][scaffold_id]['orf_names'].append(orf)
						graveyards[assembly_id][scaffold_id]['graveyard_names'].append('graveyard_{}'.format(graveyard_count))

						graveyard_count += 1

		print(' ... ... ... Found {} graveyards in {} assemblies'.format(graveyard_count, len(graveyards)))

	else:
		print(' ERROR: Not possible to annotate domains')
		print(stderr)

	return graveyards

def save_summary_table(graveyards, target_assemblies_df):

	out_table = '{}_{}_{}_target_with_graveyard_count.csv'.format(args.assemblies_csv.split('/')[-1].split('.')[0], n_cases, args.level)

	graveyard_counts = []

	for index, row in target_assemblies_df.iterrows():
		count = 0
		if row.Assemblies in graveyards:
			for scaffold_id in graveyards[row.Assemblies]:
				count += len(graveyards[row.Assemblies][scaffold_id]['starts'])
		graveyard_counts.append(count)

	target_assemblies_df['graveyard_counts'] = graveyard_counts
	target_assemblies_df.to_csv(out_table, index=False, sep='\t')
	print(target_assemblies_df.head())


# MAIN CODE

logfile = '{}_graveyards_search_{}_{}.log'.format(args.assemblies_csv.split('/')[-1].split('.')[0], n_cases, args.level)
sys.stdout = Logger(logfile)

start = time.time()

# Download and parse RefSeq and Genbank databases
print('\nDownloading and parsing RefSeq and Genbank summary tables')
refseq_gb_assembly_map = download_and_parse_refseq_and_gb_databases()

# select the target assemblies 
assemblies_df = pd.read_csv(args.assemblies_csv, sep = '\t')
target_assemblies_df = select_target_assemblies(assemblies_df, n = n_cases, level = args.level)

# separate the target assemblies through the multiple threads and collect their assembly information and the intergenic 'orfs' between the targets
individual_target_assemblies = [row.Assemblies for index, row in target_assemblies_df.iterrows()]

out_json = '{}/{}_{}_{}_target_assemblies_all_data.json'.format(args.tmp_folder, args.assemblies_csv.split('/')[-1].split('.')[0], n_cases, args.level)

if not os.path.isfile(out_json):

	if len(individual_target_assemblies) < args.n_threads:
		n_threads = len(individual_target_assemblies)
	else:
		n_threads = args.n_threads

	print('\nStarting {} threads'.format(n_threads))

	separated_jobs = chunk_list(individual_target_assemblies, n_threads)

	list_arguments = [i for i in zip(range(args.n_threads), separated_jobs, [target_assemblies_df for job in separated_jobs], [refseq_gb_assembly_map for job in separated_jobs])]

	start = time.time()
	pool = mp.Pool(args.n_threads)
	target_assemblies_data = pool.map(collect_assemblies_and_intergenic_orfs, list_arguments)
	target_assemblies_data = {key: dic[key] for dic in target_assemblies_data for key in dic.keys()}

	pool.close()
	pool.join()

	json.dump(target_assemblies_data, open(out_json, 'w'), indent = 4)

else:
	print('\nGenomic data was already collected. Will read it and continue.')
	target_assemblies_data = json.load(open(out_json, 'r'))

end = time.time()
numb_seconds = end - start
print(" ... ANALYSIS OF GENOMES: Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))


# now annotate all these sequences with GIAD
start = time.time()

out_graveyards_json = '{}_{}_{}_target_graveyards.json'.format(args.assemblies_csv.split('/')[-1].split('.')[0], n_cases, args.level)

if not os.path.isfile(out_graveyards_json):
	graveyards = annotate_domains_and_find_propellers(target_assemblies_data, rounds = rounds, max_x_content = 15, select_min_prob = 1, GIAD_min_prob = 1, min_len = 20, cpu = args.n_cpu, threads = args.n_cpu)
	json.dump(graveyards, open(out_graveyards_json, 'w'), indent = 4)
else:
	print('\nGraveyards were already found. Will read them and continue.')
	graveyards = json.load(open(out_graveyards_json, 'r'))

end = time.time()
numb_seconds = end - start
print(" ... FINDING GRAVEYARDS: Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

# now add a column with the number of graveyards in the target genomes and save table
save_summary_table(graveyards, target_assemblies_df)

main_end = time.time()
numb_seconds = main_end - main_start
print("\nEND: Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))


