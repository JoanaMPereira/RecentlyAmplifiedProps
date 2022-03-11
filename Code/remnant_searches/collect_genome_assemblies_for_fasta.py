import argparse
import time
import os
import sys
import ast
import json
import urllib.parse
import urllib.request
import requests

import subprocess as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Bio.Blast import NCBIXML
from Bio import Entrez


# Get inputs
parser = argparse.ArgumentParser(description='')
requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments with defaults')

# required inputs
requiredNamed.add_argument('-i', dest='in_file', type=str, required=True, help='input file with sequences to map to genomic assembly')
requiredNamed.add_argument('-mode', dest='mode', choices=['fasta', 'taxonomy'], type=str, required=True, help='mode/format of input file. It works with fasta files or taxonomy dictionaries')
# optional inputs
optionalNamed.add_argument('-tmp_folder', dest='tmp_folder', default = '/tmp/',type=str, help='Temporary folder (default: /tmp/)')
optionalNamed.add_argument('-user_email', dest='user_email', default = 'joana.pereira@tuebingen.mpg.de',type=str)

args = parser.parse_args()
curr_directory = os.getcwd()
Entrez.email = args.user_email

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

# HELPING ROUTINES

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
					assembly = data[0]
					link = data[19]
					if 'GCA' in assembly:
						refseq = data[17]
					else:
						refseq = assembly

					database_assembly_mapping[assembly] = {'link': link, 'refseq': refseq}

		print(' ... ... Done parsing')
		
	return database_assembly_mapping

def parse_targets(in_file, mode = args.mode):

	print(' ... Input file: {}'.format(in_file))
	print(' ... Input mode: {}'.format(mode))

	targets_list = {}
	total_count = 0

	if mode == 'fasta':
		with open(in_file, 'r') as curr_infile:
			for line in curr_infile:
				if line.startswith('>'):
					is_fasta = True
					curr_target = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()

					if curr_target not in targets_list:
						targets_list[curr_target] = {}
						total_count += 1

	elif mode == 'taxonomy':
		taxonomy = json.load(open(in_file, 'r'))
		for superkingdom in taxonomy:
			for phylum in taxonomy[superkingdom]:
				for taxclass in taxonomy[superkingdom][phylum]:
					for order in taxonomy[superkingdom][phylum][taxclass]:
						for genus in taxonomy[superkingdom][phylum][taxclass][order]:
							for species in taxonomy[superkingdom][phylum][taxclass][order][genus]:
								for curr_target in taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes']:
									if curr_target not in targets_list:
										targets_list[curr_target] = {'Superkingdom': superkingdom,
																	 'Phylum': phylum,
																	 'Class': taxclass,
																	 'Order': order,
																	 'Genus': genus,
																	 'Species': species}
										total_count += 1

	print(' ... Found {} unique ncbi codes out of {} targets'.format(len(targets_list), total_count))
	
	return targets_list

def map_uniprot_to_ncbi(uniprot_code, search_database = 'P_REFSEQ_AC'):

	if 'UniRef' in uniprot_code:
		uniprot_label = uniprot_code.split('_')[-1]
	else:
		uniprot_label = uniprot_code

	uniprot_url = 'https://www.uniprot.org/uploadlists/'
	
	params = {'from': 'ACC+ID','to': search_database,'format': 'tab','query': uniprot_label}
	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')

	ncbi_code = 'nan'

	try:
		uniprot_req = urllib.request.Request(uniprot_url, data)
		with urllib.request.urlopen(uniprot_req) as f:
			response = f.read()
		
		for line in response.decode('utf-8').split('\n'):
			if 'From' not in line and ncbi_code == 'nan':
				if len(line) > 1:
					ncbi_code = line.split('\t')[-1]

		if ncbi_code == 'nan':
			if search_database != 'EMBL':
				ncbi_code =  map_uniprot_to_ncbi(uniprot_code, search_database = 'EMBL')
				
		else:
			print(" ... {} corresponds to {} in ncbi_database".format(uniprot_code, ncbi_code))
	except:
		ncbi_code = 'nan'

	return ncbi_code

def find_ncbi_code_assembly(ncbi_code, database_assembly_mapping):

	assembly_id = 'nan'
	assembly_source = 'nan'
	assembly_link = 'nan'
	protein_accession = 'nan'
	identical_proteins = []

	if '.' not in ncbi_code:
		ncbi_code = map_uniprot_to_ncbi(ncbi_code)

	try:
		handle = Entrez.efetch(db="protein", id=ncbi_code, rettype="ipg", retmode="xml")
		record = Entrez.read(handle)

		for report in record:
			if 'ProteinList' in record[report]:
				products = record[report]['ProteinList']
				for product in products:
					if 'CDSList' in product:
						if "'assembly':" in str(product):
							cds = product['CDSList']
							cds = str(cds)
							cds = cds.replace('[','').replace(']','').replace('{', '').replace('}', '').replace('(','').replace(')','')

							assemblyId = cds.split("'assembly':")
							
							assemblyId = assemblyId[-1].replace(' ','').replace("'",'').strip()
							assemblyId = assemblyId.replace(',','')

							protein_acc = str(product).split("}, attributes={'accver': ")[1].split(',')[0].replace("'","")
							source = str(product).split("}, attributes={'accver': ")[1].split(", 'source': ")[1].split(',')[0].replace("'","")

							if protein_acc == ncbi_code:
								assembly_id = assemblyId
								assembly_source = source
								assembly_link = database_assembly_mapping[assembly_id]['link']
							else:
								identical_proteins.append([assemblyId, protein_acc])

						elif assembly_id in ['nan', 'entry_removed']:
							assembly_id = 'nan'
							assembly_link = 'nan'

			elif protein_accession == 'nan':
				assembly_id = 'entry_removed'
				assembly_link = 'entry_removed'
	except:
		print(sys.exc_info()[0])
		assembly_id = 'nan'
		assembly_source = 'nan'
		assembly_link = 'nan'

	return ncbi_code, identical_proteins, assembly_id, assembly_link

def collect_assemblies(target_codes, refseq_gb_assembly_map, out_json):

	targets = list(target_codes.keys())

	if os.path.isfile(out_json):
		with open(out_json, 'r') as fp:	 # allows to restart the searches without having to get all the info again in case it crashes
			assemblies = json.load(fp)
			already_collected = list(set([ncbi_code for assembly in assemblies for ncbi_code in assemblies[assembly]['Members']]))
			targets = [target for target in targets if target not in already_collected]

			print(' ... A partial search was already performed for {} entrezIDs. Will continue it for the remaining {}.\n'.format(len(already_collected), len(targets)))
	else:
		assemblies = {}
		already_collected = []

	curr_target_count = 0

	for curr_target_code in targets:
		curr_target_count += 1

		if curr_target_code not in already_collected:

			ncbi_code, identical_proteins, assembly_id, assembly_link = find_ncbi_code_assembly(curr_target_code, refseq_gb_assembly_map)

			if assembly_id not in assemblies:
				assemblies[assembly_id] = {'Members':[], 'identical_proteins': []}

			assemblies[assembly_id]['Members'].append(curr_target_code)
			assemblies[assembly_id]['identical_proteins'].append({curr_target_code: identical_proteins})

			if assembly_id != 'nan' and assembly_id != 'entry_removed':
				if 'Superkingdom' in target_codes[curr_target_code]:
					if 'Superkingdom' not in assemblies[assembly_id]:
						for key in target_codes[curr_target_code]:
							assemblies[assembly_id][key] = target_codes[curr_target_code][key]

					print('({}/{})'.format(curr_target_count, len(targets)), curr_target_code, assembly_id, target_codes[curr_target_code]['Superkingdom'], target_codes[curr_target_code]['Species'], identical_proteins)

			elif assembly_id == 'nan':
				print("\n ... > There is no assembly for {}. There are {} identical_proteins with assembly. ({}/{})\n".format(curr_target_code, len(identical_proteins), curr_target_count, len(targets)))
			else:
				print("\n ... > There is no assembly for {}. Entry removed! ({}/{})\n".format(curr_target_code, curr_target_count, len(targets)))
 
			with open(out_json, 'w') as fp:
			   json.dump(assemblies, fp, indent=4)

	# for assembly_id in assemblies:
	# 	assemblies[assembly_id]['Members'] = list(set(assemblies[assembly_id]['Members']))
	# with open(out_json, 'w') as fp:
	# 	json.dump(assemblies, fp, indent=4)

	return assemblies

def reorganise_assemblies_dict(assemblies, targets, refseq_gb_assembly_map, out_json):

	# will reorganise the assemblies by including the identical ones as individual entries and adding a column with the targets that brought us there

	reorganised_assemblies = {}

	for i, assembly_id in enumerate(assemblies.keys()):
		#print('({}/{})'.format(i+1, len(assemblies)), assembly_id)

		if 'String' in assembly_id:
			assembly_id = assembly_id.split('String')[0]

		if len(assembly_id) != 15:
			print(assembly_id)

		if assembly_id not in reorganised_assemblies:
			reorganised_assemblies[assembly_id] = {key: assemblies[assembly_id][key] for key in assemblies[assembly_id] if key != 'identical_proteins'}
			reorganised_assemblies[assembly_id]['source_targets'] = assemblies[assembly_id]['Members']
			if assembly_id in refseq_gb_assembly_map:
				reorganised_assemblies[assembly_id]['refseq'] = refseq_gb_assembly_map[assembly_id]['refseq']
			else:
				reorganised_assemblies[assembly_id]['refseq'] = 'na'

		for item in assemblies[assembly_id]['identical_proteins']:
			for member in item:
 				identics_to_member = item[member]

 				if len(identics_to_member) > 0:
 					for identic in identics_to_member:
 						identical_protein = identic[1]
 						identical_genome = identic[0]
 						if 'String' in identical_genome:
 							identical_genome = identical_genome.split('String')[0]


 						if identical_genome not in reorganised_assemblies:
 							reorganised_assemblies[identical_genome] = {key: targets[member][key] for key in targets[member]}
 							reorganised_assemblies[identical_genome]['source_targets'] = []
 							reorganised_assemblies[identical_genome]['Members'] = []
 							if identical_genome in refseq_gb_assembly_map:
 								reorganised_assemblies[identical_genome]['refseq'] = refseq_gb_assembly_map[identical_genome]['refseq']
 							else:
 								reorganised_assemblies[identical_genome]['refseq'] = 'na'

 						reorganised_assemblies[identical_genome]['Members'].append(identical_protein)
 						reorganised_assemblies[identical_genome]['source_targets'].append(member)

 						if member == 'OWK39889.1':
 							print(identical_genome, reorganised_assemblies[identical_genome])

	with open(out_json, 'w') as fp:
		json.dump(reorganised_assemblies, fp, indent=4)

	return reorganised_assemblies

def build_summary_table(assemblies, out_csv, n = 10):

	summary_df = {'Assemblies': [],
	              'NoMembers': []}

	for assembly_id in assemblies:
		if assembly_id not in ['nan', 'entry_removed']:
			summary_df['Assemblies'].append(assembly_id)
			summary_df['NoMembers'].append(len(assemblies[assembly_id]['Members']))
			for key in assemblies[assembly_id].keys():
				if key not in summary_df:
					summary_df[key] = []
				summary_df[key].append(assemblies[assembly_id][key])

	summary_df = pd.DataFrame(summary_df)

	if 'Superkingdom' in summary_df.columns:
		phylum_labels = ['{}_{}'.format(row.Superkingdom, row.Phylum) for index, row in summary_df.iterrows()]
		for i, label in enumerate(phylum_labels):
			if 'candidat' in label.lower():
				label = '_'.join([label.split('_')[0], 'candidate_phyla'])
				phylum_labels[i] = label

		summary_df['Phylum_label'] = phylum_labels

	print(' ... All combined (n = {})'.format(len(summary_df)))
	summary_df = summary_df.sort_values(by='NoMembers', ascending = False)
	summary_df.to_csv(out_csv, index=False, sep='\t')
	print(summary_df.head(n))

	# now remove the duplicated GCAs that have a GCF assigned
	summary_df_na = summary_df.loc[(summary_df.refseq == 'na')]
	summary_df_gcfs = summary_df.loc[(summary_df.refseq == summary_df.Assemblies)]
	summary_df = pd.concat([summary_df_na, summary_df_gcfs])

	print(' ... Filtering out repeated ones (deleting those GCA assemblies for which there is an identical entry on RefSeq) (n = {})'.format(len(summary_df)))
	summary_df = summary_df.sort_values(by='NoMembers', ascending = False)
	summary_df.to_csv('{}_filtered.csv'.format(out_csv), index=False, sep='\t')
	print(summary_df.head(n))

	return summary_df

def plot_hist(summary_df, parameter):

	values = summary_df[parameter]
	bins = max(values) - min(values)

	plt.clf()
	plt.hist(values, bins=bins)
	plt.yscale('log')
	plt.savefig('members_per_assenbly_hist.png', format = 'png')

def plot_boxplots(summary_df, parameter, cmap = 'terrain'):

	summary_df = summary_df.sort_values(by='Phylum_label')

	plt.clf()
	plt.figure(figsize=(20,7))
	sns.boxplot(x="Phylum_label", y=parameter, data=summary_df, palette = cmap)
	plt.xticks(rotation='vertical')
	plt.ylabel('Number of targets per genome assembly')
	plt.tight_layout()

	plt.savefig('members_per_assenbly_boxplot.png', format = 'png')

	return summary_df

def plot_completness_bar_plot(assemblies):

	counts = {}

	for assembly_id in assemblies:
		if assembly_id == 'entry_removed':
			label = 'Removed from\ndatabase'
		elif assembly_id == 'nan':
			label = 'Exists but\nno assembly'
		else:
			label = 'Has genome'

		if label not in counts:
			counts[label] = 0
		counts[label] += len(assemblies[assembly_id]['Members'])

	total = sum([counts[label] for label in counts])
	for label in counts:
		counts[label] = counts[label]*100/total
	
	plt.clf()
	plt.bar(counts.keys(), counts.values())
	plt.ylim(0, 100)
	plt.ylabel("% of the input target entrezIDs")
	plt.savefig('completeness_of_assembly_assignments.png', format = 'png')

# MAIN CODE
logfile = 'genome_assembly_collection.log'
sys.stdout = Logger(logfile)

# Download and parse RefSeq and Genbank databases
print('\nDownloading and parsing RefSeq and Genbank summary tables')
refseq_gb_assembly_map = download_and_parse_refseq_and_gb_databases()

# Get the list of target codes
print('\nParsing targets')
targets = parse_targets(args.in_file)

# For each target code, find the assembly to which it belongs to
print('\nCollecting assembly codes for each target')

out_label = '.'.join(args.in_file.split('/')[-1].split('.')[:-1])

out_json = '{}_assemblies.json'.format(out_label)
assemblies = collect_assemblies(targets, refseq_gb_assembly_map, out_json = out_json)

# Build summary table, sorted from the assembly with more of the targets to the one with less
print('\nReorganizing data\n')
out_json = '{}_reorganised_assemblies.json'.format(out_label)
reorganised_assemblies = reorganise_assemblies_dict(assemblies, targets, refseq_gb_assembly_map, out_json = out_json)

print('\nBuilding summary table\n')
out_csv = '{}_assemblies.csv'.format(out_label)
summary_df = build_summary_table(reorganised_assemblies, out_csv = out_csv, n = 10)

# Plot histogram of number of targets per assembly
plot_hist(summary_df, 'NoMembers')
plot_completness_bar_plot(reorganised_assemblies)
if 'Superkingdom' in summary_df.columns:
	plot_boxplots(summary_df, 'NoMembers', cmap = 'terrain')

