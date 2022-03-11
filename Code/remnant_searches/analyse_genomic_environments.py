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

# Get inputs
parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments with defaults')

# required inputs
requiredNamed.add_argument('-i', dest='assemblies_csv', type=str, required=True, help='csv file with the assemblies to which targets belong to')
# optional inputs
optionalNamed.add_argument('-n', dest='n_cases', default = -1, type=int, help='Number of cases to take (default: 1)')
optionalNamed.add_argument('-n_max_plots', dest='n_max_plots', default = 10, type=int, help='Maximum number of cases for which plots are to be made (default: 10)')
optionalNamed.add_argument('-target_group', dest='target_group', default = None, type=str, help='Target taxonomic group to focus on (default: None)')
optionalNamed.add_argument('-ecod_domains', dest='ecod_domains', default = None, type=str, help='json file with the domains annotated (default: None)')
optionalNamed.add_argument('-ecod_level', dest='ecod_level', default = 'f_name', type=str, help='ECOD level to map (default: f_name)')
optionalNamed.add_argument('-ecod_cmap', dest='ecod_cmap', default = 'Spectral', type=str, help='colormap to color the ECOD domains (default: Spectral)')
optionalNamed.add_argument('-min_pop', dest='min_pop', default = 2, type=int, help='Minimum number of targets per assembly for it to be considered. Only used for the plotting of propeller family frequencies at the end (default: 2)')
optionalNamed.add_argument('-propellers', dest='propellers', default = None, type=str, help='json file with propellers all data (default: None)')
optionalNamed.add_argument('-graveyards', dest='graveyards', default = None, type=str, help='json file with graveyards in intergenomic regions (default: None)')

optionalNamed.add_argument('-tmp_folder', dest='tmp_folder', default = '/tmp/',type=str, help='Temporary folder (default: /tmp/)')

# Define inputs
args = parser.parse_args()

# HELPING ROUTINES

# parsers

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
				curr_scaffold = '{} : {}-{}'.format(*line.split()[1:])

	return assembly	


# others

def select_target_assemblies(assemblies_df, n, focus_on = args.target_group):

	if focus_on != None:
		assemblies_df = assemblies_df[assemblies_df.apply(lambda r: r.str.contains(focus_on).any(), axis=1)]
	
	if n != -1:
		assemblies_df = assemblies_df.head(n)

	assemblies_df = assemblies_df.reset_index(drop=True)

	return assemblies_df

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

def find_scaffolds_with_targets(assembly, targets, assembly_link, tmp_folder = args.tmp_folder, label = ''):

	relevant_scaffolds = {}

	for i, ncbi_code in enumerate(assembly['ncbi_codes']):
		for target in targets:
			if ncbi_code == target:
				scaffold = assembly['scaffolds'][i]
				start = assembly['starts'][i]
				end = assembly['ends'][i]
				direction = assembly['directions'][i]

				if scaffold not in relevant_scaffolds:
					relevant_scaffolds[scaffold] = {'ncbi_codes': [],
					                                'starts': [],
					                                'ends': [],
					                                'directions':[],
					                                'length': int(scaffold.split()[-1].split('-')[-1]) - int(scaffold.split()[-1].split('-')[0])+1}

				relevant_scaffolds[scaffold]['ncbi_codes'].append(ncbi_code)
				relevant_scaffolds[scaffold]['starts'].append(start)
				relevant_scaffolds[scaffold]['ends'].append(end)
				relevant_scaffolds[scaffold]['directions'].append(direction)

	print(' ... ... Found {} scaffolds with targets'.format(len(relevant_scaffolds)))

	return relevant_scaffolds

def get_distances_of_targets_to_adjacent_neighbous(relevant_scaffolds, assembly, distances_to_neighbours = {}):

	if 'to_any_neighbour' not in distances_to_neighbours:
		distances_to_neighbours['to_any_neighbour'] = {'to_previous': [], 'to_next': []}
	if 'to_target_neigbours' not in distances_to_neighbours:
		distances_to_neighbours['to_target_neigbours'] = {'to_previous': [], 'to_next': []}

	for scaffold in relevant_scaffolds:
		for i, target in enumerate(relevant_scaffolds[scaffold]['ncbi_codes']):
			target_direction = relevant_scaffolds[scaffold]['directions'][i]

			# select only those entries in the assembly in the scaffold AND in the same strand (same direction)
			neighbours_indeces = [j for j, x in enumerate(assembly['scaffolds']) if x == scaffold and assembly['directions'][j] == target_direction]
			target_index = assembly['ncbi_codes'].index(target)
			
			# create a vector with the indeces of the target gene and the ones right flanking the it
			if len(neighbours_indeces) > 1:
				if neighbours_indeces.index(target_index) > 0 and neighbours_indeces.index(target_index) < len(neighbours_indeces)-1:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)-1], neighbours_indeces[neighbours_indeces.index(target_index)], neighbours_indeces[neighbours_indeces.index(target_index)+1]]
				elif neighbours_indeces.index(target_index) == 0:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)], neighbours_indeces[neighbours_indeces.index(target_index)+1]] 
				else:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)-1], neighbours_indeces[neighbours_indeces.index(target_index)]]

				# compute the distance between the 2 ends of 2 flanking genes
				for j, curr_gene in enumerate(genomic_context[:-1]):
					next_gene = genomic_context[j+1]

					curr_code = assembly['ncbi_codes'][curr_gene]
					next_code = assembly['ncbi_codes'][next_gene]

					curr_end = assembly['ends'][curr_gene]
					next_start = assembly['starts'][next_gene]

					distance = next_start-curr_end

					if curr_code == target:
						if target_direction == '+':
							distance_type = 'to_next'
						else:
							distance_type = 'to_previous'

					elif next_code == target:
						if target_direction == '+':
							distance_type = 'to_previous'
						else:
							distance_type = 'to_next'

					distances_to_neighbours['to_any_neighbour'][distance_type].append(distance)

					if curr_code in relevant_scaffolds[scaffold]['ncbi_codes'] and next_code in relevant_scaffolds[scaffold]['ncbi_codes']:
						distances_to_neighbours['to_target_neigbours'][distance_type].append(distance)

	return distances_to_neighbours

def get_distances_of_targets_to_closest_targets(relevant_scaffolds, assembly, distances_to_neighbours = {}):

	if 'to_close_targets' not in distances_to_neighbours:
		distances_to_neighbours['to_close_targets'] = {'to_previous': [], 'to_next': []}

	for scaffold in relevant_scaffolds:
		for i, target in enumerate(relevant_scaffolds[scaffold]['ncbi_codes']):
			target_direction = relevant_scaffolds[scaffold]['directions'][i]

			# select only those targets in the scaffold that are in the same strand (same direction)
			neighbours_indeces = [j for j, direction in enumerate(relevant_scaffolds[scaffold]['directions']) if direction == target_direction]
			target_index = relevant_scaffolds[scaffold]['ncbi_codes'].index(target)

			if len(neighbours_indeces) > 1:
				if neighbours_indeces.index(target_index) > 0 and neighbours_indeces.index(target_index) < len(neighbours_indeces)-1:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)-1], neighbours_indeces[neighbours_indeces.index(target_index)], neighbours_indeces[neighbours_indeces.index(target_index)+1]]
				elif neighbours_indeces.index(target_index) == 0:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)], neighbours_indeces[neighbours_indeces.index(target_index)+1]] 
				else:
					genomic_context = [neighbours_indeces[neighbours_indeces.index(target_index)-1], neighbours_indeces[neighbours_indeces.index(target_index)]]

				# compute the distance between the 2 ends of 2 flanking genes
				for j, curr_gene in enumerate(genomic_context[:-1]):
					next_gene = genomic_context[j+1]

					curr_code = relevant_scaffolds[scaffold]['ncbi_codes'][curr_gene]
					next_code = relevant_scaffolds[scaffold]['ncbi_codes'][next_gene]

					curr_end = relevant_scaffolds[scaffold]['ends'][curr_gene]
					next_start = relevant_scaffolds[scaffold]['starts'][next_gene]

					distance = next_start-curr_end

					if curr_code == target:
						if target_direction == '+':
							distance_type = 'to_next'
						else:
							distance_type = 'to_previous'

					elif next_code == target:
						if target_direction == '+':
							distance_type = 'to_previous'
						else:
							distance_type = 'to_next'

					distances_to_neighbours['to_close_targets'][distance_type].append(distance)


	return distances_to_neighbours

def update_graveyards_with_distances_to_targets(graveyards, assembly_id, relevant_scaffolds):

	for scaffold in relevant_scaffolds.keys():
		scaffold_id = scaffold.split(' : ')[0]

		if scaffold_id in graveyards[assembly_id]:
			curr_graveyards = graveyards[assembly_id][scaffold_id]

			graveyards[assembly_id][scaffold_id]['Closest_target'] = []
			graveyards[assembly_id][scaffold_id]['Closest_distance'] = []
			graveyards[assembly_id][scaffold_id]['Closest_position'] = []

			for i, graveyard in enumerate(curr_graveyards['graveyard_names']):
				grave_direction = curr_graveyards['directions'][i]
				grave_start = curr_graveyards['starts'][i]
				grave_end = curr_graveyards['ends'][i]

				closest_target = None
				closest_distance = None
				closest_position = None

				for j, target in enumerate(relevant_scaffolds[scaffold]['ncbi_codes']):
					target_direction = relevant_scaffolds[scaffold]['directions'][j]

					if target_direction == grave_direction:
						target_start = relevant_scaffolds[scaffold]['starts'][j]
						target_end = relevant_scaffolds[scaffold]['ends'][j]

						if grave_end < target_start:
							position = "5p"
							distance = target_start - grave_end
						elif target_end < grave_start:
							position = "3p"
							distance = grave_start - target_end
						else:
							position = 'overlap'
							distance = None

						if distance is not None and (closest_target is None or distance < closest_distance):
							closest_target = target
							closest_distance = distance
							closest_position = position

				graveyards[assembly_id][scaffold_id]['Closest_target'].append(closest_target)
				graveyards[assembly_id][scaffold_id]['Closest_distance'].append(closest_distance)
				graveyards[assembly_id][scaffold_id]['Closest_position'].append(closest_position)

	return graveyards

def add_propellers_domain_mappings_to_scaffolds(relevant_scaffolds, ecod_domains, ecod_level):

	ecod_domains = json.load(open(ecod_domains, 'r'))
	all_labels = {}

	for scaffold in relevant_scaffolds:
		relevant_scaffolds[scaffold][ecod_level] = []

		for target in relevant_scaffolds[scaffold]['ncbi_codes']:
			for entry in ecod_domains.keys():
				annotated_domains = ecod_domains[entry]['annotated_domains']
				found_propeller = False
				for domain in annotated_domains:
					if 'beta-propeller-like' in domain:
						if ecod_level == 'f_name':
							label = domain.split('F: ')[-1]
							if len(label.split('_')) > 1:
								label = '_'.join(label.split('_')[:2])
								if label.split('_')[-1].isdigit():
									label = label.split('_')[0]
						elif ecod_level == 't_name':
							label = domain.split(' F: ')[0].split('T: ')[-1]
							if len(label.split(' ')) > 4:
								if 'domain' in label.split(' ')[-1]:
									label = ' '.join(label.split(' ')[:-1])
								else:
									label = ' '.join(label.split(' ')[-4:])

						if ecod_level not in all_labels:
							all_labels[ecod_level] = []
						all_labels[ecod_level].append(label)

						#all_labels.append(label)
						found_propeller = True

				if target in entry:
					if found_propeller:
						relevant_scaffolds[scaffold][ecod_level].append(label)
					else:
						relevant_scaffolds[scaffold][ecod_level].append('UNKNOWN')

	return relevant_scaffolds, all_labels

def add_propellers_properties(relevant_scaffolds, propellers):

	print(' ... ... Adding properties to relevant scaffolds')

	propellers = json.load(open(propellers, 'r'))
	all_ecod_labels = {}

	# print(relevant_scaffolds)
	# print(propellers)

	for scaffold in relevant_scaffolds:
		for target in relevant_scaffolds[scaffold]['ncbi_codes']:
			propellers_in_target = []
			for propeller in propellers.keys():
				if target in propeller:
					propellers_in_target.append(propellers[propeller])

				for key in propellers[propeller]:
					if '_name' in key or 'Classification' in key:
						if key not in all_ecod_labels:
							all_ecod_labels[key] = []

						all_ecod_labels[key].append(propellers[propeller][key])

			selected_propeller = propellers_in_target[0]
			for propeller in propellers_in_target:
				if selected_propeller['Classification'] == 'Incomplete' and propeller['Classification'] == 'Complete':
					selected_propeller = propeller
				elif selected_propeller['Classification'] == 'Complete' and selected_propeller['blades_similarity'] < propeller['blades_similarity']:
					selected_propeller = propeller

			for key in selected_propeller.keys():
				if key not in ['Superkingdom', 'Phylum', 'Class', 'Order', 'Genus', 'Species']:
					if key not in relevant_scaffolds[scaffold]:
						relevant_scaffolds[scaffold][key] = []
					relevant_scaffolds[scaffold][key].append(selected_propeller[key])


	for key in all_ecod_labels:
		all_ecod_labels[key] = sorted(list(set(all_ecod_labels[key])))

	return relevant_scaffolds, all_ecod_labels


# for plotting

def get_domains_colors(relevant_scaffolds, all_labels, domains_level, cmap):

	colors = {}
	legend_handles = {}

	if domains_level == None or domains_level not in relevant_scaffolds[list(relevant_scaffolds.keys())[0]]:
		for scaffold in relevant_scaffolds:
			for i, target in enumerate(relevant_scaffolds[scaffold]['ncbi_codes']):
				colors[target] = {'fillcolor': 'grey',
				                  'edgecolor': 'black'}
	else:
		#print(all_labels, domains_level)
		cmap = matplotlib.cm.get_cmap(cmap)
		norm = matplotlib.colors.Normalize(vmin=0, vmax=len(all_labels[domains_level]))

		colours = [cmap(norm(i)) for i in range(len(all_labels[domains_level]))]

		for scaffold in relevant_scaffolds:
			for i, target in enumerate(relevant_scaffolds[scaffold]['ncbi_codes']):
				domain_label = relevant_scaffolds[scaffold][domains_level][i]
				color = colours[all_labels[domains_level].index(domain_label)]

				colors[target] = {'fillcolor': color,
				                  'edgecolor': color}

				handle = mpatches.Patch(color=color, label=domain_label)
				if domain_label not in legend_handles:
					legend_handles[domain_label] = handle

	legend_handles = [legend_handles[domain_label] for domain_label in sorted(legend_handles.keys())]

	return colors, legend_handles

def get_values_colors(relevant_scaffolds, propertie, cmap = 'Reds'):

	colors = {}
	all_values = [0, 1]
	if propertie == 'graveyards':
		all_values = [1 for scaffold in relevant_scaffolds for value in relevant_scaffolds[scaffold]['ncbi_codes']]
		cmap = 'binary'
	else:
		all_values = sorted([value for scaffold in relevant_scaffolds for value in relevant_scaffolds[scaffold][propertie]])
	
	cmap = matplotlib.cm.get_cmap(cmap)
	norm = matplotlib.colors.Normalize(vmin=min(all_values), vmax=max(all_values))

	colours = [cmap(norm(i)) for i in all_values]

	for scaffold in relevant_scaffolds:
		for i, target in enumerate(relevant_scaffolds[scaffold]['ncbi_codes']):
			if propertie == 'graveyards':
				color = 'black'
			else:
				curr_value = relevant_scaffolds[scaffold][propertie][i]
				color = colours[all_values.index(curr_value)]

			colors[target] = {'fillcolor': color,
			                  'edgecolor': color}

	colorbar_mappeable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

	return colors, colorbar_mappeable

def get_graveyards_colors(scaffold_graveyards, cmap = 'Reds'):

	colors = {}

	cmap = matplotlib.cm.get_cmap(cmap)
	norm = matplotlib.colors.Normalize(vmin=0, vmax=100)

	all_values = sorted([value for value in scaffold_graveyards['probabilities']])
	colours = [cmap(norm(i)) for i in all_values]

	for i, target in enumerate(scaffold_graveyards['orf_names']):

		curr_value = scaffold_graveyards['probabilities'][i]
		color = colours[all_values.index(curr_value)]

		colors[target] = {'fillcolor': color,
		                  'edgecolor': color}

	colorbar_mappeable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

	return colors, colorbar_mappeable



def plot_genomic_distributions(assembly, relevant_scaffolds, all_ecod_labels, graveyards = {}, ecod_cmap = 'Spectral', assembly_id = '', title = ''):

	out_pdf = '{}_targets_genomic_landscape.pdf'.format(assembly_id)
	pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)

	all_properties = set([key for scaffold in relevant_scaffolds for key in relevant_scaffolds[scaffold] if key not in ['Superkingdom', 'Phylum', 'Class', 'Order', 'Genus', 'Species', 'ncbi_codes', 'starts', 'ends', 'directions','length']])

	if len(all_properties) == 0:
		all_properties = ['graveyards']
	else:
		all_properties = sorted(list(all_properties))
		all_properties.append('graveyards')

	print(' ... ... Making genomic landscape picture')

	if assembly_id in graveyards:
		assembly_graveyards = graveyards[assembly_id]
	else:
		assembly_graveyards = {}

	for count, propertie in enumerate(all_properties):
		added_graveyards = False
		print( ' ... ... ... Property: {} ({}/{})'.format(propertie, count+1, len(all_properties)))

		legend_handles = None
		if '_name' in propertie or 'Classification' in propertie:
			colors, legend_handles = get_domains_colors(relevant_scaffolds, all_ecod_labels, propertie, ecod_cmap)
		else:
			colors, colorbar_mappeable = get_values_colors(relevant_scaffolds, propertie)

		curr_y_level = len(relevant_scaffolds.keys())
		all_xs = []
		yticklabels = []
		targets_count = 0

		plt.clf()

		if len(relevant_scaffolds) == 1:
			fig, ax = plt.subplots(1, 1, figsize=(20, 2))
		elif len(relevant_scaffolds) < 5:
			fig, ax = plt.subplots(1, 1, figsize=(20, len(relevant_scaffolds)))
		else:
			fig, ax = plt.subplots(1, 1, figsize=(20, int(len(relevant_scaffolds)/1.5)))

		for scaffold in sorted(list(relevant_scaffolds.keys())):
			length = relevant_scaffolds[scaffold]['length']
			all_xs.append(length)
			scaffold_id = scaffold.split()[0]

			ax.hlines(y=curr_y_level, xmin = 0, xmax=length, color='black', linestyle='-')

			# add first all genes in the scaffold as lightgrey arrows
			for i, curr_scaffold in enumerate(assembly['scaffolds']):
				if scaffold == curr_scaffold:
					gene_code = assembly['ncbi_codes'][i]
					if gene_code not in relevant_scaffolds[scaffold]['ncbi_codes']:
						gene_dx = assembly['ends'][i] - assembly['starts'][i]
						gene_direction = assembly['directions'][i]

						if gene_direction == '-':
							gene_x_tail = assembly['ends'][i]
							gene_dx = gene_dx*(-1)
							gene_y = curr_y_level-0.1
						else:
							gene_x_tail = assembly['starts'][i]
							gene_y = curr_y_level+0.1

						ax.arrow(gene_x_tail, gene_y, gene_dx, 0, width=0.125, head_width=0.125, length_includes_head = True, head_length = abs(gene_dx*0.5), facecolor = 'gainsboro', edgecolor = 'gainsboro')


			# now add the targets properly colored based on the property
			for i, target in enumerate(relevant_scaffolds[scaffold]['ncbi_codes']):
				gene_dx = relevant_scaffolds[scaffold]['ends'][i] - relevant_scaffolds[scaffold]['starts'][i]
				gene_direction = relevant_scaffolds[scaffold]['directions'][i]

				if gene_direction == '-':
					gene_x_tail = relevant_scaffolds[scaffold]['ends'][i]
					gene_dx = gene_dx*(-1)
					gene_y = curr_y_level-0.1
				else:
					gene_x_tail = relevant_scaffolds[scaffold]['starts'][i]
					gene_y = curr_y_level+0.1

				ax.arrow(gene_x_tail, gene_y, gene_dx, 0, width=0.125, head_width=0.125, length_includes_head = True, head_length = abs(gene_dx*0.5), facecolor = colors[target]['fillcolor'], edgecolor = colors[target]['edgecolor'])
				targets_count += 1

			# now, if graveyards are given for this assembly, map them
			if scaffold_id in assembly_graveyards and propertie == 'graveyards':
				graveyards_colors, colorbar_mappeable = get_graveyards_colors(assembly_graveyards[scaffold_id])
				added_graveyards = True

				for i, graveyard in enumerate(assembly_graveyards[scaffold_id]['orf_names']):
					graveyard_dx = assembly_graveyards[scaffold_id]['ends'][i] - assembly_graveyards[scaffold_id]['starts'][i]
					graveyard_direction = assembly_graveyards[scaffold_id]['directions'][i]
					graveyard_probability = assembly_graveyards[scaffold_id]['probabilities'][i]

					if graveyard_direction == '-':
						graveyard_x_tail = assembly_graveyards[scaffold_id]['ends'][i]
						graveyard_dx = graveyard_dx*(-1)
						graveyard_y = curr_y_level-0.3
						
					else:
						graveyard_x_tail = assembly_graveyards[scaffold_id]['starts'][i]
						graveyard_y = curr_y_level+0.3

					ax.arrow(graveyard_x_tail, graveyard_y, graveyard_dx, 0, width=0.125, head_width=0.125,length_includes_head = True, head_length = 0, facecolor = graveyards_colors[graveyard]['fillcolor'], edgecolor = graveyards_colors[graveyard]['edgecolor'])

			curr_y_level -= 1

			yticklabels.append(scaffold)

		yticklabels.append('')
		yticklabels.reverse()

		ax.set_yticks(np.arange(0, len(yticklabels)+1, 1.0))
		ax.set_yticklabels(yticklabels)
		ax.spines['right'].set_visible(False)
		ax.spines['left'].set_visible(False)

		ax.set_xlim(0, max(all_xs))
		
		if legend_handles != None:
			ax.legend(handles=legend_handles, loc=0, title = 'ECOD propeller {}'.format(propertie))
		else:
			cbar = plt.colorbar(colorbar_mappeable, ax = ax)
			if added_graveyards:
				propertie = 'graveyard prob.'
			cbar.ax.set_title(propertie)

		plt.title('Scaffolds from assembly {} ({}) with targets: {} targets'.format(assembly_id, title, targets_count))

		#plt.tight_layout()
		pdf.savefig(fig, bbox_inches='tight')

		plt.close(fig)

	pdf.close()

def get_distance_types_colors(distances_collected, cmap = 'RdYlBu'):

	colors = {}
	legend_handles = {}

	cmap = matplotlib.cm.get_cmap(cmap)
	norm = matplotlib.colors.Normalize(vmin=0, vmax=len(distances_collected))

	colours = [cmap(norm(i)) for i in range(len(distances_collected))]

	for i, distance_type in enumerate(sorted(list(distances_collected.keys()))):
		color = colours[i]
		colors[distance_type] = color

		handle = mpatches.Patch(color=color, label=distance_type)
		if distance_type not in legend_handles:
			legend_handles[distance_type] = handle

	legend_handles = [legend_handles[distance_type] for distance_type in sorted(legend_handles.keys())]

	return colors, legend_handles

def plot_distances_histograms(distances_collected, cmap = 'GnBu_r', nbins = 150):

	colors, legend_handles = get_distance_types_colors(distances_collected, cmap)

	plt.clf()

	fig, ax = plt.subplots(1, len(distances_collected), figsize=(len(distances_collected)*5, 5), sharey=True)

	for i, distance_type in enumerate(distances_collected):
		distances = [distance for to_whom in distances_collected[distance_type] for distance in distances_collected[distance_type][to_whom]]
		bins = list(np.linspace(min(distances), max(distances), nbins))

		ax[i].hist(distances, bins = bins, color = colors[distance_type])

		ax[i].set_xlabel('Distance (bp)')
		ax[i].set_ylabel('Counts')

		ax[i].set_yscale('log')
		ax[i].set_xlim(min(distances), max(distances))

		median = np.median(distances)
		ax[i].set_title('Median: {}'.format(round(median)))

	fig.legend(handles=legend_handles, title = 'Distance type')

	fig.suptitle('no. bins: {}'.format(nbins))
	plt.xticks(rotation = 90)
	plt.savefig('distances_to_neighbours_hist.pdf', format='pdf', bbox_inches='tight')
	plt.close('all')

# for the final boxplot

def parse_domains_from_ecod(ecod_domains, ecod_level):

	targets_domains = {}

	for entry in ecod_domains.keys():
		annotated_domains = ecod_domains[entry]['annotated_domains']
		found_propeller = False
		for domain in annotated_domains:
			if 'beta-propeller-like' in domain:
				if ecod_level == 'f_name':
					label = domain.split('F: ')[-1]
					if len(label.split('_')) > 1:
						label = '_'.join(label.split('_')[:2])
						if label.split('_')[-1].isdigit():
							label = label.split('_')[0]
				elif ecod_level == 't_name':
					label = domain.split(' F: ')[0].split('T: ')[-1]
					if len(label.split(' ')) > 4:
						if 'domain' in label.split(' ')[-1]:
							label = ' '.join(label.split(' ')[:-1])
						else:
							label = ' '.join(label.split(' ')[-4:])

				if entry not in targets_domains:
					targets_domains[entry] = []
				targets_domains[entry].append(label)

	return targets_domains

def plot_propeller_families_frequency_boxplot(assemblies_df, ecod_domains, domains_level = 'f_name', cmap = 'Spectral', min_pop = 2):

	data = {domains_level: [], 'frequency': []}
	ecod_domains = json.load(open(ecod_domains,'r'))

	selected_assemblies = assemblies_df.loc[assemblies_df.NoMembers >= min_pop]
	targets_domains = parse_domains_from_ecod(ecod_domains, domains_level)

	for index, assembly_row in selected_assemblies.iterrows():
		domains_found = []
		for target in ast.literal_eval(assembly_row.Members):
			for target_with_domains in targets_domains:
				if target in target_with_domains:
					for domain in set(targets_domains[target_with_domains]):
						domains_found.append(domain)

		unique_domains_found = sorted(list(set(domains_found)))
		counts = [domains_found.count(domain) for domain in unique_domains_found]
		for i in range(len(unique_domains_found)):
			data[domains_level].append(unique_domains_found[i])
			data['frequency'].append(counts[i]*100/sum(counts))

	data = pd.DataFrame(data)
	data = data.sort_values(by=domains_level)

	plt.clf()
	plt.figure(figsize=(20,7))
	sns.boxplot(x=domains_level, y='frequency', data=data, palette = cmap)
	plt.xticks(rotation='vertical')
	plt.title('Frequency of ECOD {} in individual assemblies (for those with at least {} targets; {} out of {})'.format(domains_level, min_pop, len(selected_assemblies), len(assemblies_df)))
	plt.ylabel('Frequency (%)')
	plt.tight_layout()

	plt.savefig('ecod_domains_frequency_per_assenbly_boxplot_min_pop_{}.pdf'.format(min_pop), bbox_inches='tight', format = 'pdf')
	plt.close('all')

def plot_propeller_families_frequency_coccurrence(assemblies_df, ecod_domains, domains_level = 'f_name', cmap = 'Blues', min_pop = 2):

	ecod_domains = json.load(open(ecod_domains,'r'))
	targets_domains = parse_domains_from_ecod(ecod_domains, domains_level)

	unique_domains = sorted(list(set([domain for target in targets_domains for domain in targets_domains[target]])))

	coccurrency_matrix = [[0 for domain in unique_domains] for domain in unique_domains]

	selected_assemblies = assemblies_df.loc[assemblies_df.NoMembers >= min_pop]

	for index, assembly_row in selected_assemblies.iterrows():
		domains_found = []
		for target in ast.literal_eval(assembly_row.Members):
			for target_with_domains in targets_domains:
				if target in target_with_domains:
					for domain in set(targets_domains[target_with_domains]):
						domains_found.append(domain)

		for i, domain in enumerate(domains_found):
			co_occurs_with = [domains_found[j] for j in range(len(domains_found)) if j != i]
			
			index_of_domain = unique_domains.index(domain)
			for co_ccur_domain in co_occurs_with:
				index_of_co_ccur = unique_domains.index(co_ccur_domain)

				coccurrency_matrix[index_of_domain][index_of_co_ccur] += 1

	for i, domain in enumerate(unique_domains):
		sum_occurrences = sum(coccurrency_matrix[i])
		if sum_occurrences > 0:
			for j in range(len(unique_domains)):
				coccurrency_matrix[i][j] = coccurrency_matrix[i][j]*100/sum_occurrences

	coccurrency_matrix = np.array(coccurrency_matrix).T

	plt.clf()

	fig, ax = plt.subplots(figsize=(int(len(unique_domains)/4), int(len(unique_domains)/4)))

	im = ax.imshow(coccurrency_matrix, cmap = cmap, vmin = 0, vmax = 100)

	ax.set_xticks(np.arange(len(unique_domains)))
	ax.set_yticks(np.arange(len(unique_domains)))
	ax.set_yticklabels(unique_domains)
	ax.set_xticklabels(unique_domains, rotation = 90)

	ax.set_ylim(len(coccurrency_matrix)-0.5, -0.5)
	ax.set_xlabel('ECOD {}'.format(domains_level))
	ax.set_ylabel('Co-occurs with')

	fig.colorbar(im,fraction=0.046, pad=0.04)

	plt.savefig('ecod_domains_coccurrence_per_assenbly_heatmap_min_pop_{}.pdf'.format(min_pop), format = 'pdf', bbox_inches='tight')
	plt.close('all')

def plot_graveyard_distribution_statistics(assemblies_df, graveyards):

	out_pdf = 'ALL_graveyards_stats.pdf'
	pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)

	# scatter plot of no. propellers vs no. graveyards

	limits = [0, max(list(assemblies_df.NoMembers) + list(assemblies_df.graveyard_counts))]
	plt.clf()
	fig, ax = plt.subplots(1, 1, figsize=(5, 5))
	ax.scatter(assemblies_df.NoMembers, assemblies_df.graveyard_counts, s = 1, color = 'black')
	ax.plot(limits, limits, lw = 1, linestyle = ':', color = 'black')
	ax.set_ylim(limits[0], limits[1])
	ax.set_xlim(limits[0], limits[1])
	ax.set_ylabel('No. propeller graveyards')
	ax.set_xlabel('No. sequences with highly repetitive propellers')

	plt.title('Number of propeller graveyards per assembly vs.\nthe number of sequences with highly repetitive propellers\n')

	pdf.savefig(fig, bbox_inches='tight')
	plt.close(fig)

	pdf.close()



# MAIN CODE

# Download and parse RefSeq and Genbank databases
print('\nDownloading and parsing RefSeq and Genbank summary tables')
refseq_gb_assembly_map = download_and_parse_refseq_and_gb_databases()

# select the target assemblies 
assemblies_df = pd.read_csv(args.assemblies_csv, sep = '\t')
target_assemblies_df = select_target_assemblies(assemblies_df, n = args.n_cases, focus_on = args.target_group)

if args.graveyards != None and os.path.isfile(args.graveyards):
	graveyards = json.load(open(args.graveyards, 'r'))
	target_assemblies_df = target_assemblies_df.loc[target_assemblies_df.apply(lambda x: x.Assemblies in list(graveyards.keys()), axis=1)]

else:
	graveyards = {}

print(target_assemblies_df)

adjacent_distances_collected = {}

# for target assembly, download its summary and plot the distribution of the targets in it
for index, assembly_row in target_assemblies_df.iterrows():
	assembly_id = assembly_row.Assemblies

	print('\n({}/{}) Assembly: {}\t{}\tSpecies: {}'.format(index+1, len(target_assemblies_df), assembly_id, assembly_row.Phylum_label, assembly_row.Species))
	print(' ... Encompasses {} targets'.format(assembly_row.NoMembers))
	assembly, assembly_link = download_and_extract_assembly(assembly_id, refseq_gb_assembly_map, label = assembly_id)

	if assembly != 'nan':
		all_ecod_labels = []
		relevant_scaffolds = find_scaffolds_with_targets(assembly, ast.literal_eval(assembly_row.Members), assembly_link, label = assembly_id)
		adjacent_distances_collected = get_distances_of_targets_to_adjacent_neighbous(relevant_scaffolds, assembly, adjacent_distances_collected)
		adjacent_distances_collected = get_distances_of_targets_to_closest_targets(relevant_scaffolds, assembly, adjacent_distances_collected)

		if args.ecod_domains != None and os.path.isfile(args.ecod_domains):
			relevant_scaffolds, all_ecod_labels = add_propellers_domain_mappings_to_scaffolds(relevant_scaffolds, args.ecod_domains, args.ecod_level)

		if args.propellers != None and os.path.isfile(args.propellers):
			relevant_scaffolds, all_ecod_labels = add_propellers_properties(relevant_scaffolds, args.propellers)

		# make figures that shows how the targets are spread through the scaffolds found
		if index < args.n_max_plots:
			plot_genomic_distributions(assembly, relevant_scaffolds, all_ecod_labels, graveyards = graveyards, ecod_cmap = args.ecod_cmap, assembly_id = assembly_id, title = '{} {} {}'.format(assembly_row.Superkingdom, assembly_row.Phylum, assembly_row.Species))

		if len(graveyards) > 0:
			graveyards = update_graveyards_with_distances_to_targets(graveyards, assembly_id, relevant_scaffolds) 

			out_graveyards = '{}_with_distances.json'.format(args.graveyards.replace('.json',''))
			json.dump(graveyards, open(out_graveyards,'w'), indent=4)

# plot the histograms of general distances between targets and their adjacent genes (or their closest, to each side, other targets)
plot_distances_histograms(adjacent_distances_collected)

# independently on the target assemblies, plot the frequency of propeller families over individual assemblies
if args.ecod_domains != None and os.path.isfile(args.ecod_domains):

	print('\nPlotting overall (all assemblies in input, no matter the targets) propeller families frequencies over assemblies')
	print(' ... Doing it for the cases with at least {} targets in the assembly'.format(args.min_pop))
	plot_propeller_families_frequency_boxplot(assemblies_df, args.ecod_domains, domains_level = args.ecod_level, cmap = args.ecod_cmap, min_pop = args.min_pop)

	print('\nPlotting overall (all assemblies in input, no matter the targets) propeller families co-occurrence over assemblies')
	print(' ... Doing it for the cases with at least {} targets in the assembly'.format(args.min_pop))
	plot_propeller_families_frequency_coccurrence(assemblies_df, args.ecod_domains, domains_level = args.ecod_level, min_pop = args.min_pop)

if 'graveyard_counts' in assemblies_df.columns:

	print('\nPlotting overall (all assemblies in input, no matter the targets) graveyard distribution statistics')
	plot_graveyard_distribution_statistics(assemblies_df, args.graveyards)


