import pandas as pd
import matplotlib.cm
import matplotlib
import json
import numpy
import os
import sys 

import matplotlib.pyplot as plt

##clans_file = sys.argv[1]
##ecod_data_json = sys.argv[2]

# for propellers
#clans_file = 'selected_full_sequences_from_ecod.latest.F70_f_name_t_name_Coverage_Number_blades_Length.clans'
#ecod_data_json = 'ecod.latest.F70_selected_propellers.json'

# for pinwheels and prisms
clans_file = 'selected_full_sequences_from_ecod.latest.F70_f_name_t_name_Coverage_Number_blades_Length.clans'
ecod_data_json = 'ecod.latest.F70_selected_propellers_manually_curated.json'

arch_data_dictionary = json.load(open(ecod_data_json,'r'))

## ROUTINES ##

def get_ncbicodes_order_in_clans(clans_file, sep = '|', position = 0):

   print("Collecting ncbiID and their corresponding orders from clans map")
   
   ncbids_ordered = []
   sequences = {}

   count = 0
   found_seq = False
   got_seqid = False
   with open(clans_file, 'r') as clans:
       for line in clans:
           if '<seq>' in line:
               found_seq = True
           elif '</seq>' in line:
               found_seq = False
           elif found_seq and line.startswith('>'):
               line = line[1:]
               line_data = line.split(sep)
               ncbi_code = line_data[position]
               ncbi_code = ncbi_code.replace('"','')
               
               ncbids_ordered.append(ncbi_code)

               got_seqid = True
               
           elif found_seq and not line.startswith('>') and got_seqid:
               sequence = line.strip()
               sequences[ncbi_code] = sequence

   return ncbids_ordered, sequences

def get_clusters_from_clans(clans_file, cluster_codes = ['cluster']):

    print("Collecting clusters from clans map which name start with '{}'".format(cluster_codes))

    ncbis_ordered, sequences = get_ncbicodes_order_in_clans(clans_file, sep = '|', position = 0)
    
    clusters = {}

    found_seqgroup = False
    with open(clans_file, 'r') as in_clans:
        for line in in_clans:

            for cluster_code in cluster_codes:
            
                if '<seqgroup' in line:
                    found_seqgroup = True
                    found_allowed_cluster = False

                elif found_seqgroup:
                    if 'name={}'.format(cluster_code) in line:
                        current_cluster = line.split('=')[-1].strip()
                        found_allowed_cluster = True
                       
                    elif 'numbers=' in line and found_allowed_cluster:
                        numbers = line.split('=')[-1].split(';')[:-1]
                        numbers = [int(i) for i in numbers]
                        numbers = [ncbis_ordered[i] for i in numbers]
                        clusters[current_cluster] = numbers

                        found_allowed_cluster = False

    return clusters, sequences, ncbis_ordered

def fix_and_save_blades_fasta(blades_fasta, out_dir):

    out_fasta = '{}/{}'.format(out_dir, blades_fasta.split('/')[-1])

    with open(out_fasta, 'w') as outfasta:
        with open(blades_fasta, 'r') as infasta:
            for line in infasta:
                if line.startswith('>'):
                    outfasta.write(line)
                elif '//' not in line:
                    outfasta.write(line.replace('-','').upper())
                    
    return out_fasta

def create_input_fasta(fasta_files, out_dir, cluster_label, delete_fastas = True):

    out_fasta = '{}/{}_representative_blades.fasta'.format(out_dir, cluster_label)
    written_sequences = []

    with open(out_fasta, 'w') as outfasta:
       for fasta_file in fasta_files:
          label = fasta_file.split('/')[-1].split('.')[0]
          with open(fasta_file, 'r') as infasta:
             for line in infasta:
                if line.startswith('>'):
                   seq_id = line.strip()
                else:
                   seq = line.strip()
                   if seq not in written_sequences:
                      outfasta.write('{}_{}\n'.format(seq_id, label))
                      outfasta.write(line)
                      written_sequences.append(seq)
          if delete_fastas:
             os.system('rm {}'.format(fasta_file))

    return out_fasta   

   
def select_clusters_representatives(clusters, arch_data_dictionary, n = 1):

    print('Selecting representative propellers for each cluster')
    
    os.system('mkdir {}/selected_cluster_representatives'.format(os.getcwd()))

    selected_representatives = []
    
    for cluster in clusters.keys():
        cluster_data = {'number_blades': [], 'domain_coverages': [], 'family_names': [], 'ecod_id': []}

        for ncbi_code in clusters[cluster]:
            cluster_data['number_blades'].append(len(arch_data_dictionary[ncbi_code]['blades']))
            cluster_data['domain_coverages'].append(arch_data_dictionary[ncbi_code]['coverage'])
            cluster_data['family_names'].append(arch_data_dictionary[ncbi_code]['f'])
            cluster_data['ecod_id'].append(ncbi_code)

        df = pd.DataFrame(cluster_data)
        df = df.sort_values(by = ['number_blades', 'domain_coverages'], ascending = False)

        representatives = df.head(n)
        row_count = 0
        out_dir = ''
        fasta_files = []

        for i, row in representatives.iterrows():             

           selected_ecodid = row['ecod_id']

           representative_label = '{}_{}_{}repeats'.format(cluster.split(':')[-1], row['family_names'], row['number_blades'])
           selected_representatives.append(selected_ecodid)

           if row_count == 0:
              out_dir = '{}/selected_cluster_representatives/{}'.format(os.getcwd(), representative_label)
              os.system('mkdir {}'.format(out_dir))

           os.system('cp {}/selected_pdbs/{}* {}'.format(os.getcwd(), selected_ecodid, out_dir))
           blades_fasta = '{}/selected_pdbs/{}_blades.fasta'.format(os.getcwd(), selected_ecodid)
           blades_fasta = fix_and_save_blades_fasta(blades_fasta, out_dir = out_dir)
           fasta_files.append(blades_fasta)

           print('\nCluster {} ({} entries):'.format(cluster, len(df)))
           print(' ...   Ecod_id: {}'.format(selected_ecodid))
           print(' ...    Family: {}'.format(row['family_names']))
           print(' ...  Topology: {}'.format(arch_data_dictionary[ncbi_code]['t']))
           print(' ... N. blades: {}'.format(row['number_blades']))

           row_count += 1

        create_input_fasta(fasta_files, out_dir, cluster_label = cluster.split(':')[-1], delete_fastas = True)
        
    return selected_representatives

def add_representatives_group_to_clans(clans_file, representatives, size = 9, group_type = 0):

   out_clans_file = '{}_with_ClusterRrepresentatives.clans'.format(clans_file[:-6])
   print("\nWriting cluster representatives to clans file '{}'".format(out_clans_file))
   out_clans = open(out_clans_file, 'w')

   found_seqgroup = False
   with open(clans_file, 'r') as in_clans:
      for line in in_clans:
         out_clans.write(line)         
         if '<seqgroup' in line:
            found_seqgroup = True
            color = "255;255;255;255" # white
            numbers = ";".join([str(numb) for numb in representatives])
            numbers = numbers + ";"
           
            out_clans.write('name=A_representatives\n')
            out_clans.write('type={}\n'.format(group_type))
            out_clans.write('size={}\n'.format(size))
            out_clans.write('color={}\n'.format(color))
            out_clans.write('numbers={}\n'.format(numbers))

      if not found_seqgroup:
         print("There is no '<seqgroup>' section in the clans file")
         print(" ... Set a group first and try again")
         print(" ... Tip: set the 'All' group")

         os.system('rm {}'.format(out_clans_file))

         
# MAIN CODE #

clusters, sequences, ncbis_ordered = get_clusters_from_clans(clans_file, cluster_codes = ['clans'])
representatives = select_clusters_representatives(clusters, arch_data_dictionary, n = 3)

representatives = [ncbis_ordered.index(representative) for representative in representatives]
add_representatives_group_to_clans(clans_file, representatives)
