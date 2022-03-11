import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import sys 
import os
import os.path
import json

tmp_folder = '/tmp/'

# define executables
symmCE = '/Applications/cesymm-2.1.0/runCESymm.sh'

# define routines

def parse_propellers_from_ECOD(ecod_table, exclude_terms = ['beta-pinwheel', 'beta-Prism']):

    propellers = {}

    with open(ecod_table, 'r') as inecod:
        for line in inecod:
            if 'beta-propeller-like' in line and not True in [term in line for term in exclude_terms]:
                propeller_data = line.split('\t')
                ecod_domain_id = propeller_data[1]
                arch_name = propeller_data[9].replace('"', '')
                x_name = propeller_data[10].replace('"', '')
                h_name = propeller_data[11].replace('"', '')
                t_name = propeller_data[12].replace('"', '')
                f_name = propeller_data[13].replace('"', '')

                propellers[ecod_domain_id] = {'arch': arch_name, 'x': x_name, 'h': h_name, 't': t_name, 'f': f_name}

    return propellers

def get_propellers_sequence(propellers, ecod_fasta):

    with open(ecod_fasta,'r') as inecod:
        for line in inecod:
            if line.startswith('>'):
                ecod_code = line.split('|')[1].replace('>', '')
                seq = ''
            else:
                seq += line.strip()
                if ecod_code in propellers.keys():
                    propellers[ecod_code]['seq'] = seq

    return propellers
    
def parse_symmCE_blades(symmCE_outfasta):

    blades = {}

    with open(symmCE_outfasta, 'r') as symmce:
        blade_count = 0
        for line in symmce:
            if line.startswith('>'):
                blade_count += 1
                label = '{}_blade_{}'.format(line.strip()[1:], blade_count)
            elif len(line) > 3:
                #sequence = line.replace('-', '').upper()
                sequence = line.upper()

                blades[label] = sequence
              
    return blades

def parse_symmCE_stats(symmCE_outsimple):

    tmscore = 0
    
    with open(symmCE_outsimple, 'r') as symmce:
        for line in symmce:
            if 'NumRepeats' not in line and len(line) > 4:
                line_data = line.strip().split('\t')
                tmscore = line_data[10]
                tmscore = float(tmscore.split(')')[0].replace(',', '.'))
                if len(line_data) > 16:
                    intervals = line_data[16]
                    if intervals != 'null':
                        intervals = intervals.split(';')
                        intervals = [interval.replace('null.', '') for interval in intervals]
                        intervals = [interval.split('_')[-1] for interval in intervals]
                        intervals = [interval.split('-') for interval in intervals]
                        real_intervals = []
                        for interval in intervals:
                            curr_interval = [interval[-2], interval[-1]]
                            if len(interval) > 2:
                                curr_interval[0] = curr_interval[0]*(-1)
                            real_intervals.append(curr_interval)
                    else:
                        real_intervals = 'null'
                else:
                    real_intervals = 'null'
                
    return tmscore, real_intervals


def get_blades_seqiD(blades):

    if len(blades) > 0:
        id_matrix = [[0 for blade in blades] for blade in blades]
        
        for i, blade_i in enumerate(sorted(list(blades.keys()))):
            for j, blade_j in enumerate(sorted(list(blades.keys()))):
                if i >= j:
                    seq_i = blades[blade_i].strip()
                    seq_j = blades[blade_j].strip()
                    
                    seqID = 0
                    aligned = 0
                    for aa in range(len(seq_i)):
                        if seq_i[aa] != '-' and seq_j[aa] != '-':
                            aligned += 1
                            
                            if seq_i[aa] == seq_j[aa]:
                                seqID += 1

                    if aligned != 0:
                        seqID = seqID*100/aligned
                    else:
                        seqID = 0

                    id_matrix[i][j] = seqID
                    id_matrix[j][i] = seqID

        id_matrix = np.array(id_matrix)
        blades_id = np.median(id_matrix)
    else:
        blades_id = 0.00
    
    return round(blades_id, 2)
        

def get_linker_lengths(blades):

    linker_lengths = []

    all_blades = sorted(list(blades.keys()))
    
    return linker_lengths
            
def mark_blades_in_pdb(pdb_file, blades_intervals):

    new_pdb = '{}_with_{}blades.pdb'.format(pdb_file[:-4], len(blades_intervals))

    previous_bfactor = 0
    coverage = 0
    rescount = 0
    
    with open(new_pdb, 'w') as outpdb:
        with open(pdb_file, 'r') as inpdb:
            for line in inpdb:
                if line.startswith('ATOM'):
                    res_num = line[22:28].strip()
                    bfactor = 0
                    for i, interval in enumerate(blades_intervals):
                        try:
                            if int(res_num) >= int(interval[0]) and int(res_num) <= int(interval[1]):
                                bfactor = i+1
                        except:
                            if res_num == interval[0] or res_num == interval[1]:
                                bfactor = i+1
                            else:
                                bfactor = previous_bfactor

                    previous_bfactor = bfactor
                    
                    if bfactor > -1:
                        outpdb.write('{} {:2}.00\n'.format(line[:60], bfactor))

                        if ' CA ' in line:
                            rescount += 1
                            if bfactor > 0:
                                coverage += 1

                else:
                    outpdb.write(line)
    
    return round(coverage*100/rescount, 2)

# MAIN CODE

def find_blades_in_ecod_folds(ecod_table, ecod_fasta, min_tmscore = 0.7, min_cov = 95.0, exclude_terms = ['beta-pinwheel', 'beta-Prism']):

    # get propellers from ecod table
    propellers = parse_propellers_from_ECOD(ecod_table, exclude_terms = exclude_terms)

    # add their sequence
    propellers = get_propellers_sequence(propellers, ecod_fasta)
    print('There are {} propellers in file {}.'.format(len(propellers.keys()), ecod_table))

    # define output fasta of all considered blades and full sequences of propellers
    #blades_outfasta = 'selected_blades_from_{}.fasta'.format(ecod_fasta[:10])
    full_outfasta = 'selected_full_sequences_from_{}.fasta'.format(ecod_fasta[:-10])

    out_pdb_folder = '{}/selected_pdbs'.format(os.getcwd())
    os.system('mkdir {}'.format(out_pdb_folder))
    
    # for each propeller, find the blades by running symmCE for the pdb file
    pdb_count = 0
    all_blade_sizes = []
    all_blade_IDs = []
    all_blade_tmscores = []
    all_coverages = []
    
    for ecod_code in list(propellers.keys()):

        curr_label = '{}/{}'.format(out_pdb_folder, ecod_code)
        pdb_file = '{}.pdb'.format(curr_label)
        pdb_link = 'http://prodata.swmed.edu/ecod/complete/structure?id={}'.format(ecod_code)

        if not os.path.isfile(pdb_file):
            os.system('wget {} -O {}'.format(pdb_link, pdb_file))
            
        pdb_count += 1
        
        print(' ... Running symmCE for {} ({}/{})'.format(ecod_code, pdb_count, len(propellers)))

        symmCE_fasta = '{}_blades.fasta'.format(curr_label)
        symmCE_stats = '{}_blades.stats'.format(curr_label)

        if not os.path.isfile(symmCE_fasta) or not os.path.isfile(symmCE_stats):
            run_symmCE = sp.Popen([symmCE, pdb_file, '--fasta={}'.format(symmCE_fasta), '--stats={}'.format(symmCE_stats), '--noshow3d'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
            stdout, stderr = run_symmCE.communicate()
    
        blades = parse_symmCE_blades(symmCE_fasta)
        blades_seqID = get_blades_seqiD(blades)
        blades_tm, blades_intervals = parse_symmCE_stats(symmCE_stats)
        
        print(' ... ... Found {} blades with TM-score {} and median sequence identity of {}%'.format(len(blades), blades_tm, blades_seqID))
        
        if blades_intervals != 'null' and blades_tm > min_tmscore:
            print(' ... ... Marking blades in pdb file')
            coverage = mark_blades_in_pdb(pdb_file, blades_intervals)

            print(' ... ... ... Identified repeats cover {}% of the entire domain'.format(coverage))

            if coverage >= min_cov:
                all_blade_sizes.append(len(blades))
                all_blade_IDs.append(blades_seqID)
                all_blade_tmscores.append(blades_tm)
                all_coverages.append(coverage)

                propellers[ecod_code]['blades'] = blades
                propellers[ecod_code]['symm_level'] = blades_tm
                propellers[ecod_code]['seqID'] = blades_seqID
                propellers[ecod_code]['coverage'] = coverage

                with open(full_outfasta, 'a+') as fullout:
                    seq_id = '>{}|{} blades|symm_level: {}|seqID: {}|domain_coverage: {}|f_name: {}'.format(ecod_code, len(blades), blades_tm, blades_seqID, coverage, propellers[ecod_code]['f'])
                    seq = propellers[ecod_code]['seq']

                    fullout.write('{}\n{}\n'.format(seq_id, seq))
            else:
                print(' ... ... No good coverage of the domain! Will eliminate it!')
                os.system('rm {}*'.format(curr_label))
        else:
            print(' ... ... No good repeats found! Will eliminate it!')
            os.system('rm {}*'.format(curr_label))

        json.dump(propellers, open('{}_selected_propellers.json'.format(ecod_fasta[:-10]), 'w'))            

    print("In total, {} propellers were selected".format(len(all_blade_IDs)))

    plt.clf()
    plt.hist(all_blade_sizes, bins = max(all_blade_sizes)-min(all_blade_sizes))
    plt.xlabel('Propeller sizes')
    plt.title('n = {}'.format(len(all_blade_sizes)))
    plt.savefig('selected_propellers_propeller_sizes_hist.pdf', format='pdf')

    plt.clf()
    plt.hist(all_blade_IDs, bins = int(max(all_blade_IDs)-min(all_blade_IDs)))
    plt.xlabel('Blades median sequence ID (%)')
    plt.title('n = {}'.format(len(all_blade_IDs)))
    plt.savefig('selected_propellers_propeller_blades_medianID_hist.pdf', format='pdf')

    plt.clf()
    plt.hist(all_blade_tmscores, bins = int((max(all_blade_tmscores)-min(all_blade_tmscores))*100))
    plt.xlabel('Blades TMscore')
    plt.title('n = {}'.format(len(all_blade_tmscores)))
    plt.savefig('selected_propellers_propeller_blades_tmscore_hist.pdf', format='pdf')

    plt.clf()
    plt.scatter(all_blade_IDs, all_blade_tmscores)
    plt.xlabel('Blades median sequence ID (%)')
    plt.ylabel('Blades TMscore')
    plt.title('n = {}'.format(len(all_blade_IDs)))
    plt.savefig('selected_propellers_propeller_blades_medianID_vs_TMscore.pdf', format='pdf')

    plt.clf()
    plt.hist(all_coverages, bins = int((max(all_coverages)-min(all_coverages))))
    plt.xlabel('Domain coverage %')
    plt.title('n = {}'.format(len(all_blade_tmscores)))
    plt.savefig('selected_propellers_propeller_blades_coverages_hist.pdf', format='pdf')
    
    return propellers

#######

#ecod_table = 'ecod.latest.F70.domains.txt'
ecod_table = 'ecod.latest.F70.domains_prisms_and_pinwheels.txt'
# ecod_table = 'ecod.develop261.F99.domains.txt'
# ecod_fasta = 'ecod.develop261.F99.fasta.txt'
#propellers = find_blades_in_ecod_folds(ecod_table, ecod_fasta, min_tmscore = 0.0, min_cov = 0.0, exclude_terms = ['beta-pinwheel', 'beta-Prism']) # for beta propellers only
propellers = find_blades_in_ecod_folds(ecod_table, ecod_fasta, min_tmscore = 0.0, min_cov = 0.0, exclude_terms = []) # for pinwheels and prisms
json.dump(propellers, open('{}_selected_propellers.json'.format(ecod_fasta[:-10]), 'w'))

