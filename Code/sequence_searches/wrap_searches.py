import subprocess as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import time
import os

# Get inputs
parser = argparse.ArgumentParser(description='')

# ... Required inputs
parser.add_argument('in_fasta', metavar='in_fasta', type=str, help='the fasta file with the seed sequences to search')
parser.add_argument('databases', metavar='databases', nargs='+', help='the databases to search')
# ... Optional inputs
parser.add_argument('-cpu', dest='cpu', default = 2,type=int, help='number of cpus to use (default: 2)')
parser.add_argument('-psiblast_rounds', dest='psiblast_rounds', nargs='+', default = 10, type=int, help='number of rounds to run psiblast (default: 10)')
parser.add_argument('-working_dir', dest='working_dir', type=str, default = '.')

args = parser.parse_args()
if args.working_dir == '.':
    curr_directory = os.getcwd()
else:
    curr_directory = args.working_dir

folder_with_dependencies = '/ebio/abt1_share/highly_repetitive_propellers/data'

## HELPING ROUTINES

def get_hit_sequence_files(run_psiblast_search_out):

    aligned_sequences = ''
    full_sequences = ''

    number_hits = 0
    number_seq = 0

    with open(run_psiblast_search_out, 'r') as psiblast_out:
        for line in psiblast_out:
                if 'aligned_only' in line:
                    aligned_sequences = line.split()[-1]
                    number_hits = int(line.split()[-5])
                elif 'full_sequence' in line:
                    full_sequences = line.split()[-1]
                    number_seq = int(line.split()[-5])

    return [aligned_sequences, number_hits], [full_sequences, number_seq]

def parse_propeller_histogram(hist_file):

    proportion_2bladed = 0
    largest_propeller = 0
    most_common_propeller = 0
    most_common_count = 0

    propeller_sizes = []
    propeller_counts = []

    with open(hist_file, 'r') as histin:
        for line in histin:
            if len(line) > 1:
                size = int(line.split('\t')[0])
                count = int(line.split('\t')[1])

                if size == 2:
                    proportion_2bladed = count

                if count > most_common_count:
                    most_common_propeller = size
                    most_common_count = count
                    
                propeller_sizes.append(size)
                propeller_counts.append(count)

    if len(propeller_sizes) > 0:
        largest_propeller = max(propeller_sizes)
        proportion_2bladed = proportion_2bladed/sum(propeller_counts)

    
    return proportion_2bladed, largest_propeller, most_common_propeller

def make_plot(data, x, out_file = 'plot.pdf', ylabel = 'value'):

    df = pd.DataFrame(data)
    df['Number of psiblast rounds'] = x
    df = pd.melt(df, id_vars = ['Number of psiblast rounds'])
    df.columns = ['Number of psiblast rounds', 'Database', 'value']

    if len(df) > 0 and set(df.value) != set([0]):
        plt.clf()
        sns.barplot(x = 'Number of psiblast rounds', y= 'value', hue = 'Database', data = df)
        plt.xlabel('Number of psiblast rounds collected')
        plt.ylabel(ylabel)
        plt.savefig(out_file, format='pdf')
    else:
        print('There is no data, will not plot anything')
    
## MAIN CODE

hard_start = time.time()

# collect the evolution of the searches in each database depending on the number of rounds collected (get the proportion of 2 bladed propellers and the size of the largest propeller)
proportion_2bladeds = {db: [] for db in args.databases}
largest_propellers = {db: [] for db in args.databases}
most_common_propellers = {db: [] for db in args.databases}

for db in args.databases:

    print("\nPREFORMING SEARCHES FOR {} DATABASE\n".format(db))

    # Make job-specific folder
    job_dir = "{}/{}".format(curr_directory, db)
    os.system("mkdir {}".format(job_dir))

    round_start = time.time()
    
    for psiblast_rounds in sorted(args.psiblast_rounds):
        
        working_directory = '{}/{}rounds'.format(job_dir, psiblast_rounds)
        os.system("mkdir {}".format(working_directory))

        psiblast_rounds_to_run = max(args.psiblast_rounds)
        psiblast_rounds_to_collect = psiblast_rounds
        
        main_start = time.time()
        
        print("\n ... DOING SEARCHES FOR {} PSIBLAST ROUNDS\n".format(psiblast_rounds))
        print(" ... ... 1. Running psiblast search over {} database and collecting hits".format(db))

        start = time.time()
        
        run_psiblast_search = sp.Popen(['nice', 'python', '{}/run_psiblast_for_fasta.py'.format(folder_with_dependencies), args.in_fasta, 'both', db,  str(psiblast_rounds_to_run), str(psiblast_rounds_to_collect), str(args.cpu), working_directory], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_psiblast_search.communicate()
        print(stderr)
        
        with open('{}/psiblast_search_{}_database.log'.format(working_directory, db), 'w') as out_log:
            out_log.write(stdout.decode('ascii'))

        [aligned_sequences, number_hits], [full_sequences, number_seq] = get_hit_sequence_files('{}/psiblast_search_{}_database.log'.format(working_directory, db))

        end = time.time()
        numb_seconds = end - start
        print(" ... ... ... Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

        if number_hits > 1:
            
            print(" ... ... 2. Searching for missed blades")

            start = time.time()
                
            find_missed_blades = sp.Popen(['nice', 'python3', '{}/search_missing_blades.py'.format(folder_with_dependencies), aligned_sequences, full_sequences, '-cpu', str(args.cpu), '-label', db, '-working_directory', working_directory], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
            stdout, stderr = find_missed_blades.communicate()
            print(stderr)
            
            with open('{}/search_for_missing_blades_{}_database.log'.format(working_directory, db), 'w') as out_log:
                out_log.write(stdout.decode('ascii'))

            end = time.time()
            numb_seconds = end - start
            print(" ... ... ... Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

            print(" ... ... 3. Finding propellers and selecting relevant sequences")

            start = time.time()
            
            find_propellers = sp.Popen(['nice', 'python3', '{}/select_relevant_sequences.py'.format(folder_with_dependencies), '{}/missing_blade_searching_{}/sequences_with_missed_terminal_blades.json'.format(working_directory, db), full_sequences, '-cpu', str(args.cpu), '-label', db, '-working_directory', working_directory], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
            stdout, stderr = find_propellers.communicate()
            print(stderr)

            with open('{}/selecting_for_relevant_sequences_{}_database.log'.format(working_directory, db), 'w') as out_log:
                out_log.write(stdout.decode('ascii'))

            end = time.time()
            numb_seconds = end - start
            print(" ... ... ... Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

            selected_fasta = '{}/selection_of_relevant_sequences_{}/selected_relevant_sequences_full_length.fasta'.format(working_directory, db)

            if os.path.isfile(selected_fasta):
                print(" ... ... 4. Mapping taxonomy to selected sequences")

                start = time.time()

                map_taxonomy = sp.Popen(['nice', 'python', '{}/map_taxonomy_to_selected.py'.format(folder_with_dependencies), selected_fasta], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
                stdout, stderr = map_taxonomy.communicate()
                print(stderr)

                with open('{}/mapping_taxonomy_to_selected_{}_database.log'.format(working_directory, db), 'w') as out_log:
                    out_log.write(stdout.decode('ascii'))

                end = time.time()
                numb_seconds = end - start
                print(" ... ... ... Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

                print(" ... ... 5. Collecting round summary measures (proportion of 2bladed propellers and size of largest propeller)")

                if os.path.isfile('{}/selection_of_relevant_sequences_{}/relevant_propellers_histogram.txt'.format(working_directory, db)):
                    proportion_2bladed, largest_propeller, most_common_propeller = parse_propeller_histogram('{}/selection_of_relevant_sequences_{}/relevant_propellers_histogram.txt'.format(working_directory, db))
                    proportion_2bladeds[db].append(proportion_2bladed)
                    largest_propellers[db].append(largest_propeller)
                    most_common_propellers[db].append(most_common_propeller)
                else:
                    proportion_2bladeds[db].append(0)
                    largest_propellers[db].append(0)
                    most_common_propellers[db].append(0)
            else:
                print("\n ... ### None of the hits in database {} is highly repetitive. Sorry... Will search on the next one".format(db))
                proportion_2bladeds[db].append(0)
                largest_propellers[db].append(0)
                most_common_propellers[db].append(0)
            
            print(" ... DONE with {} rounds (over {})!".format(psiblast_rounds, db))

            main_end = time.time()
            numb_seconds = main_end - main_start
            print(" ... ... Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

        else:
            print("\n ... ### Did not find any hits in database {}. Sorry... Will search on the next one".format(db))
            proportion_2bladeds[db].append(0)
            largest_propellers[db].append(0)
            most_common_propellers[db].append(0)
            
    round_end = time.time()
    numb_seconds = round_end - round_start
    print(" ... Finished {} database search after: {} days {}".format(db, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

if len(most_common_propellers) > 0:
    print("\nPLOTTING THE EVOLUTION OF THE SEARCHES DEPENDING ON THE NUMBER OF ROUNDS COLLECTED\n")
    make_plot(proportion_2bladeds, sorted(args.psiblast_rounds), ylabel = 'Proportion of 2-bladed propellers', out_file = '{}/proportion_of_2bladed_propellers.pdf'.format(curr_directory))
    make_plot(largest_propellers, sorted(args.psiblast_rounds),  ylabel = 'Largest propeller size', out_file = '{}/largest_propeller_per_rounds_collected.pdf'.format(curr_directory))
    make_plot(most_common_propellers, sorted(args.psiblast_rounds),  ylabel = 'Most common propeller size', out_file = '{}/most_common_propeller_per_rounds_collected.pdf'.format(curr_directory))

hard_end = time.time()
numb_seconds = hard_end - hard_start
print("\n##### DONE WITH EVERYTHING: Finished after: {} days {}".format(int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
