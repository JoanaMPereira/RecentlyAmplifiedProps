import multiprocessing as mp
import subprocess as sp
import numpy as np

import argparse
import time
import os

wrap_searches = 'wrap_searches.py'

# Get inputs
parser = argparse.ArgumentParser(description='')

# ... Required inputs
parser.add_argument('in_folder', metavar='in_folder', type=str, help='folder where the representatives are')
parser.add_argument('databases', metavar='databases', nargs='+', help='the databases to search')
# ... Optional inputs
parser.add_argument('-cpu', dest='cpu', default = 2,type=int, help='number of cpus to use (default: 2)')
parser.add_argument('-parallel_jobs', dest='parallel_jobs', default = 1,type=int, help='number of jobs to do in parallel (default: 1)')
parser.add_argument('-psiblast_rounds', dest='psiblast_rounds', nargs='+', default = 10, type=int, help='number of rounds to run psiblast (default: 10)')

args = parser.parse_args()
curr_directory = os.getcwd()

## HELPING ROUTINES

def chunk_list(l, n):
    chunks = np.array_split(np.array(l), n)
    chunks = [list(chunk) for chunk in chunks]
    return chunks

def run_wrap_searches(arguments):

    job_id = arguments[0]
    folders = arguments[1]
    databases = arguments[2]
    job_cpus = arguments[3]
    psiblast_rounds = arguments[4]

    print("Job {}. Running searches".format(job_id))
    
    for count, curr_folder in enumerate(folders):

        cycle_start = time.time()

        print(" ... Job {}. Running searches for folder {} ({}/{})".format(job_id, curr_folder.split('/')[-1], count, len(folders)))

        in_fasta = ['{}/{}'.format(curr_folder, file) for file in os.listdir(curr_folder) if 'representative_blades.fasta' in file][0]
        
        command = ['python3', wrap_searches, in_fasta]
        for db in databases:
            command.append(db)
        command.append('-cpu')
        command.append(str(job_cpus))
        command.append('-working_dir')
        command.append(curr_folder)
        command.append('-psiblast_rounds')
        for pbround in psiblast_rounds:
            command.append(str(pbround))

        run_wrap = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_wrap.communicate()

        with open('{}/searches_log.log'.format(curr_folder), 'w') as out:
            out.write(stdout.decode('ascii'))

        cycle_end = time.time()
        numb_seconds = cycle_end - cycle_start
        print(" ... Job {}. Cycle {}. Finished after: {} ".format(job_id, count, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))


## MAIN CODE

# 1. get all cases in the input folder
list_of_folders_for_jobs = ['{}/{}/{}'.format(curr_directory, args.in_folder, folder) for folder in os.listdir('{}/{}'.format(curr_directory, args.in_folder)) if os.path.isdir('{}/{}/{}'.format(curr_directory, args.in_folder, folder))]
list_of_folders_for_jobs.reverse()

print("There are a total of {} jobs to run".format(len(list_of_folders_for_jobs)))

# 2. Distribute jobs and find number of cpus per job
number_job_cpus = int(args.cpu/args.parallel_jobs)
print(" ... Will run {} jobs in parallel, each using {} cpus".format(args.parallel_jobs, number_job_cpus))
print(" ... That is a total of circa {} cycles".format(round(len(list_of_folders_for_jobs)/args.parallel_jobs, 1)))

separated_jobs = chunk_list(list_of_folders_for_jobs, args.parallel_jobs)
list_arguments = [i for i in zip(range(args.parallel_jobs), separated_jobs, [args.databases for job in separated_jobs], [number_job_cpus for job in separated_jobs], [args.psiblast_rounds for job in separated_jobs])]

start = time.time()
pool = mp.Pool(number_job_cpus)
results = pool.map(run_wrap_searches, list_arguments)

pool.close()
pool.join()

end = time.time()
numb_seconds = end - start
print(" ... Finished after: {} ".format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

