import numpy
import math
import sys 
import os
import matplotlib
import ast
import json

from Bio import SeqIO
from Bio import Entrez
from pylab import *
from pdbe import pyPDBeREST

p = pyPDBeREST()
#Entrez.email = 'joana.pereira@tuebingen.mpg.de'
#Entrez.email = 'pereira.joanam@gmail.com'

sequencefile = sys.argv[1]

## ROUTINES ##

def get_ncbids_from_fasta(fastafile):

   ncbids = []

   with open(fastafile, 'r') as infile:
      for line in infile:
         if line.startswith('>'):
            line_data = line.split('|')
            ncbid = line_data[0][1:]
            ncbids.append(ncbid)

   return ncbids

def get_dictionaries_from_fasta(fastafile):

   taxonomy = {}

   with open(fastafile, 'r') as infile:
      for line in infile:
         if line.startswith('>'):
            line_data = line.split('|')
                   
            entrezID = line_data[0][1:]
            superkingdom = line_data[1]
            phylum = line_data[2]
            taxclass = line_data[3]
            order = line_data[4]
            genus = line_data[5]
            species = line_data[6]
            prot_name = line_data[7]
            intervals = line_data[8].split(': ')[-1]
            evalues = line_data[9].split(': ')[-1]

            if superkingdom not in taxonomy.keys():
                taxonomy[superkingdom] = {}

            if phylum not in taxonomy[superkingdom].keys():
                taxonomy[superkingdom][phylum] = {}

            if taxclass not in taxonomy[superkingdom][phylum].keys():
                taxonomy[superkingdom][phylum][taxclass] = {}

            if order not in taxonomy[superkingdom][phylum][taxclass].keys():
                taxonomy[superkingdom][phylum][taxclass][order] = {}

            if genus not in taxonomy[superkingdom][phylum][taxclass][order].keys():
                taxonomy[superkingdom][phylum][taxclass][order][genus] = {}
                 
            if species not in taxonomy[superkingdom][phylum][taxclass][order][genus].keys():
                taxonomy[superkingdom][phylum][taxclass][order][genus][species] = {'ncbi_codes': [], 'Intervals': [], 'Evalues': []}

            taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Intervals'].append(intervals)
            taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Evalues'].append(evalues)
            taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes'].append(entrezID)

   return taxonomy 
   
## MAIN CODE ##

def get_sequences_taxonomy_tree(sequencefile):

    outpfile = '{}_taxonomy'.format(sequencefile)
    tax_outjsonfile = '{}_taxonomy.json'.format(sequencefile)

    if os.path.isfile(outpfile):
       ncbids_previously_collected = set(get_ncbids_from_fasta(outpfile))
       taxonomy = get_dictionaries_from_fasta(outpfile)
    else:
       ncbids_previously_collected = set([])
       taxonomy = {}

    print("Starting collecting taxonomy from sequence number {}".format(len(ncbids_previously_collected)))
    
    for protsequence in SeqIO.parse(sequencefile, u"fasta"):

        if 'pdb' not in protsequence.description and '|sp|' not in protsequence.description:
           
           seq_header = str(protsequence.description).split("|")

           entrezID = seq_header[0]           
           species = seq_header[1]
           if species.split()[0] == 'uncultured':
              genus = species.split()[1]
           else:
              genus = species.split()[0]
              
           prot_name = seq_header[2]        
           medians = seq_header[3].split(':')[-1][1:]
           propellers = seq_header[4].split(':')[-1][1:]
           intervals = seq_header[5].split(':')[-1][1:]

           seq_length = len(protsequence.seq)
           
           if entrezID not in ncbids_previously_collected:
              
              try:
                  if '_#' in entrezID:
                     ncbi_code = entrezID.split('_#')[0]
                  else:
                     ncbi_code = entrezID
                     
                  handle = Entrez.esummary(db="protein", id=ncbi_code)   
                  record = Entrez.read(handle)

                  #try:
                  taxID = record[0]["TaxId"]
                  taxsearch = Entrez.efetch(id = taxID, db = "taxonomy", retmode = "xml")
                  taxrecords = Entrez.parse(taxsearch)

                  for taxrecord in taxrecords:
                      taxrecord = taxrecord

                  superkingdom = 'na'
                  phylum = 'na'
                  taxclass = 'na'
                  order = 'na'
                   
                  for level in taxrecord[u"LineageEx"]:
                     if level[u"Rank"] == "superkingdom":
                         superkingdom = level[u"ScientificName"]
                     elif level[u"Rank"] == "phylum":
                         phylum = level[u"ScientificName"]
                     elif level[u"Rank"] == "class":
                         taxclass = level[u"ScientificName"]
                     elif level[u"Rank"] == "order":
                         order = level[u"ScientificName"]

                  if superkingdom not in taxonomy.keys():
                      taxonomy[superkingdom] = {}
                       
                  if phylum not in taxonomy[superkingdom].keys():
                      taxonomy[superkingdom][phylum] = {}

                  if taxclass not in taxonomy[superkingdom][phylum].keys():
                      taxonomy[superkingdom][phylum][taxclass] = {}

                  if order not in taxonomy[superkingdom][phylum][taxclass].keys():
                      taxonomy[superkingdom][phylum][taxclass][order] = {}

                  if genus not in taxonomy[superkingdom][phylum][taxclass][order].keys():
                      taxonomy[superkingdom][phylum][taxclass][order][genus] = {}
                       
                  if species not in taxonomy[superkingdom][phylum][taxclass][order][genus].keys():
                      taxonomy[superkingdom][phylum][taxclass][order][genus][species] = {'ncbi_codes': [], 'Intervals': [], 'Medians': [], 'Propellers': []}

                  taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Intervals'].append(intervals)
                  taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Medians'].append(medians)
                  taxonomy[superkingdom][phylum][taxclass][order][genus][species]['Propellers'].append(propellers)
                  taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes'].append(entrezID)
                   
                  print('{} {} {} {} {} {} {}'.format(entrezID, superkingdom, phylum, taxclass, order, genus, species))

                  protID = '>{}|{}|{}|{}|{}|{}|{}|{}| Intervals: {} | Medians: {} | Propellers: {}'.format(entrezID, superkingdom, phylum, taxclass, order, genus, species, prot_name, intervals, medians, propellers)

                  outp = open(outpfile, 'a+')
                  outp.write('{}\n{}\n'.format(protID, str(protsequence.seq)))
                  outp.close()

                  outjson = open(tax_outjsonfile, "w")
                  json.dump(taxonomy, outjson, indent = 4)
          
   ##               except IndexError:
   ##                   print("Can't get taxonomy for '{}'; can't get taxrecord".format(entrezID))

              except RuntimeError:
                  print("Can't get taxonomy for '{}'; can't get handle.".format(entrezID))


    return taxonomy

def build_taxonomy_distribution_summary(taxonomy):

    summary_distribution = {}

    for superkingdom in taxonomy.keys():
       if superkingdom not in summary_distribution.keys():
          summary_distribution[superkingdom] = {}
          
       for phylum in taxonomy[superkingdom].keys():
          if phylum not in summary_distribution[superkingdom].keys():
             summary_distribution[superkingdom][phylum] = {}

          for taxclass in taxonomy[superkingdom][phylum].keys():
             if taxclass not in summary_distribution[superkingdom][phylum].keys():
                summary_distribution[superkingdom][phylum][taxclass] = {}

             for order in taxonomy[superkingdom][phylum][taxclass].keys():
                if order not in summary_distribution[superkingdom][phylum][taxclass].keys():
                   summary_distribution[superkingdom][phylum][taxclass][order] = {}

                for genus in taxonomy[superkingdom][phylum][taxclass][order].keys():
                   if genus not in summary_distribution[superkingdom][phylum][taxclass][order].keys():
                      summary_distribution[superkingdom][phylum][taxclass][order][genus] = {'Species': [], 'ncbi_codes': []}

                   for species in taxonomy[superkingdom][phylum][taxclass][order][genus].keys():

                      summary_distribution[superkingdom][phylum][taxclass][order][genus]['Species'].append('{} ({})'.format(species, len(taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes'])))

                      for ncbi_code in taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes']:
                         summary_distribution[superkingdom][phylum][taxclass][order][genus]['ncbi_codes'].append(ncbi_code)                   

    outjsonfile = '{}_taxonomy_summary_distribution.json'.format(sequencefile)
    outjson = open(outjsonfile, "w")
    json.dump(summary_distribution, outjson, indent = 4)           
             
    return summary_distribution

      
# RUN MAIN ROUTINES

taxonomy_mapped_tree = get_sequences_taxonomy_tree(sequencefile)
summary_distribution = build_taxonomy_distribution_summary(taxonomy_mapped_tree)
