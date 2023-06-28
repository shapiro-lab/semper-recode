

from os import listdir
from os.path import isfile, join
from os import path
import pandas as pd
import numpy as np
import re

from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import parse, write
from Bio import AlignIO
from Bio.SeqFeature import ExactPosition, SeqFeature, FeatureLocation

#################PACKAGES OF INTEREST ##############################
""" 
HMMER
HOMER
SAMTOOLS
Salmon
Trinity
Bedtools
Clustalo
Mappy
Minimap
Fuzzy
MMSEQS2
anti-smash
Chop Chop

CookieCutter
venv
virtualenv
"""



def get_paths(folder_path):
    """
    Given a folder as path, outputs lits of all paths in that folder

    Parameters
    ----------
    folder_path : str
        path of folder for which you want.

    Returns
    -------
    list of paths of items in that folder
    """
    if folder_path[-1]=='/':
        return [folder_path + f for f in listdir(folder_path) if isfile(join(folder_path, f))]
    else:
        return [folder_path + '/' + f for f in listdir(folder_path) if isfile(join(folder_path, f))]



def dataframe_from_fasta(fasta_path, case):
    """
    takes in a fasta file and then outputs the fasta file and metadata as a dataframe.
    
    Parameters
    ----------
    fasta_path : str
        path of the fasta file you want to parse to dataframe
    case : str
        "uniprot"
    
    Returns
    -------
    dataframe for parsed fasta file
    """
    
    # raise errors
    if case not in ['uniprot', 'general', 'uniref', 'swissprot']:
        raise ValueError(
                "Second input arugment (case) must be \"uniprot\" or \"general\" "
            )

    seqs_df = pd.DataFrame()
    for hit in SeqIO.parse(fasta_path, "fasta"):
        if case in ['uniprot', 'swissprot']:
            line_item = {'fasta-id' : hit.id,
                         'accession' : re.findall('[\d\w]{1,4}\|([\d\w]{1,})\|', hit.id)[0],
                         'sequence' : str(hit.seq),
                         'length' : int(len(str(hit.seq))),
                         'name' : hit.name,
                         'description' : re.findall('[_\d\w\s\|] (.*)', hit.description)[0],
                         'OS':re.findall('OS=(.*) OX', hit.description)[0] if re.findall('OS=(.*) OX', hit.description) else None ,
                         'OX':re.findall('OX=(.*) GN', hit.description)[0] if re.findall('OX=(.*) GN', hit.description) else None,
                         'GN':re.findall('GN=(.*) PE', hit.description)[0] if re.findall('GN=(.*) PE', hit.description) else None,
                         'PE':re.findall('PE=(.*) SV', hit.description)[0] if re.findall('PE=(.*) SV', hit.description) else None,
                         'SV':re.findall('SV=(.*)', hit.description)[0] if re.findall('SV=(.*)', hit.description) else None,
                        }
            
        elif case =="uniref":
            #print(hit.description)
            line_item = {'fasta-id' : hit.id,
                         'accession' : hit.id,
                         'sequence' : str(hit.seq),
                         'length' : int(len(str(hit.seq))),
                         'name' : hit.name,
                         'description' : re.findall('UniRef100_[\d\w]{1,} (.*) n=', hit.description)[0],
                         'Tax':re.findall('Tax=(.*) TaxID', hit.description)[0] if re.findall('Tax=(.*) TaxID', hit.description) != [] else None,
                         'TaxID':re.findall('TaxID=(.*) RepID', hit.description)[0] if re.findall('TaxID=(.*) RepID', hit.description) != [] else None,
                         'RepID':re.findall('RepID=(.*)', hit.description)[0] if re.findall('RepID=(.*)', hit.description) !=[] else None,
                        }
            
        
        elif case == 'general':
            line_item = {'fasta-id' : hit.id,
                        'sequence' : str(hit.seq),
                        'name' : hit.name,
                        'description' : hit.description ,
                            }
            
        seqs_df = seqs_df.append(line_item, ignore_index = True)
    return seqs_df


def dataframe_to_fasta(df, id_column, sequence_column, description_columns, output_path):
    """converts a dataframe of sequence information to a fasta file

    df : pd.DataFrame
        pandas DataFrame containing id, sequence, and description information

    id_column : str
        the name of the column that has the sequences identifications in it.

    sequence_column : str
        the name of the column that has the sequences in it. sequences should be strings

    description_columns : list
        list of column names of useful information you would like to be in the fasta description for the sequence
    """
    # get the seqrecords
    seq_records = []
    for index, row in df.iterrows():
        seq_records.append(SeqRecord(
                                    Seq(row[sequence_column]),
                                    id=row[id_column],
                                    description=" ".join([str(x + "=" + str(row[x])) for x in description_columns])
                                    )
                          )

    #write the seqrecords to file
    with open(output_path, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")
    return



##### HMMER ######
'''
def run_jackhmmer(seeds_paths, database_name, database_path, output_path):
    for s in seeds_paths:

        file_name = re.findall('/([\d\s\w]{1,}).fasta', s)[0] + "_jackhmmer_summary_" + database_name
        summary_out_path = output_path +file_name + ".out"
        tblout_out_path = output_path +file_name + ".tblout"

        jackhmmer_command = "jackhmmer -o " + summary_out_path + " --tblout " + tblout_out_path + " --acc --noali --notextw -E 0.1 --cpu 6 " + s + " " + database_path
        print(jackhmmer_command)
        ! $jackhmmer_command 
    return


##### MMSEQS2 ########

def mmseqs2_test_min_seq(input_path, output_directory, file_handle, threads, v, cov_mode, min_seq_ids):

    """
    tests different min_seq_id values and plots them for viewing. 
    
    Parameters
    ----------
    input_path : str
        the path to the fasta file with your sequences to cluster
    
    output_directory : str
        the output directory path where you want all these documents to live. must end with "/"
        
    
    file_handle : str
        the name you want all the resulting files to share
    
    threads : int
        number of threads you want to use to run this
    
    v : int
        how verbose do you want the printed command line output to be
        verbosity level: 0=nothing, 1: +errors, 2: +warnings, 3: +info
    
    cov_mode : int
        0: coverage of query and target - used to cluster full length protein sequences
        1: coverage of target - can be used to cluster protein fragments
        2: coverage of query - The query coverage mode can be used while searching e.g. to assure a certain level of coverage.
        3: target seq. length needs be at least x% of query length - 
        4: query seq. length needs be at least x% of target length - 
        
    (NOT INCLUDED) c : int
        list matches above this fraction of aligned (covered) residues (see --cov-mode)
    
    min_seq_id : list
        list of floats from [0.0,1.0] to test: min_seq_id in mmseqs is used list matches above this sequence identity (for clustering)
        
        
    Returns
    -------
    a plot of the number of cluster for each --min-seq-id in in min-seq-ids
    """
    
    if output_directory[-1] != "/": 
        raise ValueError(
                        "output_directory parameter must end with \"/\"")
    
    no_clusters = []
    for msi in min_seq_ids:
        # make a sub directory for each msi we will use
        output_sub_dir = output_directory + "cluster_" + str(msi) + "/"
        ! mkdir $output_sub_dir
        
        # define some vars for the functions
        path = input_path
        db_out_path = output_sub_dir + file_handle + ".db"
        cluster_out_path = output_sub_dir + file_handle + "_clusters.db"
        tsv_out_path =  output_sub_dir + file_handle+ "_clusters.tsv"
        
        # run mmseqs2 for each msi
        ! mmseqs createdb $path $db_out_path
        ! mmseqs cluster $db_out_path $cluster_out_path /tmp --v $v --threads $threads --cov-mode $cov_mode --min-seq-id $msi
        ! mmseqs createtsv --full-header true $db_out_path $db_out_path $cluster_out_path $tsv_out_path
        
        print(mmseqs cluster $db_out_path $cluster_out_path /tmp -v $v --threads $threads --cov-mode $cov_mode --min-seq-id $msi)

        # collect useful clustering information
        db = pd.read_csv(tsv_out_path, sep="\t", header = None)
        new_db = db.groupby(by=0)[1].apply(list).to_frame().reset_index()
        new_db = new_db.rename(columns={0: "centroid", 1: "members"})
        new_db['count-members'] = ''
        
        for index, row in new_db.iterrows():
            new_db.at[index,'count-members'] = len(row['members'])
        new_db = new_db[[ 'centroid', 'members', 'count-members']].sort_values(by=['count-members'], ascending = False)
        new_db.to_csv(output_sub_dir + file_handle + "_csv_summary.csv")
        
        # add the number of clusters to our list
        no_clusters.append(len(new_db))
        
    plt.scatter(min_seq_ids, no_clusters)
    plt.xlabel("min-seq-id")
    plt.ylabel("number of clusters")
    
    return 


def parse_hmmer_results(tblout_path, output_path, regex):
    """ Function to parse the tblout files from hmmer into csvs and then save the csvs. 
    
    Parameters
    ----------
    tblout_path: str 
        output path of tblout file
    
    output_path: str
        path to where you want to store the parsed csv
        
        
    Returns
    -------
    
    df_hits: pd.DataFrame
        Dataframe object with the parsed csvs. 
    """
    hit_data = []
    df_cols = [
        "protein-id",
        "e-value",
        "bitscore",
        "description",
        "taxon-id"
    ]


    results = SearchIO.to_dict(SearchIO.parse(open(tblout_path, "r"), "hmmer3-tab"))
    for queryset in results.values():
        for hit in queryset:
            # Remove hits with more than 1 hit domain per protein
            #if hit.domain_included_num != 1: continue
            
            if re.findall('TaxID=([\d]{1,})',hit.description) != []:
                taxon_id = re.findall('TaxID=([\d]{1,})',hit.description)[0]
            else:
                taxon_id = ''
            
            hit_data.append({
                "protein-id": hit.id,
                "e-value": hit.evalue,
                "bitscore": hit.bitscore,
                "description": hit.description,
                "taxon-id": taxon_id
                
            })
    df_hits = pd.DataFrame(hit_data, columns=df_cols)
    
    df_hits['protein-accession'] = ''
    
    for index, row in df_hits.iterrows():
        accession = re.findall(regex, row['protein-id'])[0]
        df_hits.at[index, 'protein-accession'] = accession
        
    df_hits.to_csv(output_path)
    print("Found " + str(len(df_hits)) +" hits.")
    
    return df_hits
'''