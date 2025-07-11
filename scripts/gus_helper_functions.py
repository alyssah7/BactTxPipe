#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing helper scripts for GUS
"""
from Bio import SeqIO
from itertools import compress
import os
import numpy as np
import pickle
import csv
import glob
import subprocess
import gzip

def read_samples_CSV(spls):
    hdr_check = ['Path','Sample','FileName','Reference','Group','Outgroup']
    switch = "on"
    file = open(spls, 'r')
    #initialize lists
    list_path = []; list_splID = []; list_fileN = []; list_refG = []
    list_group=[]; list_outgroup=[]
    for line in file:
        line = line.strip('\n').split(',')
        # Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
        if switch == "on":
            if (line == hdr_check):
                print("Passed CSV header check")
            else:
                Warning("CSV did NOT pass header check! Code continues, but first line ignored")
            switch = "off"
            continue
        # build lists
        list_path.append(line[0])
        list_splID.append(line[1])
        list_fileN.append(line[2])
        try:
            list_refG.append(line[3])
        except Exception as e:
            list_refG.append(None)
        try:
            list_group.append(line[4])
        except Exception as e:
            list_group.append(None)
        try:
            list_outgroup.append(line[5])
        except Exception as e:
            list_outgroup.append(None)

    return [list_path,list_splID,list_fileN,list_refG,list_group,list_outgroup]


def split_samplesCSV(PATH_ls,SAMPLE_ls,FILENAME_ls,REF_Genome_ls,GROUP_ls,OUTGROUP_ls):
    '''Take info from samples.csv, concat by sample name & save each line as sample_info.csv in data/{sampleID}'''
    
    #Loop through unique samples
    for sample in set(SAMPLE_ls):
        # Concat info for this sample
        sample_info_bool=[s==sample for s in SAMPLE_ls]
        sample_paths=list(compress(PATH_ls,sample_info_bool))
        sample_filenames=list(compress(FILENAME_ls,sample_info_bool))
        sample_references=list(compress(REF_Genome_ls,sample_info_bool))
        sample_groups=list(compress(GROUP_ls,sample_info_bool))
        sample_outgroups=list(compress(OUTGROUP_ls,sample_info_bool))
        
        sample_info_csv=list(zip(sample_paths,[sample]*sum(sample_info_bool),sample_filenames,sample_references,sample_groups,sample_outgroups))
        
        # make data directory for this sample if it doesn't already exist
        if not(os.path.isdir('data/' + sample)):
            os.makedirs('data/' + sample, exist_ok=True)
        # check to see if this mini csv with sample info already exists
        if os.path.isfile('data/' + sample + '/sample_info.csv'):
            # if so, read file
            old_file = open('data/' + sample + '/sample_info.csv','r')
            old_info_read = csv.reader(old_file)
            old_info = list(map(tuple, old_info_read))
            old_file.close()
            
            # check to see if the existing file is consistent with samples.csv
            if not(old_info == sample_info_csv):
                # if not, remove the old file and save sample info in a new file
                # print('Information file must be updated.')
                os.remove('data/' + sample + '/sample_info.csv')
                with open('data/' + sample + '/sample_info.csv', "w") as f:
                    writer = csv.writer(f)
                    for row in sample_info_csv:
                        writer.writerow(row)
        else: # if mini csv with sample info does not already exist
            # save sample info in mini csv
            with open('data/' + sample + '/sample_info.csv', "w") as f:
                writer = csv.writer(f)
                for row in sample_info_csv:
                    writer.writerow(row)


def read_sample_info_CSV(path_to_sample_info_csv):
    with open(path_to_sample_info_csv,'r') as f:
        this_sample_info = f.readline() # only one line to read
    this_sample_info = this_sample_info.strip('#').split(',')
    path = this_sample_info[0] # remember python indexing starts at 0
    paths = path.split(' ')
    sample = this_sample_info[1]
    filename = this_sample_info[2]
    reference = this_sample_info[3]
    
    return paths, sample, reference, filename


def findfastqfile(dr,ID,filename):
    fwd=[]
    rev=[]
    #search for filename as directory first
    # AH added wildcard in '/*' because of weird vedanta file naming thing (2024/04/09)
    potentialhits_forward=glob.glob(dr + '/*' + filename +'/*1.fq.gz')
    potentialhits_reverse=glob.glob(dr + '/*' + filename +'/*2.fq.gz')
    if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
        fwd=potentialhits_forward[0]
        rev=potentialhits_reverse[0]
    #then search for filename as file.gz
    elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
        potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fq.gz')
        potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fq.gz')
        if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
            fwd=potentialhits_forward[0]
            rev=potentialhits_reverse[0]
        #then search as unzipped file
        elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
            potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq')
            potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq')
            if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                subprocess.run("gzip " + potentialhits_forward[0], shell=True)  
                subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
                fwd=potentialhits_forward[0]+'.gz'
                rev=potentialhits_reverse[0]+'.gz'
            else:
                foldername=glob.glob(dr + '/' + filename + '*')
                if foldername and os.path.isdir(foldername[0]):
                    foldername=foldername[0]
                    potentialhits_forward=glob.glob(foldername + '/*' + filename + '*1*.fastq.gz')
                    potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq.gz')
                    if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                        fwd=potentialhits_forward[0]
                        rev=potentialhits_reverse[0]
                    elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
                        print(foldername + '/*' + filename + '*2*.fastq.gz')
                        potentialhits_forward=glob.glob(foldername +  '/*' + filename + '*1*.fastq')
                        potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq')
                        if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                            subprocess.run("gzip " + potentialhits_forward[0], shell=True)  
                            subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
                            fwd=potentialhits_forward[0]+'.gz'
                            rev=potentialhits_reverse[0]+'.gz'
    if not(fwd) or not(rev):
        raise ValueError('Either no file or more than 1 file found in ' + dr + 'for ' + ID)
    ##zip fastq files if they aren't already zipped
    subprocess.run("gzip " + fwd, shell=True)   
    subprocess.run("gzip " + rev, shell=True)   
    return [fwd, rev]
    
#Jonathan new code
# def findfastqfile(dr,ID,filename):
#     fwd=[]
#     rev=[]
#     potentialhits_forward=glob.glob(dr + '/' + filename +'/*1.fastq.gz')
#     potentialhits_reverse=glob.glob(dr + '/' + filename +'/*2.fastq.gz')
#     if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
#         fwd=potentialhits_forward[0]
#         rev=potentialhits_reverse[0]
#     elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
#         potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq.gz')
#         potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq.gz')
#         if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
#             fwd=potentialhits_forward[0]
#             rev=potentialhits_reverse[0]
#         elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
#             potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq')
#             potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq')
#             if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
#                 subprocess.run("gzip " + potentialhits_forward[0], shell=True)
#     # wildcards: sampleID=20221
#     # resources: mem_mb=1000, disk_mb=1000, tmpdir=/tmp
#                 subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
#                 fwd=potentialhits_forward[0]+'.gz'
#                 rev=potentialhits_reverse[0]+'.gz'
#             else:
#                 foldername=glob.glob(dr + '/' + filename + '*')
#                 if foldername and os.path.isdir(foldername[0]):
#                     foldername=foldername[0]
#                     potentialhits_forward=glob.glob(foldername + '/*' + filename + '*1*.fastq.gz')
#                     potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq.gz')
#                     if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
#                         fwd=potentialhits_forward[0]
#                         rev=potentialhits_reverse[0]
#                     elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
#                         print(foldername + '/*' + filename + '*2*.fastq.gz')
#                         potentialhits_forward=glob.glob(foldername +  '/*' + filename + '*1*.fastq')
#                         potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq')
#                         if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
#                             subprocess.run("gzip " + potentialhits_forward[0], shell=True)
#                             subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
#                             fwd=potentialhits_forward[0]+'.gz'
#                             rev=potentialhits_reverse[0]+'.gz'
#     if not(fwd) or not(rev):
#         raise ValueError('Either no file or more than 1 file found in ' + dr + 'for ' + ID)
#     ##make read only --- addition 20220415 by Jonathan (jdg)
#     subprocess.run("chmod 444 " + fwd, shell=True)
#     subprocess.run("chmod 444 " + rev, shell=True)
#     ##zip fastq files if they aren't already zipped
#     subprocess.run("gzip " + fwd, shell=True)
#     subprocess.run("gzip " + rev, shell=True)
#     return [fwd, rev]

def makelink(path,sample,filename,output_dir):
    #When sample is run on a single lane
    #File name can be either a COMPLETE directory name or a file name in batch(called path in this fx)
    [fwd_file, rev_file]=findfastqfile(path,sample, filename)
    # subprocess.run('ln -s -T ' + fwd_file + ' data/' + sample + '/R1.fq.gz', shell=Trueq)    
    # subprocess.run('ln -s -T ' + rev_file + ' data/' + sample + '/R2.fq.gz', shell=True)    
    
    subprocess.run(f"ln -s -T {fwd_file} {output_dir}/{sample}/R1.fq.gz", shell=True)    
    subprocess.run(f"ln -s -T {rev_file} {output_dir}/{sample}/R2.fq.gz", shell=True)    


def cp_append_files(paths,sample,filename,output_dir):
    #When sample is run on multiple lanes with same barcode
    fwd_list=''
    rev_list=''
    for path in paths:
        #Provider name can be either a COMPLETE directory name or a file name in batch(called path in this fx)
        [fwd_file, rev_file]=findfastqfile(path,sample, filename)
        fwd_list=fwd_list+ ' ' + fwd_file
        rev_list=rev_list+ ' ' + rev_file
        print(rev_list)
        print(fwd_list)
    subprocess.run("zcat " + fwd_list + ' | gzip > ' + output_dir + '/' +  sample + '/R1.fq.gz', shell=True)
    subprocess.run("zcat " + rev_list + ' | gzip > ' + output_dir + '/' +  sample + '/R2.fq.gz', shell=True)
    

def read_fasta(REFGENOME_DIR): 
    '''Reads in fasta file. If directory is given, reads in dir/genome.fasta
    Args:
        REFGENOME_DIR (str): Path to reference genome.

    Returns: SeqIO object for reference genome.
    '''
    fasta_file = glob.glob(REFGENOME_DIR + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(REFGENOME_DIR + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOME_DIR)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    
    return refgenome

def genomestats(REFGENOMEFOLDER):
    '''Parse genome to extract relevant stats

    Args:
        REFGENOMEFOLDER (str): Directory containing reference genome file.

    Returns:
        ChrStarts (arr): DESCRIPTION.
        Genomelength (arr): DESCRIPTION.
        ScafNames (arr): DESCRIPTION.

    '''

    refgenome = read_fasta(REFGENOMEFOLDER)
    
    Genomelength = 0
    ChrStarts = []
    ScafNames = []
    for record in refgenome:
        ChrStarts.append(Genomelength) # chr1 starts at 0 in analysis.m
        Genomelength = Genomelength + len(record)
        ScafNames.append(record.id)
    # close file
    #refgenome.close() # biopy update SeqIO has no close attribute anymore.
    # turn to np.arrys!
    ChrStarts = np.asarray(ChrStarts,dtype=int)
    Genomelength = np.asarray(Genomelength,dtype=int)
    ScafNames = np.asarray(ScafNames,dtype=object)
    return ChrStarts,Genomelength,ScafNames

def p2chrpos(p, ChrStarts):
    '''Convert 1col list of pos to 2col array with chromosome and pos on chromosome

    Args:
        p (TYPE): DESCRIPTION.
        ChrStarts (TYPE): DESCRIPTION.

    Returns:
        chrpos (TYPE): DESCRIPTION.

    '''
        
    # get chr and pos-on-chr
    chromo = np.ones(len(p),dtype=int)
    if len(ChrStarts) > 1:
        for i in ChrStarts[1:]:
            chromo = chromo + (p > i) # when (p > i) evaluates 'true' lead to plus 1 in summation. > bcs ChrStarts start with 0...genomestats()
        positions = p - ChrStarts[chromo-1] # [chr-1] -1 due to 0based index
        chrpos = np.column_stack((chromo,positions))
    else:
        chrpos = np.column_stack((chromo,p))
    return chrpos


# def get_clade_wildcards(cladeID):
#     is_clade = [int(i == cladeID) for i in GROUP_ls]
#     sampleID_clade = list(compress(SAMPLE_ls,is_clade))
#     reference_clade = list(compress(REF_Genome_ls,is_clade))
#     outgroup_clade = list(compress(OUTGROUP_ls,is_clade))
#     return sampleID_clade,reference_clade,outgroup_clade
    
# def get_sampleID_names(wildcards):  
#     sampleID_clade,_,_ = get_clade_wildcards(wildcards.cladeID)
#     return sampleID_clade

# def get_outgroup_bool(wildcards):  
#     _,_,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     return outgroup_clade

# def get_positions_prep(wildcards):
#     sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     mat_positions_prep=expand("Case/temp/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.pickle",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
#     return mat_positions_prep

# def get_diversity(wildcards):
#     sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     diversity_mat = expand("Mapping/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.diversity.pickle.gz",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
#     return diversity_mat   

# def get_quals(wildcards):
#     sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     quals_mat = expand("Mapping/quals/{sampleID}_ref_{reference}_outgroup{outgroup}.quals.pickle.gz",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
#     return quals_mat 

# def get_ref_genome(wildcards):
#     sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     return reference_clade
