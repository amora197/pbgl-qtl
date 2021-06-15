import os
import matplotlib
import yaml
import pandas as pd
import numpy as np
import seaborn as sns
import allel
import re

from pathlib import Path
from natsort import natsorted
from collections import OrderedDict
from matplotlib import pyplot as plt

sns.set()

pd.options.display.max_columns = None # default is 20
pd.options.display.max_rows = 60 # default is 60


###########################################################################################################
################################################ Functions ################################################
###########################################################################################################

# This function generates a tidy/longform table, which makes is easier to plot with seaborn.
# Will have one row per sample per variant, sample name in "sample", and rows for control/mutant observations and 
# frequency of mutant alleles
def getTidyTable(raw_table, samples, chromosome_lengths):
    
    # collects new rows
    row_accumulator = []

    def splitListToRows(row):
        '''Split one row into long form, one new row per sample'''
        if row['CHROM'] not in chromosome_lengths:
            # Only keep markers on the included chromosomes 
            return
        
        split_gt = row['GT'].split(',')
        if '.' in split_gt:
            # Remove rows where not all samples are genotyped
            return
        split_gq = row['GQ'].split(',')
        split_ro = row['SampleRO'].split(',')
        split_ao = row['SampleAO'].split(',')

        control_ref = split_gt[0] == '0/0'
        new_rows = []
        
        for i in range(4):
            new_row = row.to_dict()
            
            # Remove rows unnecessary for plotting
            new_row.pop('REF')
            new_row.pop('ALT')
            new_row.pop('RO')
            new_row.pop('AO')
            
            new_row['GT'] = split_gt[i]
            new_row['GQ'] = int(split_gq[i])
            new_row['SampleRO'] = split_ro[i]
            new_row['SampleAO'] = split_ao[i]
            new_row['sample'] = samples[i]
            
            # Get observations of WT/mutant alleles
            if control_ref:
                new_row['controlO'] = int(new_row['SampleRO'])
                new_row['mutantO'] = int(new_row['SampleAO'])
            else:
                new_row['controlO'] = int(new_row['SampleAO'])
                new_row['mutantO'] = int(new_row['SampleRO'])

            if (new_row['mutantO'] + new_row['controlO']) == 0:
                # Remove rows where there are no observations for a sample
                break
            else:
                # Get mutant allele frequency
                new_row['mutant_freq'] = new_row['mutantO'] / (new_row['mutantO'] + new_row['controlO'])
            new_rows.append(new_row)

        if len(new_rows) == 4:
            # only keep rows if there is one for each sample
            row_accumulator.extend(new_rows)

    raw_table.apply(splitListToRows, axis=1)
    table = pd.DataFrame(row_accumulator)
    table = table[['CHROM', 'POS', 'sample', 'GT', 'GQ', 
                   'SampleRO', 'SampleAO', 'controlO', 'mutantO', 'mutant_freq']]
    return table



# Create a table with averaged mutant frequencies for the F2 pools with a sliding window
def get_window_table(table, sample_tt, sample_TT, windowsize, stepsize, chromosome_lengths):
    
    # positions to include in window before/after the current pos, i.e. half the window size
    w = int(windowsize/2)
    
    rows = []
    for chrom in chromosome_lengths:
        
        # get F2 pool genotypes for current chromosomes
        chrom_table = table[(table['CHROM'] == chrom) 
                               & ((table['sample'] == sample_TT) 
                                  | (table['sample'] == sample_tt))].reset_index(drop=True)

        for i in range(1, chromosome_lengths[chrom]+1, stepsize):
            # Loop over windows
            
            wstart = max(1, i-w)
            wend = min(chromosome_lengths[chrom], i+w)
            
            # Get table for window
            wtable = chrom_table[(chrom_table['POS'] >= wstart) 
                                 & (chrom_table['POS'] <= wend)] 

            freqs_tt = wtable[wtable['sample'] == sample_tt]['mutant_freq']
            freqs_TT = wtable[wtable['sample'] == sample_TT]['mutant_freq']
            varcount = int(wtable.shape[0] / 2) # because of two samples per marker
            avg_TT = np.mean(freqs_TT)
            avg_tt = np.mean(freqs_tt)
            std_TT = np.std(freqs_TT)
            std_tt = np.std(freqs_tt)

            wsize = wend - wstart
            row_tt = {'CHROM': chrom,
                      'POS': i,
                      'sample': sample_tt,
                      'avg_mutant_freq': avg_tt,
                      'std': std_tt,
                      'varcount': varcount,
                      'window_size': wsize}
            rows.append(row_tt)
            row_TT = {'CHROM': chrom,
                      'POS': i,
                      'sample': sample_TT,
                      'avg_mutant_freq': avg_TT,
                      'std': std_TT,
                      'varcount': varcount,
                      'window_size': wsize}
            rows.append(row_TT)
    window_table = pd.DataFrame(rows)
    window_table = window_table[['CHROM', 'POS', 'sample', 'avg_mutant_freq', 'std', 'varcount', 'window_size']]
    return window_table

# Chromosome length extraction from VCF file
def chrom_length_extraction(vcf_file, chromosome_list):
    # naturally (human) sort chromosomes 
    chromosomes = natsorted(chromosome_list)
    
    ########## Headers Extraction ##########  
    # extract VCF headers
    # headers contains a list of lists/dictionaries
    # headers[0] - list of all the header lines
    # headers[1] - dict of FILTER field
    # headers[2] - dict of INFO fields
    # headers[3] - dict of FORMAT fields
    # headers[4] - list of samples
    headers = allel.read_vcf_headers(vcf_file)
    
    # extract chromosome lengths from 'headers[0]' by using 'contigs' keyword
    contigs = []
    for line in headers[0]:
        if 'contig' in line:
            for chromosome in chromosomes:
                if chromosome in line:
                    contigs.append(line)
                    break
                    
    # create dictionary of chromosome's name and length
    chrom_lengths = {}
    keyword_id = 'ID=(.*?),'
    keyword_len = 'length=(.+?)>'
    for contig in contigs:
        # search for keyword_id inside contig line
        ID = re.search(keyword_id, contig)
        for chromosome in chromosomes:
            if chromosome == ID.group(1):
                # search for keyword_len inside contig line
                length = re.search(keyword_len, contig)
                chrom_lengths[chromosome] = int(length.group(1))

    return OrderedDict(chrom_lengths)

# plot raw allele frequencies
def plot_allele_frequencies_raw(config_file):
    # open and read config file
    with open(config_file) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
        
    # extract VCF file 
    vcf_file = config['vcf_file']
    vcf_file_path = vcf_file['path']
    vcf_file_name = vcf_file['name']
    vcf_file_ext = vcf_file['extension']
    
    # step and window sizes
    window_size = config['window_size']
    step_size = config['step_size']
    
    # chromosomes extraction
    chromosomes = config['chromosomes']
    vcf_file_full_path = "%s/%s.%s" % (vcf_file_path, vcf_file_name, vcf_file_ext)
    chromosome_lengths = chrom_length_extraction(vcf_file_full_path, chromosomes)
    
    # samples ordered as: control, mutant, dwarf/mutant F2 pool, WT F2 pool
    control = config['samples']['control']
    mutant = config['samples']['mutant']
    sample_F2_WT = config['samples']['F2_wild_type']
    sample_F2_mutant = config['samples']['F2_mutant']
    samples = [control, mutant, sample_F2_mutant, sample_F2_WT]
    
    # TSV file creation with bash commands
    cmd_create_tsv = 'printf "#Samples:%s,%s,%s,%s\nCHROM\\tPOS\\tREF\\tALT\\tRO\\tAO\\tGT\\tGQ\\tSampleRO\\tSampleAO\n" > "Allele_Frequency_Plots_Computomics/%s.tsv"' % (
        samples[0], samples[1], samples[2], samples[3], vcf_file_name)
    os.system(cmd_create_tsv)
    
    data_flag = r'-f "%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\t[,%GT]\t[,%GQ]\t[,%RO]\t[,%AO]\n"'
    cmd_bcftools = 'bcftools view %s/%s.%s -s %s,%s,%s,%s | bcftools query %s >> Allele_Frequency_Plots_Computomics/%s.tsv' % (
        vcf_file_path, vcf_file_name, vcf_file_ext,
        samples[0], samples[1], samples[2], samples[3],
        data_flag, vcf_file_name)
    os.system(cmd_bcftools)
    
    # TSV file paths
    tsv_file = "Allele_Frequency_Plots_Computomics/%s.tsv" % vcf_file_name
    raw_output = "Allele_Frequency_Plots_Computomics/%s.pdf" % vcf_file_name
    raw_output_window = "Allele_Frequency_Plots_Computomics/%s_window.pdf" % vcf_file_name
    
    # Define input table paths
    tsv_file = Path(tsv_file)
    raw_output_pdf = Path(raw_output)
    window_output_pdf = Path(raw_output_window)
    
    # Load TSV table
    raw_table = pd.read_csv(tsv_file, sep='\t', na_values=['.'], comment='#')
    
    # Remove leading commas
    for col in ['GT', 'GQ', 'SampleRO', 'SampleAO']:
        raw_table.loc[:,col] = raw_table[col].str[1:]
        
    table = getTidyTable(raw_table, samples, chromosome_lengths)
    
    # Plot raw mutant allele frequencies of F2 pools
    plot_table = table[(table['sample'] == sample_F2_mutant) | (table['sample'] == sample_F2_WT)]
    plot = sns.relplot(data=plot_table, x='POS', y='mutant_freq', row='CHROM', style='sample',hue='sample', aspect=7.0, height=4., kind='line', markers=True, dashes=False, linewidth=0.5)
    plot.savefig(raw_output_pdf)

    
# plot weighted allele frequencies
def plot_allele_frequencies_weighted(config_file):    
    # open and read config file
    with open(config_file) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
        
    # extract VCF file 
    vcf_file = config['vcf_file']
    vcf_file_path = vcf_file['path']
    vcf_file_name = vcf_file['name']
    vcf_file_ext = vcf_file['extension']
    
    # step and window sizes
    window_size = config['window_size']
    step_size = config['step_size']
    
    # chromosomes extraction
    chromosomes = config['chromosomes']
    vcf_file_full_path = "%s/%s.%s" % (vcf_file_path, vcf_file_name, vcf_file_ext)
    chromosome_lengths = chrom_length_extraction(vcf_file_full_path, chromosomes)
    
    # samples ordered as: control, mutant, dwarf/mutant F2 pool, WT F2 pool
    control = config['samples']['control']
    mutant = config['samples']['mutant']
    sample_F2_WT = config['samples']['F2_wild_type']
    sample_F2_mutant = config['samples']['F2_mutant']
    samples = [control, mutant, sample_F2_mutant, sample_F2_WT]
    
    # TSV file creation with bash commands
    cmd_create_tsv = 'printf "#Samples:%s,%s,%s,%s\nCHROM\\tPOS\\tREF\\tALT\\tRO\\tAO\\tGT\\tGQ\\tSampleRO\\tSampleAO\n" > "Allele_Frequency_Plots_Computomics/%s.tsv"' % (
        samples[0], samples[1], samples[2], samples[3], vcf_file_name)
    os.system(cmd_create_tsv)
    
    data_flag = r'-f "%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\t[,%GT]\t[,%GQ]\t[,%RO]\t[,%AO]\n"'
    cmd_bcftools = 'bcftools view %s/%s.%s -s %s,%s,%s,%s | bcftools query %s >> Allele_Frequency_Plots_Computomics/%s.tsv' % (
        vcf_file_path, vcf_file_name, vcf_file_ext,
        samples[0], samples[1], samples[2], samples[3],
        data_flag, vcf_file_name)
    os.system(cmd_bcftools)
    
    # TSV file paths
    tsv_file = "Allele_Frequency_Plots_Computomics/%s.tsv" % vcf_file_name
    raw_output = "Allele_Frequency_Plots_Computomics/%s.pdf" % vcf_file_name
    raw_output_window = "Allele_Frequency_Plots_Computomics/%s_window.pdf" % vcf_file_name
    
    # Define input table paths
    tsv_file = Path(tsv_file)
    raw_output_pdf = Path(raw_output)
    window_output_pdf = Path(raw_output_window)
    
    # Load TSV table
    raw_table = pd.read_csv(tsv_file, sep='\t', na_values=['.'], comment='#')
    
    # Remove leading commas
    for col in ['GT', 'GQ', 'SampleRO', 'SampleAO']:
        raw_table.loc[:,col] = raw_table[col].str[1:]
        
    table = getTidyTable(raw_table, samples, chromosome_lengths)
        
    # Plot the per-window average allele frequencies
    # Create one subplot per chromosome
    wtable = get_window_table(table, sample_F2_mutant, sample_F2_WT, window_size, step_size, chromosome_lengths)
    wgrid = sns.FacetGrid(wtable, row='CHROM', aspect=7.0, height=4., hue='sample', legend_out=True)

    # Plot lines
    wgrid = wgrid.map(sns.lineplot, 'POS', 'avg_mutant_freq', linewidth=0.5).add_legend()

    # Overlay with scatterplot with densities per chromosome
    for i, chrom in enumerate(chromosome_lengths):
        plotdata = wtable[wtable['CHROM']==chrom]
        mincount = min(plotdata['varcount'])
        maxcount = max(plotdata['varcount'])
        sizes=(np.log2(mincount+1)*100,np.log2(maxcount+1)*100)
        sns.scatterplot(data=wtable[wtable['CHROM']==chrom], 
                        x='POS', y='avg_mutant_freq', size='varcount', hue='sample', legend=False,
                        ax=wgrid.axes[i,0], sizes=sizes)
    # save figure
    wgrid.savefig(window_output_pdf)