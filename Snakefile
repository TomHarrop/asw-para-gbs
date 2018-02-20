#!/usr/bin/env python3

import os
import pandas
import pickle
import re


#############
# FUNCTIONS #
#############

def read_key_and_write_config(key_file, outdir):
    key_data = pandas.read_csv(key_file, delimiter='\t')
    grouped_key_data = key_data.groupby(['flowcell', 'lane'])
    for name, group in grouped_key_data:
        prefix = '_'.join([str(x) for x in name])
        config_file = os.path.join(outdir, '%s.config' % prefix)
        subset = group[['barcode', 'sample']]
        subset['sample'] = list(re.sub('\s+', '_', x)
                                for x in subset['sample'])
        if len(subset) > 0:
            subset.to_csv(config_file,
                          sep='\t',
                          header=False,
                          index=False)


def filter_combined_covstats(cutoff, covstats_file, filtered_popmap):
    covstats = pandas.read_csv(covstats_file)
    passed_filter = covstats[covstats['primary'] > cutoff]
    passed_filter['population'] = passed_filter['individual'].transform(
        lambda x: re.sub("[\d|_]+", "", x))
    popmap = passed_filter[['individual', 'population']]
    popmap.to_csv(filtered_popmap,
                  sep='\t',
                  header=False,
                  index=False)


def write_flag_files(filtered_popmap, outdir):
    passed_filter = pandas.read_csv(filtered_popmap, sep='\t')
    for index, row in passed_filter.iterrows():
        my_flag = os.path.join(outdir, row[0])
        with open(my_flag, 'wt') as f:
            f.write('')

###########
# GLOBALS #
###########

# files and folders
key_file = 'data/C6JRFANXX/SQ2532.txt'
read_file = 'data/C6JRFANXX/SQ2532_NoIndex_L007_R1_001.fastq.gz'

#########
# SETUP #
#########


#########
# RULES #
#########

rule target:
    input:
        'output/03_ustacks/batch_1.catalog.tags.tsv.gz'

# 05 generate the catalog (takes ~ 1 week)
rule cstacks:
    input:
        map = 'output/01_config/filtered_popmap.tsv',
        flag = dynamic('output/01_config/flags/{filtered_indiv}'),
    output:
        'output/03_ustacks/batch_1.catalog.tags.tsv.gz',
        'output/03_ustacks/batch_1.catalog.snps.tsv.gz',
        'output/03_ustacks/batch_1.catalog.alleles.tsv.gz',
    params:
        ustacks_dir = 'output/03_ustacks'
    threads:
        75
    log:
        'output/logs/05_cstacks.log'
    shell:
        'cstacks '
        '-p {threads} '
        '-P {params.ustacks_dir} '
        '-M {input.map} '
        '-n 3 '
        '&> {log} '
  
rule cstacks_flags:
    input:
        filtered_popmap = 'output/01_config/filtered_popmap.tsv'
    output:
        dynamic('output/01_config/flags/{filtered_indiv}')
    params:
        outdir = 'output/01_config/flags'
    run:
        write_flag_files(
            filtered_popmap=input.filtered_popmap,
            outdir=params.outdir)

# 04 calculate coverage stats per individual and filter
rule filter_by_coverage:
    input:
        combined_covstats = ('output/04_covstats/'
                             'individual_covstats_combined.csv')
    output:
        filtered_popmap = 'output/01_config/filtered_popmap.tsv'
    params:
        cutoff = 15
    run:
        filter_combined_covstats(
            cutoff=params.cutoff,
            covstats_file=input.combined_covstats,
            filtered_popmap=output.filtered_popmap)


rule combine_individual_covstats:
    input:
        dynamic('output/04_covstats/{individual}.csv')
    output:
        combined = 'output/04_covstats/individual_covstats_combined.csv'
    script:
        'src/combine_csvs.R'

rule individual_covstats:
    input:
        tags_file = 'output/03_ustacks/{individual}.tags.tsv.gz'
    output:
        covstats = 'output/04_covstats/{individual}.csv'
    log:
        log = 'output/logs/04_covstats/{individual}.log'
    threads:
        1
    script:
        'src/calculate_mean_coverage.R'

# 03 assemble loci
rule ustacks:
    input:
        fastq = 'output/02_demux/{individual}.fq.gz',
        individual_i_pickle = 'output/01_config/individual_i.p'
    params:
        wd = 'output/03_ustacks'
    output:
        'output/03_ustacks/{individual}.alleles.tsv.gz',
        'output/03_ustacks/{individual}.snps.tsv.gz',
        'output/03_ustacks/{individual}.models.tsv.gz',
        'output/03_ustacks/{individual}.tags.tsv.gz'
    threads:
        15
    log:
        'output/logs/03_ustacks/{individual}.log'
    run:
        # open the pickled dictionary and look up the sample_i
        with open(input.individual_i_pickle, 'rb') as f:
            individual_i = pickle.load(f)
        sample_i = individual_i[wildcards.individual]
        shell('ustacks '
              '-p {threads} '
              '-t gzfastq '
              '-f {input.fastq} '
              '-o {params.wd} '
              '-i {sample_i} '
              '-m 3 '
              '-M 3 '
              '&> {log}')

# 02b. make a dictionary of sample:i for cstacks
rule enumerate_filtered_samples:
    input:
        dynamic('output/02_demux/{individual}.fq.gz')
    output:
        pickle = 'output/01_config/individual_i.p'
    run:
        # read the filtered popmap
        my_files = [re.sub('\.fq\.gz$', '', os.path.basename(x))
                    for x in input]
        my_individuals = enumerate(sorted(set(my_files)))
        individual_i = {y: x for x, y in my_individuals}
        # pickle the individual_i dict for other rules to use
        with open(output.pickle, 'wb+') as f:
            pickle.dump(individual_i, f)

# 02a demux
rule process_radtags:
    input:
        read_file = read_file,
        config_file = 'output/01_config/C6JRFANXX_7.config'
    output:
        dynamic('output/02_demux/{individual}.fq.gz')
    params:
        outdir = 'output/02_demux'
    log:
        'output/logs/02_demux.log'
    threads:
        1
    shell:
        'process_radtags '
        '-f {input.read_file} '
        '-i gzfastq -y gzfastq '
        '-b {input.config_file} '
        '-o {params.outdir} '
        '-c -q '
        '-t 91 '                    # truncate output to 91 b
        '-w 0.1 '                   # window: approx. 9 bases
        '-s 15 '                    # minimum avg PHRED in window
        '--inline_null '
        '--renz_1 pstI '
        '&> {log}'

# 01 generate a key file for demuxing
rule generate_demux_key:
    input:
        key_file = key_file
    output:
        'output/01_config/C6JRFANXX_7.config'
    params:
        outdir = 'output/01_config'
    threads:
        1
    run:
        read_key_and_write_config(input.key_file, params.outdir)
