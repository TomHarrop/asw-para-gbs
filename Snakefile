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


def generate_full_popmap(key_file, full_popmap):
    key_data = pandas.read_csv(key_file, delimiter='\t')
    subset = key_data['sample']
    individual = list(re.sub('\s+', '_', x) for x in subset)
    population = list(re.sub("[\d|_]+", "", x) for x in individual)
    popmap = pandas.DataFrame({'individual': individual,
                               'population': population})
    popmap.to_csv(full_popmap,
                  sep='\t',
                  header=False,
                  index=False)


def filter_combined_covstats(cutoff, covstats_file, filtered_popmap):
    covstats = pandas.read_csv(covstats_file)
    passed_filter = covstats[covstats['final_coverage_mean'] > cutoff]
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
        'output/04_stacks/gstacks.vcf.gz',
        'output/03_params/compare_defaults/optimised_samplestats_combined.csv',
        expand(('output/11_stacks-populations/r{r}/'
                'populations.sumstats_summary.tsv'),
               r=[0.8])

# 11 generate filtered SNP table for adegenet
rule populations:
    input:
        'output/04_stacks/catalog.fa.gz',
        'output/04_stacks/catalog.calls',
        map = 'output/01_config/filtered_popmap.tsv'
    output:
        'output/11_stacks-populations/r{r}/populations.sumstats_summary.tsv',
        'output/11_stacks-populations/r{r}/populations.markers.tsv',
        'output/11_stacks-populations/r{r}/populations.hapstats.tsv',
        'output/11_stacks-populations/r{r}/populations.sumstats.tsv',
        'output/11_stacks-populations/r{r}/populations.haplotypes.tsv',
        'output/11_stacks-populations/r{r}/populations.snps.genepop',
        'output/11_stacks-populations/r{r}/populations.snps.vcf'
    params:
        stacks_dir = 'output/04_stacks',
        outdir = 'output/11_stacks-populations/r{r}'
    threads:
        50
    log:
        'output/logs/11_stacks-populations/r{r}.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-O {params.outdir} '
        '-t {threads} '
        '-r {wildcards.r} '
        '--genepop --vcf '
        '&> {log}'

rule link_catalog_output2:
    input:
        'output/04_stacks/gstacks.vcf.gz'
    output:
        'output/04_stacks/catalog.calls'
    shell:
        'ln -s '
        '"$(readlink -f {input})" '
        '"$(readlink -f {output})"'

rule link_catalog_output:
    input:
        'output/04_stacks/gstacks.fa.gz'
    output:
        'output/04_stacks/catalog.fa.gz'
    shell:
        'ln -s '
        '"$(readlink -f {input})" '
        '"$(readlink -f {output})"'


# 10 integrate the alignment information
rule integrate_alignments:
    input:
        bam = 'output/09_bwa/gstacks_filtered.bam',
        fa = 'output/04_stacks/gstacks.fa.gz',
        vcf = 'output/04_stacks/gstacks.vcf.gz'
    output:
        stacks_indir = temp('output/10_stacks'),
        tmp_fa = temp('output/10_stacks/catalog.fa.gz'),
        tmp_vcf = temp('output/10_stacks/catalog.calls'),
        fa = 'output/10_aligned-catalog/gstacks.fa.gz',
        vcf = 'output/10_aligned-catalog/gstacks.vcf.gz'
    params:
        stacks_outdir = 'output/10_aligned-catalog'
    log:
        'output/logs/10_integrate-alignments.log'
    shell:
        'cp {input.fa} {output.tmp_fa} & '
        'cp {input.vcf} {output.tmp_vcf} & '
        'wait ; '
        'stacks-integrate-alignments '
        '-P {output.stacks_indir} '
        '-B {input.bam} '
        '-O {params.stacks_outdir} '
        '&> {log}'


# 09 map the catalog to the reference genome
rule sam_to_bam:
    input:
        sam = 'output/09_bwa/gstacks.sam'
    output:
        bam = 'output/09_bwa/gstacks_filtered.bam'
    threads:
        75
    log:
        'output/log/09_sam-to-bam.log'
    shell:
        'samtools view '
        '-S '
        # '-f 16 '
        '-b -h '
        '-@ {threads} '
        '-o {output.bam} '
        '{input.sam} '
        '2> {log}'


rule bwa_mem:
    input:
        ref = 'output/09_bwa/mh.fa',
        bwt = 'output/09_bwa/mh.fa.bwt',
        fa = 'output/04_stacks/gstacks.fa.gz'
    output:
        fa = temp('output/09_bwa/gstacks.fa'),
        sam = 'output/09_bwa/gstacks.sam'
    threads:
        75
    log:
        'output/log/09_bwa-mem.log'
    shell:
        'zcat {input.fa} > {output.fa} ; '
        'bwa mem '
        '-t {threads} '
        '-S -P '
        '{input.ref} '
        '{output.fa} '
        '> {output.sam} '
        '2> {log}'

rule bwa_index:
    input:
        fa = 'data/mh/final.scaffolds.fa'
    output:
        fa = 'output/09_bwa/mh.fa',
        bwt = 'output/09_bwa/mh.fa.bwt',
        sa = 'output/09_bwa/mh.fa.sa'
    threads:
        1
    log:
        'output/logs/09_bwa-index.log'
    shell:
        'cp {input.fa} {output.fa} ; '
        'bwa index {output.fa} &> {log}'    


# 08 generate catalog
rule gstacks:
    input:
        dynamic('output/04_stacks/{tsv2bam_indiv}.bam'),
        filtered_popmap = 'output/01_config/filtered_popmap.tsv'
    output:
        'output/04_stacks/gstacks.fa.gz',
        'output/04_stacks/gstacks.vcf.gz'
    params:
        stacks_dir = 'output/04_stacks'
    threads:
        75
    log:
        'output/logs/08_gstacks.log'
    shell:
        'gstacks '
        '-P {params.stacks_dir} '
        '-M {input.filtered_popmap} '
        '-t {threads} '
        '&> {log}'

# 07 match samples to the catalog and convert to bam
rule tsv2bam:
    input:
        dynamic('output/04_stacks/{sstacks_indiv}.matches.tsv.gz'),
        filtered_popmap = 'output/01_config/filtered_popmap.tsv'
    output:
        dynamic('output/04_stacks/{tsv2bam_indiv}.bam')
    params:
        stacks_dir = 'output/04_stacks'
    threads:
        75
    log:
        'output/logs/07_sstacks/tsv2bam.log'
    shell:
        'tsv2bam '
        '-P {params.stacks_dir} '
        '-M {input.filtered_popmap} '
        '-t {threads} '
        '&> {log} '

rule sstacks:
    input:
        catalog = 'output/04_stacks/batch_1.catalog.tags.tsv.gz',
        filtered_popmap = 'output/01_config/filtered_popmap.tsv'
    output:
        dynamic('output/04_stacks/{sstacks_indiv}.matches.tsv.gz')
    params:
        stacks_dir = 'output/04_stacks'
    threads:
        75
    log:
        'output/logs/07_sstacks/sstacks.log'
    shell:
        'sstacks '
        '-P {params.stacks_dir} '
        '-M {input.filtered_popmap} '
        '-p {threads} '
        '&> {log}'

# 06 generate the catalog
rule cstacks:
    input:
        dynamic('output/04_stacks/{individual}.alleles.tsv.gz'),
        map = 'output/01_config/filtered_popmap.tsv'
    output:
        'output/04_stacks/batch_1.catalog.tags.tsv.gz',
        'output/04_stacks/batch_1.catalog.snps.tsv.gz',
        'output/04_stacks/batch_1.catalog.alleles.tsv.gz'
    params:
        stacks_dir = 'output/04_stacks'
    threads:
        75
    log:
        'output/logs/06_cstacks.log'
    shell:
        'cstacks '
        '-p {threads} '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-n 3 '
        '&> {log} '

# 05 calculate coverage stats per individual and filter
rule filter_by_coverage:
    input:
        combined_covstats = ('output/05_covstats/'
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
        dynamic('output/05_covstats/{individual}.csv')
    output:
        combined = 'output/05_covstats/individual_covstats_combined.csv'
    script:
        'src/combine_csvs.R'

rule individual_covstats:
    input:
        tags_file = 'output/04_stacks/{individual}.tags.tsv.gz'
    output:
        covstats = 'output/05_covstats/{individual}.csv'
    log:
        log = 'output/logs/05_covstats/{individual}.log'
    threads:
        1
    script:
        'src/calculate_mean_coverage.R'

# 04 assemble loci
rule ustacks:
    input:
        fastq = 'output/02_demux/{individual}.fq.gz',
        individual_i_pickle = 'output/01_config/individual_i.p'
    params:
        wd = 'output/04_stacks'
    output:
        'output/04_stacks/{individual}.alleles.tsv.gz',
        'output/04_stacks/{individual}.snps.tsv.gz',
        'output/04_stacks/{individual}.models.tsv.gz',
        'output/04_stacks/{individual}.tags.tsv.gz'
    threads:
        10
    log:
        'output/logs/04_ustacks/{individual}.log'
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
              '-m 5 '
              '-M 2 '
              '&> {log}')   


# 03 parameter optimisation
rule optim_compare:
    input:
        dynamic('output/02_demux/{individual}.fq.gz'),
        full_popmap = 'output/01_config/full_popmap.txt',
        mm_results = 'output/03_params/stats_Mm/samplestats_combined.csv',
        n_results = 'output/03_params/stats_n/samplestats_combined.csv'
    output:
        'output/03_params/compare_defaults/optimised_samplestats_combined.csv'
    params:
        sample_dir = 'output/02_demux',
        outdir = 'output/03_params'
    threads:
        50
    log:
        'output/logs/03_params/compare.log'
    shell:
        'stacks_parameters '
        '--mode compare_defaults -M 2 -m 5 -n 2 '
        '-o {params.outdir} '
        '--threads {threads} '
        '{input.full_popmap} '
        '{params.sample_dir} '
        '&> {log}'


rule optim_n:
    input:
        dynamic('output/02_demux/{individual}.fq.gz'),
        full_popmap = 'output/01_config/full_popmap.txt',
        mm_results = 'output/03_params/stats_Mm/samplestats_combined.csv'
    output:
        'output/03_params/stats_n/samplestats_combined.csv'
    params:
        sample_dir = 'output/02_demux',
        outdir = 'output/03_params'
    threads:
        50
    log:
        'output/logs/03_params/optim-n.log'
    shell:
        'stacks_parameters '
        '--mode optim_n -M 2 -m 5 '
        '-o {params.outdir} '
        '--threads {threads} '
        '{input.full_popmap} '
        '{params.sample_dir} '
        '&> {log}'

rule optim_mm:
    input:
        dynamic('output/02_demux/{individual}.fq.gz'),
        full_popmap = 'output/01_config/full_popmap.txt',
        setup_popmap = 'output/03_params/filtering/replicate_1_popmap.txt'
    output:
        'output/03_params/stats_Mm/samplestats_combined.csv'
    params:
        sample_dir = 'output/02_demux',
        outdir = 'output/03_params'
    threads:
        50
    log:
        'output/logs/03_params/optim-mm.log'
    shell:
        'stacks_parameters '
        '--mode optim_Mm '
        '-o {params.outdir} '
        '--threads {threads} '
        '{input.full_popmap} '
        '{params.sample_dir} '
        '&> {log}'

rule optim_setup:
    input:
        dynamic('output/02_demux/{individual}.fq.gz'),
        full_popmap = 'output/01_config/full_popmap.txt'
    output:
        'output/03_params/filtering/replicate_1_popmap.txt'
    params:
        sample_dir = 'output/02_demux',
        outdir = 'output/03_params'
    threads:
        50
    log:
        'output/logs/03_params_setup.log'
    shell:
        'stacks_parameters '
        '--mode setup '
        '-o {params.outdir} '
        '--threads {threads} '
        '{input.full_popmap} '
        '{params.sample_dir} '
        '&> {log}'

rule full_popmap:
    input:
        key_file = key_file
    output:
        full_popmap = 'output/01_config/full_popmap.txt'
    run:
        generate_full_popmap(input.key_file, output.full_popmap)


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
