'''fasterq-dump for all runs'''
from pandas import read_table
from snakemake.utils import validate, min_version
from os.path import abspath, join
##### set minimum snakemake version #####
min_version("6.5.1")

##### load config and sample sheets #####
configfile: "config.yaml"

samples = read_table(config["samples"]).set_index(["run"], drop=False)
temp_dir = abspath("./temp_gz")
##### target rules #####

shell.executable("/bin/bash")
shell.prefix("export LC_ALL=en_US.utf-8 && export LANG=en_US.utf-8 &&")

rule all:
    input:
        expand([join(temp_dir, "{accession}_1.fastq.gz"),
                join(temp_dir, "{accession}_2.fastq.gz")],
                accession=samples["run"])


##### load rules #####

rule get_fastq:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        join(temp_dir, "{accession}_1.fastq"),
        join(temp_dir, "{accession}_2.fastq"),
    params:
        extra=" --split-files --include-technical --outdir " + temp_dir,
        run="{accession}",

    threads: 10  # defaults to 6
    resources:
        mem_mb=20000
    shell:
        ("fasterq-dump -p -x --threads {threads} --mem {resources.mem_mb}M {params.extra} {params.run}")

rule get_fastq_gz:
    input:
       join(temp_dir, "{accession}_1.fastq"),
       join(temp_dir, "{accession}_2.fastq")
    output:
       join(temp_dir, "{accession}_1.fastq.gz"),
       join(temp_dir, "{accession}_2.fastq.gz")
    threads: 12
    shell:
        ("pigz -p {threads} {input}")