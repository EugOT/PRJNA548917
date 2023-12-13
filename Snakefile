"""Romanov et al., 2020"""
import pandas as pd
from os import listdir, rename, getcwd
from os.path import join, basename, dirname, abspath
from pathlib import Path
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("7.20.0")

##### load config and sample sheets #####
configfile: "config.yaml"
samples = pd.read_table(config["samples"]).set_index("Run", drop=False)
resolutions = [0.001, 0.01, 0.05, 0.1]
bprj = "PRJNA548917"
prj = "romanov2020_Hypoth-dev"

def plots_doublets_raw(wildcards):
    x = "output/figures/{wildcards.run}_raw/doublets_call_FPR_{wildcards.res}".format(wildcards=wildcards)
    return x.replace("\.", "_")


def get_mem_mb(wildcards, attempt):
    return attempt * 500000


##### target rules #####

shell.executable("/bin/bash")

rule all:
    input:
        expand("cellbender/{run}/{run}_output_FPR_{res}_filtered.h5",
                run=samples["Run"], res=resolutions),
        expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("cellranger/{run}/outs/filtered_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("scrublet/{run}/{run}_initial_annotation_FPR_{res}.h5ad",
                run=samples["Run"], res=resolutions),
        expand("piscem_spliceu/{run}/af_quant/alevin/quants_mat.mtx",
                run=samples["Run"]),
        expand(["output/figures/combined-top5_logreg-umap-whole_dataset-fpr_{res}.pdf",
            "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_{res}.pdf",
            "output/tables/01A-eda-whole_dataset-fpr_{res}/parameters.json"], res=resolutions),
        expand(["data/{bprj}-whole_dataset-fpr_{res}-clusters.h5Seurat",
            "data/{bprj}-whole_dataset-fpr_{res}-clusters.h5ad",
            "data/class_cello/{bprj}-astrocytes_dataset-{res}-initial_selection.h5ad"], bprj=bprj, res=resolutions),
        expand(["output/tables/01A-eda-whole_dataset-fpr_{res}/{prj}_all_mrk-MAST_sct-combined-whole_dataset-fpr_{res}.csv",
            "output/tables/01A-eda-whole_dataset-fpr_{res}/{prj}_all_mrk-logreg_sct-combined-whole_dataset-fpr_{res}.csv"], prj=prj, res=resolutions),
        expand(["data/{bprj}-whole_dataset-nc-clusters.h5Seurat",
            "data/{bprj}-whole_dataset-nc-clusters.h5ad",
            "data/class_cello/{bprj}-astrocytes_dataset-nc-initial_selection.h5ad"], bprj=bprj),
        expand(["output/tables/01-eda-whole_dataset-nc/{prj}_all_mrk-MAST_sct-combined-whole_dataset-nc.csv",
            "output/tables/01-eda-whole_dataset-nc/{prj}_all_mrk-logreg_sct-combined-whole_dataset-nc.csv"], prj=prj),
        "output/figures/combined-top5_logreg-umap-whole_dataset-nc.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-nc.pdf",
        "output/tables/01-eda-whole_dataset-nc/parameters.json",
        expand(["output/tables/{bprj}-astrocytes_dataset-0.001-graph-regulons.graphml"], bprj=bprj),

##### load rules #####

CELLRANGER="source /home/etretiakov/src/cellranger-7.1.0/sourceme.bash && cellranger "

rule cellranger_count:
    output:
        raw="cellranger/{run}/outs/raw_feature_bc_matrix.h5",
        filtered="cellranger/{run}/outs/filtered_feature_bc_matrix.h5",
        summary="cellranger/{run}/outs/web_summary.html",
        bam="cellranger/{run}/outs/possorted_genome_bam.bam",
    params:
        ids="cellranger/{run}",
        sample="{run}",
        idx="../mm10_optimized"
    threads: 32
    resources:
        mem_mb=64000
    shell:
        ("{run} count --include-introns true \
            --id={params.ids} \
            --sample={params.sample} \
            --transcriptome={params.idx} \
            --fastqs=fastq \
            --jobmode=local \
            --localcores={threads} ")

rule cellbender:
    input:
        "cellranger/{run}/outs/raw_feature_bc_matrix.h5"
    output:
        expand(["cellbender/{{run}}/{{run}}_output_FPR_{res}.h5", "cellbender/{{run}}/{{run}}_output_FPR_{res}_filtered.h5"], res=resolutions)
    params:
        ndroplets=lambda wildcards: samples["NTotalDropletsIncluded"][wildcards.run],
        ncells=lambda wildcards: samples["NTotalCells"][wildcards.run],
        h5="cellbender/{run}/{run}_output.h5"
    container:
        "docker://etretiakov/cellbender:v0.0.1"
    threads: 4
    resources:
        nvidia_gpu=1,
        mem_mb=10000
    shell:
        ("cellbender remove-background \
            --input {input} \
            --output {params.h5} \
            --cuda \
            --expected-cells {params.ncells} \
            --total-droplets-included {params.ndroplets} \
            --fpr 0.001 0.01 0.05 0.1 \
            --epochs 150")

rule doublets_call:
    input:
        filt_h5="cellbender/{run}/{run}_output_FPR_{res}_filtered.h5"
    output:
        scrublet_calls="scrublet/{run}/{run}_scrublet_calls_FPR_{res}.tsv",
        dr="cellbender/{run}/{run}_latent_gene_expression_FPR_{res}.csv",
        h5ad="scrublet/{run}/{run}_initial_annotation_FPR_{res}.h5ad"
    params:
        expected_dblt=lambda wildcards: samples["NExpectedDoubletRate"][wildcards.run],
        sample_run_name="{run}",
        plots=plots_doublets_raw
    container:
        "docker://etretiakov/scrna-seq:jammy-2022.12.09-v0.0.1"
    threads: 8
    resources:
        mem_mb=20000
    script:
        "../scrublet_cb-z.py"


rule exploratory_data_analysis_0_001:
    input:
        rmd="analysis/01A-eda-whole_dataset-fpr_0.001.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        cellbender=expand("cellbender/{run}/{run}_output_FPR_0.001_filtered.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.001.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.001.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.001.pdf",
        f"output/tables/01A-eda-whole_dataset-fpr_0.001/{prj}_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.001.csv",
        f"output/tables/01A-eda-whole_dataset-fpr_0.001/{prj}_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.001.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.001/parameters.json",
        f"data/{bprj}-whole_dataset-fpr_0.001-clusters.h5Seurat"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 6
    resources:
        mem_mb=get_mem_mb
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


rule exploratory_data_analysis_0_01:
    input:
        rmd="analysis/01A-eda-whole_dataset-fpr_0.01.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        cellbender=expand("cellbender/{run}/{run}_output_FPR_0.01_filtered.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.01.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.01.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.01.pdf",
        f"output/tables/01A-eda-whole_dataset-fpr_0.01/{prj}_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.01.csv",
        f"output/tables/01A-eda-whole_dataset-fpr_0.01/{prj}_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.01.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.01/parameters.json",
        f"data/{bprj}-whole_dataset-fpr_0.01-clusters.h5Seurat"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 6
    resources:
        mem_mb=get_mem_mb
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


rule exploratory_data_analysis_0_05:
    input:
        rmd="analysis/01A-eda-whole_dataset-fpr_0.05.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        cellbender=expand("cellbender/{run}/{run}_output_FPR_0.05_filtered.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.05.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.05.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.05.pdf",
        f"output/tables/01A-eda-whole_dataset-fpr_0.05/{prj}_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.05.csv",
        f"output/tables/01A-eda-whole_dataset-fpr_0.05/{prj}_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.05.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.05/parameters.json",
        f"data/{bprj}-whole_dataset-fpr_0.05-clusters.h5Seurat"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 6
    resources:
        mem_mb=get_mem_mb
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


rule exploratory_data_analysis_0_1:
    input:
        rmd="analysis/01A-eda-whole_dataset-fpr_0.1.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        cellbender=expand("cellbender/{run}/{run}_output_FPR_0.1_filtered.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.1.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.1.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.1.pdf",
        f"output/tables/01A-eda-whole_dataset-fpr_0.1/{prj}_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.1.csv",
        f"output/tables/01A-eda-whole_dataset-fpr_0.1/{prj}_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.1.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.1/parameters.json",
        f"data/{bprj}-whole_dataset-fpr_0.1-clusters.h5Seurat"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 6
    resources:
        mem_mb=get_mem_mb
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


rule map_spliceu:
    input:
        r1="fastq/{run}_S1_L001_R1_001.fastq.gz",
        r2="fastq/{run}_S1_L001_R2_001.fastq.gz",
        barcodes="3M-february-2018.txt"
    output:
        map="piscem_spliceu/{run}/piscem_map/map.rad"
    params:
        prefix="piscem_spliceu/{run}/piscem_map",
        act="piscem_spliceu/{run}/.afhome"
    container:
        "docker://combinelab/usefulaf:0.9.3"
    threads: 32
    resources:
        mem_mb=32000
    benchmark:
        "benchmarks/piscem_spliceu_map/{run}.tsv"
    shell:
        ("export ALEVIN_FRY_HOME={params.act} \
        && simpleaf set-paths \
        && piscem map-sc \
        --index /data/GRCm39/index/piscem_idx \
        --threads {threads} \
        -o {params.prefix} \
        -1 {input.r1} \
        -2 {input.r2} \
        --geometry chromium_v3 ")


rule quant_spliceu:
    input:
        map="piscem_spliceu/{run}/piscem_map/map.rad"
    output:
        result="piscem_spliceu/{run}/af_quant/alevin/quants_mat.mtx"
    params:
        prefix="piscem_spliceu/{run}/",
        map="piscem_spliceu/{run}/piscem_map",
        act="piscem_spliceu/{run}/.afhome"
    container:
        "docker://combinelab/usefulaf:0.9.3"
    threads: 32
    resources:
        mem_mb=32000
    benchmark:
        "benchmarks/piscem_spliceu/{run}.tsv"
    shell:
        ("export ALEVIN_FRY_HOME={params.act} \
        && simpleaf set-paths \
        && simpleaf quant \
        -c 10xv3 \
        -o {params.prefix} \
        -t {threads} \
        --map-dir {params.map} \
        -r cr-like -u \
        -m /data/GRCm39/index/t2g_3col.tsv ")


rule get_h5ad:
    input:
        gene="piscem_spliceu/{run}/af_quant/alevin/quants_mat.mtx"
    output:
        knee="output/figures/{run}_raw/knee-plot.pdf",
        h5ad="piscem_spliceu/{run}/{run}.h5ad"
    params:
        sample_run_name="{run}",
        expected_num_cells=lambda wildcards: samples["NTotalCells"][wildcards.run],
        path="piscem_spliceu/{run}/af_quant"
    container:
        "docker://combinelab/usefulaf:0.9.3"
    threads: 4
    resources:
        mem_mb=16000
    script:
        "../code/pyroe.py"


rule cellbender_velocity:
    input:
        "piscem_spliceu/{run}/{run}.h5ad"
    output:
        expand(["cellbender_velocyto/{{run}}/{{run}}_output_FPR_{res}.h5", "cellbender_velocyto/{{run}}/{{run}}_output_FPR_{res}_filtered.h5"], res=resolutions)
    params:
        ndroplets=lambda wildcards: samples["NTotalDropletsIncluded"][wildcards.run],
        ncells=lambda wildcards: samples["NTotalCells"][wildcards.run],
        h5="cellbender_velocyto/{run}/{run}_output.h5"
    container:
        "docker://etretiakov/cellbender:v0.0.1"
    threads: 4
    resources:
        nvidia_gpu=1,
        mem_mb=10000
    shell:
        ("cellbender remove-background \
            --input {input} \
            --output {params.h5} \
            --cuda \
            --expected-cells {params.ncells} \
            --total-droplets-included {params.ndroplets} \
            --fpr 0.001 0.01 0.05 0.1 \
            --epochs 150")


rule doublets_call_af:
    input:
        filt_h5ad="piscem_spliceu/{run}/{run}.h5ad"
    output:
        scrublet_calls="scrublet/{run}/{run}_scrublet_calls_nc.tsv",
        h5ad="scrublet/{run}/{run}_initial_annotation_nc.h5ad"
    params:
        expected_dblt=0.1,
        sample_run_name="{run}",
        plots="output/figures/{run}_raw/doublets_call_nc"
    container:
        "docker://etretiakov/scrna-seq:jammy-2022.12.09-v0.0.1"
    threads: 8
    resources:
        mem_mb=20000
    script:
        "../code/scrublet_cb.py"


rule exploratory_data_analysis:
    input:
        rmd="analysis/01-eda-whole_dataset-nc.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.001.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-nc.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-nc.pdf",
        f"output/tables/01-eda-whole_dataset-nc/{prj}_all_mrk-MAST_sct-combined-whole_dataset-nc.csv",
        f"output/tables/01-eda-whole_dataset-nc/{prj}_all_mrk-logreg_sct-combined-whole_dataset-nc.csv",
        "output/tables/01-eda-whole_dataset-nc/parameters.json",
        f"data/{bprj}-whole_dataset-nc-clusters.h5Seurat"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 6
    resources:
        mem_mb=get_mem_mb
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")



rule convert_seurat_to_h5ad:
    input:
        h5ad="data/{bprj}-whole_dataset-fpr_{res}-clusters.h5Seurat"
    output:
        h5ad="data/{bprj}-whole_dataset-fpr_{res}-clusters.h5ad"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 2
    resources:
        mem_mb=30000
    script:
        "../code/convert_h5ad.R"

rule convert_nc_seurat_to_h5ad:
    input:
        h5ad=f"data/{bprj}-whole_dataset-nc-clusters.h5Seurat"
    output:
        h5ad=f"data/{bprj}-whole_dataset-nc-clusters.h5ad"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 2
    resources:
        mem_mb=30000
    script:
        "../code/convert_h5ad.R"


rule subset_astrocytes:
    input:
        "data/{bprj}-whole_dataset-fpr_{res}-clusters.h5ad"
    output:
        h5ad_annotations_all="data/class_cello/{bprj}-whole_dataset-{res}-cello_annotation.h5ad",
        tables_annotations_all="output/tables/class_cello/{bprj}-whole_dataset-{res}-CellO_output.tsv",
        h5ad_annotations_astrocytes="data/class_cello/{bprj}-astrocytes_dataset-{res}-initial_selection.h5ad",
        tables_annotations_astrocytes="output/tables/class_cello/{bprj}-astrocytes_dataset-{res}-initial_selection.tsv"
    params:
        prj=prj,
        bioprj=bprj,
        res="{res}"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 2
    resources:
        mem_mb=30000
    script:
        "../code/class_cello.py"


rule subset_nc_astrocytes:
    input:
        "data/{bprj}-whole_dataset-nc-clusters.h5ad"
    output:
        h5ad_annotations_all="data/class_cello/{bprj}-whole_dataset-nc-cello_annotation.h5ad",
        tables_annotations_all="output/tables/class_cello/{bprj}-whole_dataset-nc-CellO_output.tsv",
        h5ad_annotations_astrocytes="data/class_cello/{bprj}-astrocytes_dataset-nc-initial_selection.h5ad",
        tables_annotations_astrocytes="output/tables/class_cello/{bprj}-astrocytes_dataset-nc-initial_selection.tsv"
    params:
        prj=prj,
        bioprj=bprj,
        res="nc"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 2
    resources:
        mem_mb=30000
    script:
        "../code/class_cello.py"


rule correlations:
    input:
        h5ad="data/class_cello/{bprj}-astrocytes_dataset-0.001-initial_selection.h5ad"
    output:
        expr_mtx="data/{bprj}-astrocytes_dataset-0.001-expr-mtx.csv",
        pdf="output/figures/stat-corrmatrix-plt_{bprj}_.pdf"
    params:
        prj=prj,
        bioprj=bprj
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2023.04.30-custom-11.7"
    threads: 32
    resources:
        mem_mb=30000
    shell:
        ("quarto render analysis/03-correlations.qmd --to html")


rule grn:
    input:
        expr_mtx="data/{bprj}-astrocytes_dataset-0.001-expr-mtx.csv"
    output:
        modules="data/{bprj}-astrocytes_dataset-0.001.adjacencies.tsv"
    params:
        prj=prj,
        bioprj=bprj
    container:
        "docker://aertslab/pyscenic:0.12.1"
    threads: 32
    resources:
        mem_mb=120000
    shell:
        ("pyscenic grn {input.expr_mtx} /data/data/mm_tfs.csv -o {output.modules} --num_workers {threads}")


rule motifs:
    input:
        expr_mtx="data/{bprj}-astrocytes_dataset-0.001-expr-mtx.csv",
        modules="data/{bprj}-astrocytes_dataset-0.001.adjacencies.tsv"
    output:
        motifs="data/{bprj}-astrocytes_dataset-0.001.motifs.csv"
    params:
        prj=prj,
        bioprj=bprj
    container:
        "docker://aertslab/pyscenic:0.12.1"
    threads: 64
    resources:
        mem_mb=200000
    shell:
        ("pyscenic ctx {input.modules} /data/data/mm9-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather /data/data/mm9-tss-centered-5kb-10species.mc9nr.genes_vs_motifs.rankings.feather /data/data/mm9-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname /data/data/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname {input.expr_mtx} --output {output.motifs} --num_workers {threads}")


rule aucell:
    input:
        h5ad="data/class_cello/{bprj}-astrocytes_dataset-0.001-initial_selection.h5ad",
        motifs="data/{bprj}-astrocytes_dataset-0.001.motifs.csv"
    output:
        regulons="data/{bprj}-astrocytes_dataset-0.001.regulons.dat",
        auc_mtx="data/{bprj}-astrocytes_dataset-0.001-auc-mtx.csv",
        h5ad_scenic="data/{bprj}-astrocytes_dataset-0.001-scenic_plus.h5ad",
        h5ad_regulons="data/{bprj}-astrocytes_dataset-0.001-regulons.h5ad",
        gephi="output/tables/{bprj}-astrocytes_dataset-0.001-graph-regulons.graphml"
    params:
        prj=prj,
        bioprj=bprj
    container:
        "docker://aertslab/pyscenic_scanpy:0.12.1_1.9.1"
    threads: 32
    resources:
        mem_mb=200000
    script:
        "../code/aucell.py"
