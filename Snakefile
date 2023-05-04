from os.path import join, isfile, dirname, basename
from glob import glob
from collections import defaultdict

import pandas as pd

conf_file = "config.yaml"

if isfile(conf_file):

    configfile: conf_file


metadata = pd.read_csv("data/metadata.csv", sep="\t", index_col="Date Processed at CG")

celltype_annotations = {
    "celltype_cluster": {
        "CD4+ regulatory T cell",
	"Follicular helper CD4 T cell",
	"Naive CD4+ T cells",
	"CD4+ T-Helper 17 Cell",
	"CD4+ transitional memory T cell",
	"CD4+ central memory T cell",
	"CD4+ T-Helper 2 Cell",
	"CD4+ T-Helper 1 Cell",
	"Terminal effector CD4+ T cell",
	"Gamma delta T cell",
	"CD4+ effector memory T cell",
	"Mucosal associated invariant T cell",
	"CD8+ Terminal effector T cell",
	"CD8+ stem cell memory T cell",
	"CD8+ effector memory T cell",
	"Naive CD8+ T cell",
	"CD8+ central memory T cell",
	"Proliferating CD8+ T cell",
	"IgM Memory B cell",
	"Double negative 1 B cell",
	"Double negative 2 B cell",
	"Classical Memory B cell",
	"Naive B cell",
	"B cell",
	"Transitional B cell",
	"Double negative 3 B cell",
	"Conventional Dendritic cell Type 2",
	"CD16+ Monocyte",
	"Megacaryocytes",
	"Proliferating Monocytes",
	"CD14+ Monocyte",
	"Plasmacytoid Dendritic cell",
	"Conventional Dendritic cell Type 1",
	"NKT-like cell",
	"Proliferating CD4+ T cell",
	"Plasmablasts",
	"CD56- Natural Killer cell",
	"CD56+ Natural Killer cell"
    },
    "celltype": {
        "CD4+ T-Helper Cell",
	"Naive CD4+ T cell",
	"CD4+ Memory T cell",
	"Gamma delta T cell",
	"Mucosal associated invariant T cell",
	"CD8+ Memory T cell",
	"Naive CD8+ T cell",
	"Proliferating CD8+ T cell",
	"Memory B cell",
	"Double negative B cell",
	"Naive B cell",
	"B cell",
	"Transitional B cell",
	"Conventional Dendritic cell",
	"CD16+ Monocyte",
	"Megacaryocytes",
	"Proliferating Monocyte",
	"CD14+ Monocyte",
	"Plasmacytoid Dendritic cell",
	"NKT-like cell",
	"Proliferating CD4+ T cell",
	"Plasmablasts",
	"CD56-Dim NK cell",
	"CD56-Bright NK cell"
    }
}

celltype_annotations["celltype_w_batch_and_sex"] = celltype_annotations["celltype"]

CATEGORIES = config["enrichment"]["categories"]

COMPARISONS = ["NDvsHC", "CIvsHC", "MCIvsHC", "ADvsHC", "PDvsHC", "PDMCIvsHC",
               "ADvsPD", "AD_PDvsHC", "MCI_PDMCIvsHC", "ADvsMCI", "PDvsMCI",
               "AD_PDvsMCI", "ADvsMCI_PDMCI", "PDvsMCI_PDMCI", "AD_PDvsMCI_PDMCI"
]

SHAP_COMPARISONS = ["MCIvsHC", "ADvsHC", "PDvsHC", "PDMCIvsHC",
               "AD_PDvsHC", "MCI_PDMCIvsHC", "ADvsMCI", "PDvsMCI",
               "AD_PDvsMCI", "ADvsMCI_PDMCI", "PDvsMCI_PDMCI", "AD_PDvsMCI_PDMCI"
]

SETTING_CELL_POP_UP_ANNOT = defaultdict(list)
SETTING_CELL_POP_DOWN_ANNOT = defaultdict(list)
SETTING_CELL_POP_ALL_ANNOT = defaultdict(list)
SETTING_CELL_POP_GSEA_ANNOT = defaultdict(list)
SETTING_CELL_POP_GSEA_SHAP_ANNOT = defaultdict(list)
for c in COMPARISONS:
    for cannot, celltypes in celltype_annotations.items():
        for ctype in celltypes:
            SETTING_CELL_POP_UP_ANNOT[cannot].append("{}|{}|{}_up".format(cannot, c, ctype))
            SETTING_CELL_POP_DOWN_ANNOT[cannot].append("{}|{}|{}_down".format(cannot, c, ctype))
            SETTING_CELL_POP_ALL_ANNOT[cannot].append("{}|{}|{}_all".format(cannot, c, ctype))
            SETTING_CELL_POP_GSEA_ANNOT[cannot].append("{}|{}|{}".format(cannot, c, ctype))
            if c in SHAP_COMPARISONS:
                SETTING_CELL_POP_GSEA_SHAP_ANNOT[cannot].append("{}|{}|{}".format(cannot, c, ctype))

SETTING_CELL_POP = ["{}|{}|{}_{}".format(cannot, c, ctype, regdir) for c in COMPARISONS for cannot, celltypes in celltype_annotations.items() for ctype in celltypes for regdir in ["up", "down", "all"]]

TARGET_DBS = ["GO_-_Biological_Process", "GO_-_Cellular_Component", "GO_-_Molecular_Function", "KEGG_-_Pathways", "Reactome_-_Pathways", "Pfam_-_Protein_families", "MSigDB_v7_C6_Oncogenic_Signatures", "MSigDB_v7_C7_Immunologic_Signatures", "MSigDB_v7_H_Hallmark"]

TOOLS = ("limma_voom")

### No pre-Processing ! ####
rule all:
    input:
        expand(join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_all_comparisons_{celltype}.all.pdf"), celltype=("celltype_cluster"), tool=TOOLS), #, "celltype_w_batch_and_sex", "celltype_cluster", "celltype_cluster_w_batch_and_sex"



include: "rules/diff_exp.smk"


def get_raw_matrix_folder(wildcards):
    sample = wildcards.sample.split("_")[0]
    folder = glob(f"data/unfiltered_with_new/{sample}")
    #    if len(folder) == 0:
    #        folder = glob("data/unfiltered_with_new/{}".format(basename(dirname(dirname(metadata.loc[wildcards.sample, "matrix"])))))
    assert len(folder) == 1
    return folder[0]



rule clean_soup:
    input:
        raw=get_raw_matrix_folder,
        filtered=lambda w: metadata.loc[w.sample, "matrix"],
    output:
        counts=directory(join(config["results_dir"], "soupx/{sample}")),
        stats=join(config["results_dir"], "soupx/{sample}.stats.json"),
    conda:
        "envs/rstudio_env.yml"
    params:
        cont_fraction=lambda w: soup_cont.get(w.sample, None),
        min_genes_per_cell=300,
    threads: 1
    script:
        "scripts/create_soupX_cleaned_matrix.R"


rule filter_sample:
    input:
        filtered=join(config["results_dir"], "soupx/{sample}"),
    output:
        filtered=join(config["results_dir"], "filtered", "{sample}.rds"),
        stats=join(config["results_dir"], "filtered", "{sample}.stats.csv"),
    conda:
        "envs/rstudio_env.yml"
    params:
        min_cells_per_sample=300,
        max_mito_perc=7,
        min_genes_per_cell=300,
        max_genes_per_cell=4000,
    threads: 1
    script:
        "scripts/filter_sample_and_remove_doublets.R"


rule add_info_2_metadata:
    input:
        metadata="data/metadata.csv",
        filtered=expand(
            join(config["results_dir"], "filtered", "{sample}.rds"),
            sample=metadata.index,
        ),
        soupx_stats=expand(
            join(config["results_dir"], "soupx/{sample}.stats.json"),
            sample=metadata.index,
        ),
    output:
        metadata=join(config["results_dir"], "metadata.csv"),
    params:
        min_cells_per_sample=300,
        max_soup=0.4,
    threads: 16
    conda:
        "envs/rstudio_env.yml"
    script:
        "scripts/add_cells_and_soup_2_metadata.R"

rule create_per_sample_umap_and_ref_mapped_annot:
    input:
        rds=join(config["results_dir"], "filtered", "{sample}.rds"),
        ref="data/seurat/pbmc_multimodal.rds"
    output:
        rds=join(config["results_dir"], "azimuth_annot", "{sample}.rds"),
        png=join(config["results_dir"], "azimuth_annot", "{sample}.png")
    threads: 8
    script:
        "scripts/create_umap_w_azimuth_celltypes.R"

rule merge_samples:
    input:
        metadata="results/metadata.csv",
        geneid2symbol="data/annotations/geneid2symbol.gencode_v32.csv",
        filtered=expand(
            join(config["results_dir"], "filtered", "{sample}.rds"),
            sample=metadata.index,
        ),
    output:
        rds=join(config["results_dir"], "processed", "merged.rds"),
        h5ad=join(config["results_dir"], "processed", "merged.h5ad"),
        rds_unfiltered=join(config["results_dir"], "processed", "merged.unfiltered_genes.rds"),
        h5ad_unfiltered=join(config["results_dir"], "processed", "merged.unfiltered_genes.h5ad"),
    conda:
        "envs/rstudio_env.yml"
    threads: 16
    script:
        "scripts/merge_samples.R"


