import numpy as np
import pandas as pd
import anndata as ad
import scvelo as scv
import scanpy as sc
from scipy import stats


def remove_na(adata: ad.AnnData, layers: list=None):
    """remove NA values from given anndata in allele specific data

    Args:
        adata (ad.AnnData): Anndata to be modified
        layers (list): List of layers containing allele specific data
    """
    if layers is None:
        allele_1_layer = 'allele_1'
        allele_2_layer = 'allele_2'
    else:
        allele_1_layer = layers[0]
        allele_2_layer = layers[1]

    
    
    allele_1 = pd.DataFrame(adata.layers[allele_1_layer], columns=adata.var.index)
    allele_2 = pd.DataFrame(adata.layers[allele_2_layer], columns=adata.var.index)

    allele_1 = allele_1.dropna(axis=1)
    allele_2 = allele_2.dropna(axis=1)

    common_genes = allele_1.columns.intersection(allele_2.columns)

    adata._inplace_subset_var(common_genes)
    
    adata.layers[allele_1_layer] = allele_1[common_genes]
    adata.layers[allele_2_layer] = allele_2[common_genes]

    return


def find_ratios_sum(adata: ad.AnnData, layers: list=None):
    """find ratio sum of counts in allele specific data

    Args:
        adata (ad.AnnData): Anndata to be modified
        layers (list): List of layers containing allele specific data
    """
    if layers is None:
        allele_1_layer = 'allele_1'
        allele_2_layer = 'allele_2'
    else:
        allele_1_layer = layers[0]
        allele_2_layer = layers[1]

    allele_1 = pd.DataFrame(adata.layers[allele_1_layer], columns=adata.var.index)
    allele_2 = pd.DataFrame(adata.layers[allele_2_layer], columns=adata.var.index)

    allele_1_sum = allele_1.sum(axis=0)
    allele_2_sum = allele_2.sum(axis=0)

    allele_1_sum_ratio = allele_1_sum / (allele_1_sum + allele_2_sum + np.finfo(float).eps)
    allele_2_sum_ratio = allele_2_sum / (allele_1_sum + allele_2_sum + np.finfo(float).eps)

    allele_1_sum_ratio = allele_1_sum_ratio.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    allele_2_sum_ratio = allele_2_sum_ratio.replace([np.inf, -np.inf], np.nan).dropna(axis=0)

    common_genes = allele_1_sum_ratio.index.intersection(allele_2_sum_ratio.index)

    adata._inplace_subset_var(common_genes)

    adata.var["sum_" + allele_1_layer] = allele_1_sum[common_genes]
    adata.var["sum_" + allele_2_layer] = allele_2_sum[common_genes]
    adata.var["sum_ratio_" + allele_1_layer] = allele_1_sum_ratio[common_genes]
    adata.var["sum_ratio_" + allele_2_layer] = allele_2_sum_ratio[common_genes]

    return


def find_ratios_std(adata: ad.AnnData, layers: list=None):
    """find standard deviationratios of count ratios in allele specific data

    Args:
        adata (ad.AnnData): Anndata to be modified
        layers (list): List of layers containing allele specific data
    """
    if layers is None:
        allele_1_layer = 'allele_1'
        allele_2_layer = 'allele_2'
    else:
        allele_1_layer = layers[0]
        allele_2_layer = layers[1]

    allele_1 = pd.DataFrame(adata.layers[allele_1_layer], columns=adata.var.index)
    allele_2 = pd.DataFrame(adata.layers[allele_2_layer], columns=adata.var.index)

    allele_1_ratio = allele_1 / (allele_1 + allele_2 + np.finfo(float).eps)
    allele_2_ratio = allele_2 / (allele_1 + allele_2 + np.finfo(float).eps)

    allele_1_ratio_std = allele_1_ratio.std(axis=0)
    allele_2_ratio_std = allele_2_ratio.std(axis=0)

    allele_1_ratio_std = allele_1_ratio_std.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    allele_2_ratio_std = allele_2_ratio_std.replace([np.inf, -np.inf], np.nan).dropna(axis=0)

    common_genes = allele_1_ratio_std.index.intersection(allele_2_ratio_std.index)
    adata._inplace_subset_var(common_genes)


    adata.layers["ratio_" + allele_1_layer] = allele_1_ratio[common_genes]
    adata.layers["ratio_" + allele_2_layer] = allele_2_ratio[common_genes]

    adata.var["ratio_sum_" + allele_1_layer] = allele_1_ratio.sum(axis=0)[common_genes]
    adata.var["ratio_sum_" + allele_2_layer] = allele_2_ratio.sum(axis=0)[common_genes]

    adata.var["ratio_mean_" + allele_1_layer] = allele_1_ratio.mean(axis=0)[common_genes]
    adata.var["ratio_mean_" + allele_2_layer] = allele_2_ratio.mean(axis=0)[common_genes]

    adata.var["ratio_std_" + allele_1_layer] = allele_1_ratio_std[common_genes]
    adata.var["ratio_std_" + allele_2_layer] = allele_2_ratio_std[common_genes]

    return


def get_p_values(adata: ad.AnnData, layers: list=None):
    """generates P-values for ks statistics using allele specific data

    Args:
        adata (ad.AnnData): Anndata to be modified
        layers (list): List of layers containing allele specific data
    """
    if layers is None:
        allele_1_layer = 'allele_1'
        allele_2_layer = 'allele_2'
    else:
        allele_1_layer = layers[0]
        allele_2_layer = layers[1]

    p_list = []
    x, y = adata.shape
    for n in range(0, y):
        a = adata.layers[allele_1_layer][:, n]
        b = adata.layers[allele_2_layer][:, n]
        ks_test = stats.ks_2samp(a, b)
        p_list.append(ks_test[1])

    adata.var['allele_p_value'] = p_list
    return

def get_gtf_data(adata: ad.AnnData, file_path:str=""):
    """get information regarding gene structure using gtf file

    Args:
        adata (ad.AnnData): Anndata to be modified
        file_path: location of gtf file
    """

    if file_path == "":
        print("file path for the GTF file is required")
        return

    # get GTF data
    gtf = pd.read_csv(file_path, sep='\t', skiprows=5, header=None)
    gtf.columns = ['chr', 'database', 'type', 'start', 'end', '.', 'strand', '.', 'other']
    gtf = gtf[gtf['type'] == "gene"]
    gtf = gtf.reset_index()

    dic_gene_names = {}
    index = 0
    for line in gtf['other'].tolist():
        line = line.split(';')
        dic_gene_names[line[0][9:-1]] = [gtf['chr'][index], line[2][12:-1]]
        index += 1

    adata.var['gene_name'] = [dic_gene_names[x][1] if x in dic_gene_names.keys() else 'NA' for x in adata.var.Accession]
    adata.var['chromosome'] = [dic_gene_names[x][0] if x in dic_gene_names.keys() else 'NA' for x in adata.var.Accession]

    adata.var['chromosome'] = adata.var['chromosome'].astype('str') # need to specify as string to save data

    autosome_label = []
    nuclear_label = []
    ribosomal_protein_label=[]
    for chrom in adata.var['chromosome']:
        if chrom == "X":
            autosome_label.append('X')
        elif chrom == "Y":
            autosome_label.append('Y')
        else:
            autosome_label.append('autosome')

        if chrom == "MT":
            nuclear_label.append('MT')
        else:
            nuclear_label.append('nuclear')

    for genes in adata.var['gene_name']:
        if genes[:2] == 'Rp':
            ribosomal_protein_label.append('Ribosomal protein')
        else:
            ribosomal_protein_label.append('other protein')

    adata.var['autosomes'] = autosome_label
    adata.var['nuclear'] = nuclear_label
    adata.var['Ribosomal_prot'] = ribosomal_protein_label

    return

def get_pca_and_umap_data(adata: ad.AnnData):
    """run and calculate the pca, umap and get cluster using leiden algorithm

    Args:
        adata (ad.AnnData): Anndata to be modified
    """

    sc.pp.neighbors(adata, n_neighbors=8, n_pcs=31, use_rep='X', method='gauss')
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    adata.obs['clusters'] = adata.obs['leiden'].copy()
