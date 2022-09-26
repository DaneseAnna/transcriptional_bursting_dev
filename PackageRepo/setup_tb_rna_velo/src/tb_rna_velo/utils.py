import numpy as np
import pandas as pd
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import scvelo as scv


def test():
    print("this is a test function")
    

def remove_na(adata: ad.AnnData):
    """remove NA values from given anndata in slpiced and unspliced layers

    Args:
        adata (ad.AnnData): Anndata to be modified
    """
    allele_1 = pd.DataFrame(adata.layers['spliced'], columns=adata.var.index)
    allele_2 = pd.DataFrame(adata.layers['unspliced'], columns=adata.var.index)
    total = pd.DataFrame(adata.X, columns=adata.var.index, index=adata.obs.index)
    
    allele_1 = allele_1.dropna(axis=1)
    allele_2 = allele_2.dropna(axis=1)
    total = total.dropna(axis=1)
    
    common_genes1 = allele_1.columns.intersection(allele_2.columns)
    common_genes2 = allele_2.columns.intersection(total.columns)
    common_genes = common_genes1.intersection(common_genes2)
    
    adata._inplace_subset_var(common_genes)
    
    adata.X = total[common_genes]
    adata.layers['spliced'] = allele_1[common_genes]
    adata.layers['unspliced'] = allele_2[common_genes]
    
    return


def find_ratios_sum(adata: ad.AnnData):
    """find ratio sum of counts in spliced nd unspliced layers

    Args:
        adata (ad.AnnData): Anndata to be modified
    """
    allele_1 = pd.DataFrame(adata.layers['spliced'], columns=adata.var.index)
    allele_2 = pd.DataFrame(adata.layers['unspliced'], columns=adata.var.index)
    
    allele_1_sum = allele_1.sum(axis=0)
    allele_2_sum = allele_2.sum(axis=0)
    
    allele_1_sum_ratio = allele_1_sum/(allele_1_sum + allele_2_sum)
    allele_2_sum_ratio = allele_2_sum/(allele_1_sum + allele_2_sum)
    
    allele_1_sum_ratio = allele_1_sum_ratio.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    allele_2_sum_ratio = allele_2_sum_ratio.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    
    common_genes = allele_1_sum_ratio.index.intersection(allele_2_sum_ratio.index)
    
    adata._inplace_subset_var(common_genes)
    
    adata.var["sum_allele_1"] = allele_1_sum[common_genes]
    adata.var["sum_allele_2"] = allele_2_sum[common_genes]
    adata.var["ratio_allele_1"] = allele_1_sum_ratio[common_genes]
    adata.var["ratio_allele_2"] = allele_2_sum_ratio[common_genes]
    
    return


def find_ratios_std(adata: ad.AnnData):
    """find standard deviationratios of counts in spliced and unslpiced layers

    Args:
        adata (ad.AnnData): Anndata to be modified
    """
    allele_1 = pd.DataFrame(adata.layers['spliced'], columns=adata.var.index)
    allele_2 = pd.DataFrame(adata.layers['unspliced'], columns=adata.var.index)
    
    allele_1_ratio = allele_1 / (allele_1 + allele_2)
    allele_2_ratio = allele_2 / (allele_1 + allele_2)

    
    allele_1_ratio_std = allele_1_ratio.std(axis=0)
    allele_2_ratio_std = allele_2_ratio.std(axis=0)
    
    allele_1_ratio_std = allele_1_ratio_std.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    allele_2_ratio_std = allele_2_ratio_std.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    
    common_genes = allele_1_ratio_std.index.intersection(allele_2_ratio_std.index)

    
    adata._inplace_subset_var(common_genes)
    
    #adata.layers["ratio_allele_1"] = allele_1_ratio.loc[:, [common_genes]]
    #adata.layers["ratio_allele_2"] = allele_2_ratio.loc[:, [common_genes]]
    
    adata.layers["ratio_allele_1"] = allele_1_ratio[common_genes]
    adata.layers["ratio_allele_2"] = allele_2_ratio[common_genes]
    
    adata.var["ratio_sum_allele_1"] = allele_1_ratio.sum(axis=0)[common_genes]
    adata.var["ratio_sum_allele_2"] = allele_2_ratio.sum(axis=0)[common_genes]
    adata.var["ratio_mean_allele_1"] = allele_1_ratio.mean(axis=0)[common_genes]
    adata.var["ratio_mean_allele_2"] = allele_2_ratio.mean(axis=0)[common_genes]
    adata.var["ratio_std_allele_1"] = allele_1_ratio_std[common_genes]
    adata.var["ratio_std_allele_2"] = allele_2_ratio_std[common_genes]
    
    return


def find_count_per_gene(adata: ad.AnnData):
    adata.obs['n_counts'] = adata.X.sum(1)
    adata.var['UMI_counts_per_gene'] = adata.X.sum(0)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.var['UMI_log_counts_per_gene'] = np.log(adata.var['UMI_counts_per_gene'])


def find_max_counts_per_gene_per_allele(adata: ad.AnnData):
    allele_1 = adata.to_df(layer="spliced").transpose()
    allele_2 = adata.to_df(layer="unspliced").transpose()
    max_count_allele_1 = []
    max_count_allele_2 = []
    for x in adata.var.index:
        max_count_allele_1.append(np.max(allele_1.loc[x]))
        max_count_allele_2.append(np.max(allele_2.loc[x]))

    adata.var['max_count_allele_1'] = max_count_allele_1
    adata.var['max_count_allele_2'] = max_count_allele_2
    return


def sns_scatter(
    x:str=None,
    y:str=None,
    data=None,
    is_log_scale:bool=False,
    selected_genes=None,
    **kwargs
):
    """creates sns scatter plot, modified to possibly highlight given set of genes and display in log scale

    Args:
        x (_type_, optional): _description_. Defaults to None.
        y (str, optional): _description_. Defaults to None.
        data (_type_, optional): _description_. Defaults to None.
        is_log_scale (bool, optional): _description_. Defaults to False.
        selected_genes (_type_, optional): _description_. Defaults to None.
    """
    sns.scatterplot(x=x, y=y, data=data, **kwargs)
    if selected_genes != None:
        df = pd.DataFrame(data, index=selected_genes, columns=data.columns)
        sns.scatterplot(x=x, y=y, data=df, **kwargs)
    if is_log_scale == True:
        plt.yscale('log')
        
    return

def get_p_values(adata: ad.AnnData):    
    """generates P-values for ks statistics using spliced and unspliced layers

    Args:
        adata (ad.AnnData): _description_
    """
    p_list = []
    x,y = adata.shape
    for n in range(0, y):
        a = adata.layers['spliced'][:,n]
        b = adata.layers['unspliced'][:,n]
        ks_test = stats.ks_2samp(a, b)
        p_list.append(ks_test[1])
    
    adata.var['p_value'] = p_list
    return


def scatter_plot(adata=None, x=None, xlabel=None, ylabel=None, gene_name=None, UMI_counts=None, title=None, color=None):
    if adata is None:
        print("Requires data field")
        return

    if x is None:
        print("Requires x field")
        return

    if gene_name is None:
        gene_name = adata.var.loc[x]['gene_name']

    if UMI_counts is None:
        UMI_counts = adata.var.loc[x]['sum_allele_1'] + adata.var.loc[x]['sum_allele_2']

    if title is None:
        title = f"{gene_name} ({x}): UMI counts sum = {UMI_counts}"

    # add some boolean for equal scale
    max_count = np.max([adata.var.loc[x]['max_count_allele_1'], adata.var.loc[x]['max_count_allele_2']])
    pad = int(0.2 * max_count)
    x_lim = (-pad, max_count + pad)
    y_lim = (-pad, max_count + pad)

    scv.pl.scatter(adata, x, xlabel=xlabel, ylabel=ylabel, title=title, xlim=x_lim, ylim=y_lim, color = color)

    
    









    