import numpy as np
import pandas as pd
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats


def test():
    print("this is a test function")
    

def remove_na(adata):
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
    
    return None


def find_ratios_sum(adata):
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
    
    return None


def find_ratios_std(adata):
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
    
    return None


def sns_scatter(
    x=None,
    y=None,
    data=None,
    is_log_scale=False,
    selected_genes=None,
    **kwargs
):
    sns.scatterplot(x=x, y=y, data=data, **kwargs)
    if selected_genes != None:
        df = pd.DataFrame(data, index=selected_genes, columns=data.columns)
        sns.scatterplot(x=x, y=y, data=df, **kwargs)
    if is_log_scale == True:
        plt.yscale('log')
        
    return None

def get_p_values(adata):    
    p_list = []
    x,y = adata.shape
    for n in range(0, y):
        a = adata.layers['spliced'][:,n]
        b = adata.layers['unspliced'][:,n]
        ks_test = stats.ks_2samp(a, b)
        p_list.append(ks_test[1])
    
    adata.var['p_value'] = p_list
    return None
    
    









    