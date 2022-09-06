import scvelo as scv
import tb_rna_velo.utils as utils


'''
ScVelo class
adata
gtf
'''
class ScVelo:
    def __init__(self, adata):
        self.adata = adata
        self.initial_filtering()
        self.show_proportions()

    def initial_filtering(self):
        scv.pp.filter_and_normalize(self.adata, min_shared_counts=10, n_top_genes=20000)
        scv.pp.moments(self.adata, n_pcs=30, n_neighbors=30, method="gauss")

    def show_proportions(self):
        self.adata.var_names_make_unique()
        scv.pl.proportions(self.adata)
    
    def get_adata(self):
        return self.adata
        
        
        