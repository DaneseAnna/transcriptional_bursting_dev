
- might need to create new conda env with tensor-flow dependency (needed for Mac with silicon chips: https://towardsdatascience.com/installing-tensorflow-on-the-m1-mac-410bb36b776)
>conda env create --file=Documents/environment.yml --name scvelo


- install the package from PyPi
> pip install tb-rna-velo==<version>


- might need to install Leiden algorithm dependency using the following command or pip alternative:
> conda install -c conda-forge leidenalg