from setuptools import setup, find_packages

setup(
    name='tb_rna_velo',
    version='0.8',
    license='TEST',
    author="RijeshShrestha",
    author_email='craaabby@gmail.com.com',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    keywords='example project',
    install_requires=[
        'setuptools==59.8.0',
        'tensorflow',
        'anndata',
        'matplotlib==3.5.1',
        'pandas',
        'scikit_learn',
        'scipy',
        'scvelo',
        'seaborn'        
    ]
)