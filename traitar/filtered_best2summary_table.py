#!/usr/bin/env python
"""script to create a summary matrix and a gene2hmm mapping from filtered and aggregated hmmer output files"""
from .PhenotypeCollection import PhenotypeCollection 
import pandas as ps
def gene2hmm(domtblout_fs, pt_models, gene2hmm_out = None, is_gene2hmm = False):
    """function to create a summary matrix and a gene2hmm mapping from filtered and aggregated hmmer output files"""
    #read accession file 
    accs = pt_models.get_pf2desc()
    f = open(domtblout_fs, 'r')
    domtblout_list = [i.strip() for i in f.readlines()] 
    f.close()
    sum_df = ps.DataFrame(ps.np.zeros((len(domtblout_list), accs.shape[0])))
    #set index to the actual files by getting rid of the preceding path
    sum_df.index = [i.split("/")[-1].replace("_filtered_best.dat", "") for i in domtblout_list]
    sum_df.columns = accs.index 
    #create a gene2hmm dictionary
    gene2hmm = {}  
    for f in domtblout_list:
        gene_df = ps.read_csv(f, sep = "\t")
        if not f in gene2hmm:
            f_mod = f.split("/")[-1].replace("_filtered_best.dat", "")
            gene2hmm[f_mod] = {}
        for i in gene_df.index:
            #append genes to gene2hmm dict
            gene = gene_df.loc[i, 'target name']
            if pt_models.get_hmm_name() == "pfam":
                query = gene_df.iloc[i, 4].split('.')[0]
            if pt_models.get_hmm_name() == "dbcan":
                query = gene_df.iloc[i, 3].split('.')[0]
            if is_gene2hmm:
                if gene not in gene2hmm[f_mod]:
                    gene2hmm[f_mod][gene] = [query] 
                else:
                    gene2hmm[f_mod][gene].append(query)
            #add annotations to summary file
            try:
                sum_df.loc[f_mod, query] += 1
            except KeyError:
                print query
    #write annotation summary to disk
    if not is_gene2hmm:
        sum_df.to_csv(gene2hmm_out, sep = "\t")
    #write gene2hmm to disk 
    #for i in gene2hmm:
    #    with open(gene2hmm_out, 'w') as out: 
    #        out.write("%s\t%s\n" (i, ",".join(gene2hmm[i])))
    return sum_df, gene2hmm
