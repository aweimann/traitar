import pandas as ps
def flatten_df(df1, df2, name1, name2, out):
    with open(out, 'w') as f:
        f.write("sample\tphenotype\tphenotype model\n")
        for i in set(df1.index).union(set(df2.index)):
            for j in set(df1.columns).union(set(df2.columns)):
                if i in df1.index and j in df1.columns and  not df1.loc[i,j] == 0:
                    f.write("\t".join([str(i), str(j), str(df1.loc[i,j]), name1]) + "\n")
                if i in df2.index and j in df2.columns and not df2.loc[i,j] == 0:
                    f.write("\t".join([str(i), str(j), str(df2.loc[i,j]), name2]) + "\n")

def comb_preds(phypat_dir, phypat_PGL_dir, primary_name, secondary_name, out_dir, k):
    #phypat only predictions
    m1_maj = ps.read_csv("%s/predictions_majority-vote.txt"%phypat_dir, index_col = 0, sep = "\t")
    m1 = ps.read_csv("%s/predictions_single-votes.txt"%phypat_dir, index_col = 0, sep = "\t")
    m2_maj = ps.read_csv("%s/predictions_majority-vote.txt"%phypat_PGL_dir, index_col = 0, sep = "\t")
    m2 = ps.read_csv("%s/predictions_single-votes.txt"%phypat_PGL_dir, index_col = 0, sep = "\t")
    #write to disk a single vote version of the predictions
    #phypat plus ml predictions
    #set NAs to zero
    m1[ps.isnull(m1)] = 0
    m2[ps.isnull(m2)] = 0
    #collect phenotypes that are shared / different in the two data sets
    #combine predictions of phypat and phypat+PGL into one matrix 
    m_columns = set(m1.columns).union(set(m2.columns))
    m = ps.DataFrame(ps.np.zeros((m1.shape[0], len(m_columns))))
    m_maj = ps.DataFrame(ps.np.zeros((m1.shape[0], len(m_columns))))
    m_maj.columns = m.columns = m_columns
    m_maj.index = m1.index
    m.index = m1.index
    #combine single voters
    m.loc[:, m1.columns] = m1
    m.loc[:, m2.columns] = m.loc[:, m2.columns] + m2
    #combine majority votes
    m_maj.loc[:, set(m1.columns).difference(set(m2.columns))] = m1_maj.loc[:, set(m1.columns).difference(set(m2.columns))] * 2
    m_maj.loc[:, set(m2.columns).difference(set(m1.columns))] = m2_maj.loc[:, set(m2.columns).difference(set(m1.columns))] 
    m_temp = m_maj.loc[:, set(m2.columns).intersection(set(m1.columns))]
    #sys.stderr.write(str((m1_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1) & (m2_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1)))
    if not set(m2.columns).intersection(set(m1.columns)) == set():
        m_temp[(m1_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1) & (m2_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1)] = 3
        m_temp[(m1_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1) & (m2_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 0)] = 1
        m_temp[(m1_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 0) & (m2_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1)] = 2
        m_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] = m_temp
    m.to_csv("%s/predictions_single-votes_combined.txt"%out_dir, sep = "\t")
    m_maj.to_csv("%s/predictions_majority-vote_combined.txt"%out_dir, index_label = None, sep = "\t")
    flatten_df(m1_maj, m2_maj, primary_name, secondary_name, "%s/predictions_flat_majority-votes_combined.txt"%out_dir)
    flatten_df(m1, m2, primary_name, secondary_name, "%s/predictions_flat_single-votes_combined.txt"%out_dir)
