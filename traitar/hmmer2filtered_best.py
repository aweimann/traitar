"""parse hmmer output file and generated filtered best file"""
import sys
import pandas as ps
from StringIO import StringIO
hmmer_colnames = ['target name','target accession','tlen','query name','accession','qlen','E-value','score_overall','bias_overall','#','of','c-Evalue','i-Evalue','score_domain','bias_domain','from_hmm','to_hmm','ali_from','ali_to','env_from','env_to','acc','description of target']

def filter_dbcan(m):
    """apply thresholds as suggested by the dbCAN authors i.e. a bit score threshold of 25. If the alignment is less than 80 bps long use an e-value threshold of 0.001 and or an evalue threshold of 0.00001 if the alignment is longer than 80 bps """
    return ((m.loc[:,"score_domain"] >= 25) & ((m.loc[:,"i-Evalue"] <= 0.001) & (m.loc[:,"ali_to"] - m.loc[:, "ali_from"] <= 80) | (m.loc[:,"i-Evalue"] <= 0.00001) & (m.loc[:,"ali_to"] - m.loc[:, "ali_from"] > 80))) 

def filter_any(m):
    """default bit score threshold of 25 and evalue threshold of 0.01"""
    #TODO check thresholds
    return (m.loc[:,"i-Evalue"] <= 0.01) & (m.loc[:,"score_domain"] >= 25)

def apply_thresholds(infile_f, hmm_name, out_filt_f, out_excl_f):
    """parse HMMER output file and apply thresholds for e-valu and bit score)"""
    #preparse lines by replacing white space delimitation by tabs
    #skip header
    infile_f.readline()
    infile_f.readline()
    infile_f.readline()
    #replace whitespace by tabs and skip lines which start with a # char
    cleaned = "".join(filter(lambda x: not x.startswith("#"), ["\t".join(i.split(None, 22)) for i in infile_f.readlines()]))
    #read  tab delimited hmmer output file with pandas via stringIO
    #account for cases where hmmer didn't return any hits
    try:
        m = ps.read_csv(StringIO(cleaned), sep = "\t",  header = None)
    except ValueError:
        m = ps.DataFrame(columns = hmmer_colnames)
    m.columns = hmmer_colnames
    if hmm_name == "dbcan":
        keep = filter_dbcan(m)
    else:
        keep = filter_any(m)
    #apply eval threshold
    m_eval = m.loc[keep, :]
    m_eval = m.loc[keep, :]
    if not out_filt_f is None: 
        m_eval.to_csv(out_filt_f, sep = "\t")
    m_eval_excl = m.loc[~keep, :]    
    if not out_excl_f is None:
        m_eval_excl.to_csv(out_excl_f, sep = "\t")
    return m_eval
    
    
def aggregate_domain_hits(filtered_df, out_f):
    """look for best scoring domain in each gene"""
    #sort by gene identifier and Pfam
    with open(out_f, 'w') as out_fo:
        ps.DataFrame(filtered_df.columns).T.to_csv(out_f, sep = "\t", index = False, header = False, mode = 'a')
        filtered_df.sort_index(by = ["target name", "query name"], inplace = True)
        if filtered_df.shape[0] > 0:
            current_max = filtered_df.iloc[0,] 
        else:
            #nothing todo
            return
        for i in range(1, filtered_df.shape[0]):
            if current_max.loc["query name"] != filtered_df.iloc[i,].loc["query name"] or current_max.loc["target name"] != filtered_df.iloc[i,].loc["target name"]:
                ps.DataFrame(current_max).T.to_csv(out_f, sep = "\t", index = False, header = False, mode = 'a')
                current_max = filtered_df.iloc[i,]
            else: 
                if current_max.iloc[13] < filtered_df.iloc[i,13]:
                    current_max = filtered_df.iloc[i,]
        #write last domain hit to disk
        ps.DataFrame(current_max).T.to_csv(out_f, sep = "\t", index = False, header = False, mode = 'a')
