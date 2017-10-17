import pandas as pd
"""script to generate summary annotation table from individual tables"""
def make_summary_table(out_table, fs):
    """read and combine tables"""
    with open(fs, 'r') as sum_tables:
        tables = []
        for f in sum_tables:
            t = pd.read_csv(f.strip(), index_col = 0)
            tables.append(t)
    summary_table = pd.concat(tables)
    summary_table.to_csv(out_table)
