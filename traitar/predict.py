import os.path
import pandas as pd
import sys
from .PhenotypeCollection import PhenotypeCollection

"""predict new samples"""

def filter_pred(scores, is_majority, k):
    """either do majority vote aggregation or conservative all or nothing vote"""
    if is_majority:
        if (scores > 0).sum() >= (k/2 + 1):
            return scores[scores >= 0].mean()
        else:
            return pd.np.NaN
    else:
        if (scores > 0).all():
            return scores.mean()
        else:
            return pd.np.NaN

def aggregate(pred_df, k, out_dir, pt2acc):
    """employ different prediction strategies"""
    out = ["majority-vote", "conservative-vote", "single-votes"]
    maj_pred_dfs = pd.DataFrame(pd.np.zeros(shape = (pred_df.shape[0], pred_df.shape[1] / k)), columns = ["%s" % i for i in range(pred_df.shape[1]/k)])
    maj_pred_dfs.index = pred_df.index
    maj_pred_dfs_columns = maj_pred_dfs.columns.tolist()
    for i in range(pred_df.shape[1] / k):
        maj_pred_dfs_columns[i] = pred_df.columns.values[i * k].split("_")[0] 
        maj_pred_dfs.iloc[:, i] = pred_df.iloc[:, (i * k) : (i * k + k)].apply(lambda x: (x > 0).sum(), axis = 1).astype('int')
    maj_pred_dfs.columns = pt2acc.loc[maj_pred_dfs_columns, :].iloc[:, 0]
    #majority vote
    maj_pred_dfs.to_csv("%s/predictions_%s.txt"%(out_dir, out[2]), sep = "\t", float_format='%.0f', encoding = "utf-8")
    (maj_pred_dfs >= k/2 + 1).astype('int').to_csv("%s/predictions_%s.txt"%(out_dir, out[0]), sep = "\t", float_format='%.0f', encoding = "utf-8")
    (maj_pred_dfs == k).astype('int').to_csv("%s/predictions_%s.txt"%(out_dir, out[1]), sep = "\t", float_format='%.0f', encoding = "utf-8")   #conservative vote
    return maj_pred_dfs
    

    
def majority_predict(pt, pt_models, test_data, k, bias_weight = 1):
    """predict the class label based on a committee vote of the models in models""" 
    #TODO if the classification model is trained on non binary data we would require the same normalization applied to the training data 
    #binarize
    if not pt_models.is_standardized:
        test_data_n = (test_data > 0).astype('int')
    #check if classification model exists for the requested phenotype
    try: 
        bias = pt_models.get_bias(pt)
        predictors = pt_models.get_predictors(pt)
    except KeyError:
        return pd.DataFrame()
    #build prediction matrix with the k best models
    preds = pd.np.zeros((test_data.shape[0], k))
    for i in range(k):
        if pt_models.is_standardized:
            scale, mean = pt_models.get_scale_and_mean(pt)
            from sklearn import preprocessing
            scaler = preprocessing.StandardScaler()
            scaler.mean_ = mean.T
            scaler.scale_ = scale.T
            test_data_n = pd.DataFrame(scaler.transform(test_data))
            test_data_n.index = test_data.index
            test_data_n.columns = test_data.columns
        preds[:, i] = bias.iloc[i, 0] * bias_weight +  predictors.iloc[:, i].dot(test_data_n.loc[:, predictors.iloc[:, i].index].T)
        pred_df = pd.DataFrame(preds, index = test_data.index)
    #set column names
    pred_df.columns = ["%s_%s" %(pt, predictors.columns[i].split("_")[0]) for i in range(k)]
    return pred_df

def annotate_and_predict(pt_models, summary_f, out_dir, k):
    """load annotation previously generated with HMMER and predict phenotypes with phypat and phypat+GGL models"""
    pred_df = pd.DataFrame()
    #read pfam to description file
    pfam_mapping = pt_models.get_pf2desc()
    pt_mapping = pt_models.get_pt2acc()
    #read annotation file
    m = pd.read_csv(summary_f, sep="\t", index_col = 0)
    #restrict to those pfams contained in the model
    #m_red = m.loc[:, pfam_mapping.index ].astype('int')
    #m.columns = pfam_mapping.index
    for pt in pt_mapping.index:
        #predict an append to prediction df
        preds = majority_predict(pt, pt_models, m, k)
        pred_df = pd.concat([pred_df, preds], axis = 1)
    pred_df.index = m.index
    pred_df.to_csv("%s/predictions_raw.txt"%out_dir, sep = "\t", float_format='%.3f')
    #aggregate predictions
    aggr_dfs = aggregate(pred_df, k, out_dir, pt_mapping)
    return pred_df
