import pandas as pd
import os
import sys
from traitar.PhenotypeCollection import PhenotypeCollection

class evaluate:
    
    @staticmethod
    def evaluate(out, gold_standard_f, traitar_pred_f, min_samples, are_pt_ids = True, phenotype_archive = None):
        """compare traitar predictions with a given gold standard"""
        #read in gold standard
        gs = pd.read_csv(gold_standard_f, index_col = 0, sep = "\t", na_values = "?", encoding = "utf-8" )
        gs.replace(["-", "+"], [0, 1], inplace = True)
        gs.replace(["ND"], [pd.np.nan], inplace = True)
        #check if gold_standard uses phenotype ids and replace with accessions in that case
        if are_pt_ids: 
            #read in phenotype mapping
            pc = PhenotypeCollection.PhenotypeCollection(phenotype_archive)
            pt_id2acc = pc.get_pt2acc()
            gs.columns = pt_id2acc.loc[gs.columns, :].iloc[:, 0]
        #read in traitar preds
        tp = pd.read_csv(traitar_pred_f, index_col = 0, sep = "\t", encoding = "utf-8")
        #get pts that are in gold standard and in the pt models 
        pts = list(set(gs.columns.tolist()).intersection(set(tp.columns.tolist())))
        #get list of samples that are in gold standard and in the pt models
        samples = list(set(gs.index.tolist()).intersection(set(tp.index.tolist())))
        gs = gs.loc[samples, :]
        tp = tp.loc[samples, :]
        #combine positive predictions
        #tp[(tp == 2.0) | (tp == 1.0) | (tp == 3.0)] = 1
        #conservative classification
        #tp[(tp == 2.0) | (tp == 1.0) ] = 0
        #either both or primary 
        #tp[(tp == 3.0)] = 1
        #tp[(tp == 2.0) ] = 0
        #either both or secondary 
        #tp[(tp == 1.0) ] = 0
        #tp[(tp == 3.0) ] = 1
        #tp[(tp == 2.0) ] = 1

        pt_gold_too_few_samples = gs.apply(lambda x: pd.Series(((x[~pd.isnull(x) & (x == 0)].sum() >= min_samples) & (x[~pd.isnull(x) & (x == 1)].sum() >= min_samples))))
        #print pd.concat([(gs > 0).sum(), (gs == 0).sum()], axis = 1)
        #print pt_gold_too_few_samples
        if len(pts) == 0:
            sys.exit("No phenotypes shared between traitar predictions and gold standard")
        #confusion matrix per phenotype
        conf_per_pt = pd.DataFrame(pd.np.zeros(shape = (len(pts), 4)))
        conf_per_pt.index = pts
        #performance measures per phenotype 
        perf_per_pt = pd.DataFrame(pd.np.zeros(shape = (len(pts), 4)))
        perf_per_pt.index = pts
        perf_per_pt.columns = ["recall pos", "recall neg", "macro accuracy", "precision"]
        for pt in pts:
            not_null = gs.loc[~pd.isnull(gs.loc[:, pt]),].index
            conf_per_pt.loc[pt, ] = evaluate.confusion_m(gs.loc[not_null, pt], tp.loc[not_null, pt])
            perf_per_pt.loc[pt, ] = evaluate.get_performance(conf_per_pt.loc[pt, ]) 
            miscl = evaluate.get_miscl(gs.loc[:, pt], tp.loc[:, pt])
            if not len(miscl) == 0:
                miscl_m = pd.concat([gs.loc[miscl, pt], tp.loc[miscl, pt]], axis = 1)
                miscl_m.columns = ["gold standard", "traitar predictions"]
                miscl_m.to_csv('%s.txt' % os.path.join(out, pt), sep = "\t")
            
        #sum up confusion
        overall_conf = conf_per_pt.sum(axis = 0)
        #micro averaged performance measures
        micro_perf = evaluate.get_performance(overall_conf)
        #macro averaged performance
        macro_perf = perf_per_pt.mean(axis = 0)
        #write to disk perf2pt and overall accuracy measures
        freq_per_pt = pd.concat([(gs > 0).sum(), (gs == 0).sum()], axis = 1)
        pd.concat([freq_per_pt, perf_per_pt], axis = 1).to_csv(os.path.join(out, "perf_per_pt.txt"), sep = "\t", encoding = "utf-8")
        res = pd.concat([micro_perf, macro_perf], axis = 1)
        res.columns = ["micro",  "macro"]
        res.to_csv(os.path.join(out, "perf_overall.txt"), sep = "\t")
        #write to disk confusion matrix
        conf_per_pt.columns = ["True negatives", "False positives", "False negatives", "True positives"] 
        conf_per_pt.to_csv(os.path.join(out, "confusion_matrix_per_pt.txt"), sep = "\t", encoding = "utf-8")
        
  
    @staticmethod
    def get_miscl(y, y_pred):
        """get misclassified samples""" 
        #print y[y == 1]
        #print y[y == 0]
        #print y_pred[y == 1]
        #print y_pred[y == 0]
        #print y[y == 1] != y_pred[y == 1]
        #print y[y == 0] != y_pred[y == 0]
        FN =  y[y == 1].loc[y[y == 1] != y_pred[y == 1],].index
        FP =  y[y == 0].loc[y[y == 0] != y_pred[y == 0], ].index
        return FN.tolist() + FP.tolist()

    @staticmethod 
    def get_performance(conf):
        """get different performance measures from a confusion matrix"""
        recall_pos = evaluate.recall_pos_conf(conf)
        recall_neg = evaluate.recall_neg_conf(conf)
        bacc = evaluate.bacc(recall_pos, recall_neg)
        prec = evaluate.precision_conf(conf)
        s = pd.Series([recall_pos, recall_neg, bacc, prec])
        s.index = ["recall pos", "recall neg", "macro accuracy", "precision"]
        return s

    @staticmethod 
    def confusion_m(y,y_pred):
        """get confusion matrix TN FP /n FN TP"""
        TP = (y[y == 1] == y_pred[y == 1]).sum()
        FN = (y[y == 1] != y_pred[y == 1]).sum()
        TN = (y[y == 0] == y_pred[y == 0]).sum()
        FP = (y[y == 0] != y_pred[y == 0]).sum()
        return pd.np.array([TN, FP, FN, TP])

    @staticmethod 
    def recall_pos(y,y_pred):
        """compute recall of the positive class"""
        return (y[y == 1] == y_pred[y==1]).sum()/float((y==+1).sum())
    
    @staticmethod 
    def recall_pos_conf(conf):
        """compute recall of the positive class"""
        TN, FP, FN, TP = conf
        if (TP + FN) == 0:
            float('nan') 
        return TP/float(TP+FN)
    
    @staticmethod 
    def recall_neg_conf(conf):
        """compute recall of the positive class"""
        TN, FP, FN, TP = conf
        if (TN + FP) == 0:
            return float('nan')
        return TN/float(TN+FP)

    @staticmethod 
    def recall_neg(y, y_pred):
        """compute recall of the negative class"""
        return (y[y == 0] == y_pred[y==0]).sum()/float((y==0).sum())
    
    @staticmethod 
    def precision(y, y_pred):
        """compute precision"""
        TP = (y[y == 1] == y_pred[y == 1]).sum()
        FP = (y[y == 0] != y_pred[y == 0]).sum()
        if (TP + FP) == 0:
            return 0
        return TP / float(TP + FP)   
    @staticmethod 
    def precision_conf(conf):
        """compute precision"""
        TN, FP, FN, TP = conf
        if (TP + FP) == 0:
            return float('nan')
        return TP / float(TP + FP)

    @staticmethod
    def bacc(pos_acc, neg_acc):
        """compute balanced accuracy"""
        return float(pos_acc + neg_acc)/2
