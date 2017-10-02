import pandas as pd
import tarfile
import sys
import json

class PhenotypeCollection:

    def __init__(self, archive_f):
        self.archive_f = archive_f
        self.tar = tarfile.open(archive_f, mode = "r:gz")
        info = self.tar.extractfile("config.txt")
        cf = json.load(info)
        self.name = cf["archive_name"]
        self.hmm_name = cf["annot_name"]
        if "hmm_f" in cf:
            self.hmm_f = cf.loc["hmm_f"]
        if not "is_standardized" in cf:
            self.is_standardized = False 
        else:
            #otherwise assume the feature used to train the model were not standardized
            self.is_standardized = True if cf["is_standardized"] == "True" else False 

    def get_pt2acc(self):
        """get and parse phenotype to phenotype id to accession mapping"""
        pt2acc_f = self.tar.extractfile("pt2acc.txt")
        pt2acc = pd.read_csv(pt2acc_f, sep = "\t", index_col = 0,  encoding='utf-8')
        pt2acc.index = pt2acc.index.astype('U')
        return pt2acc
    
    def get_scale_and_mean(self, pt):
        """get and parse scale and mean values used for standardizing the original data"""
        scale_f = self.tar.extractfile("%s_scale.txt" % pt)
        scale = pd.read_csv(scale_f, sep = "\t", index_col = 0)
        scale.index = scale.index.astype('U')
        mean_f = self.tar.extractfile("%s_mean.txt" % pt)
        mean = pd.read_csv(mean_f, sep = "\t", index_col = 0,  encoding='utf-8')
        mean.index = scale.index.astype('U')
        return scale, mean 

    def get_pt2id(self):
        """get and parse phenotype to phenotype accession to id mapping"""
        pt2id_f = self.tar.extractfile("pt2acc.txt")
        pt2id = pd.read_csv(pt2id_f, sep = "\t", index_col = 1,  encoding='utf-8')
        pt2id.index = pt2id.index.astype('U')
        return pt2id

    def get_acc2pt(self):
        """get and parse phenotype to phenotype accession to id mapping"""
        pt2acc_f = self.tar.extractfile("pt2acc.txt")
        pt2acc = pd.read_csv(pt2acc_f, sep = "\t", index_col = 1,  encoding='utf-8')
        pt2acc.index = pt2acc.index.astype('U')
        return pt2acc

    def get_pf2desc(self):
        """get and parse pfam (or any feature) to description mapping"""
        pfam_mapping = pd.read_csv(self.tar.extractfile("pf2acc_desc.txt"), index_col = 0, sep = "\t")
        return pfam_mapping

    def get_name(self):
        return self.name

    def get_archive_f(self):
        return self.archive_f

    def get_hmm_f(self):
        return self.hmm_f

    def get_hmm_name(self):
        return self.hmm_name

    def get_bias(self, pt):
        """get SVM bias term"""
        bias_f = self.tar.extractfile("%s_bias.txt"%(pt))
        bias = pd.read_csv(bias_f, sep = "\t", index_col = 0, header = None)
        return bias

    def get_predictors(self, pt):
        """get SVM weight matrix"""
        extracted_f = self.tar.extractfile("%s_feats.txt"%(pt))
        predictors = pd.read_csv(extracted_f, sep = "\t",  index_col = 0 )
        return predictors

    def get_selected_features(self, pt, strategy, include_negative):
        """get all or subset of non-zero features in the model"""
        pt2id = self.get_pt2id()
        try:
            pt_id = pt2id.loc[pt,].iloc[0]
        except:
            pt_id = pt
        try:
            extracted_f = self.tar.extractfile("%s_non-zero+weights.txt" % pt_id)
        except KeyError:
            sys.stderr.write("target phenotype %s has no associated model in the phenotype collection\n" % pt)
            sys.exit(1)
        feats = pd.read_csv(extracted_f, sep = "\t", index_col = 0)
        #use the 5 best models
        feats = feats.loc[:, feats.columns[0:6].tolist() + feats.columns[-2:].tolist()]
        pos = feats.apply(lambda x: (x.iloc[1:5] > 0).sum(), axis = 1)
        neg = feats.apply(lambda x: (x.iloc[1:5] < 0).sum(), axis = 1)
        if strategy == "non-zero":
            if include_negative:
                out_feats = feats.loc[(pos > 0) | (neg > 0),]
            else: 
                out_feats = feats.loc[(pos > 0),]
        else:
            if include_negative:
                out_feats = feats.loc[((pos > 3) | (neg > 3) ), ]
            else:
                out_feats =  feats.loc[pos > 3, ]
        return out_feats
          

