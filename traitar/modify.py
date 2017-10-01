import pandas as pd
import tarfile
import os
import StringIO
from . import PhenotypeCollection
import sys
import json

mfs = ["%s_bias.txt", "%s_feats.txt","%s_non-zero+weights.txt", "%s_perf.txt"]

def validate(model_dir, pts):
    """validate that there is a model for each phenotype"""
    for i in pts.index:
        for j in mfs:
            if not os.path.exists(os.path.join(model_dir, j % i)):
                sys.stderr.write("%s does not exists" % os.path.join(model_dir, j % i))
                raise Exception

def remove(archive_f, phenotypes, out_f, keep = False):
    """remove given phenotypes from the pt archive and write a new archive"""
    pts = pd.read_csv(phenotypes, index_col = 0, header = None)
    pts.index = pts.index.values.astype('string')
    ptc = PhenotypeCollection.PhenotypeCollection(archive_f)
    pt2acc = ptc.get_pt2acc()
    pt2acc.index = pt2acc.index.values.astype('string')
    if not keep:
        pt2acc.drop(pts.index, inplace = True)
    else:
        pt2acc = pt2acc.loc[pts.index,]
    pt2acc_s = StringIO.StringIO(pt2acc.to_csv(sep = "\t", encoding = 'utf8'))
    pt2acc_tarinfo = tarfile.TarInfo("pt2acc.txt")
    pt2acc_tarinfo.size = len(pt2acc_s.buf)
    t = ptc.tar
    members = [] 
    members.append((t.extractfile("pf2acc_desc.txt"), t.getmember("pf2acc_desc.txt")))
    members.append((t.extractfile("config.txt"), t.getmember("config.txt")))
    for i in pt2acc.index:
        members.append((t.extractfile("%s_non-zero+weights.txt" % i),t.getmember("%s_non-zero+weights.txt" % i)))
        members.append((t.extractfile("%s_bias.txt" % i), t.getmember("%s_bias.txt" % i)))
        members.append((t.extractfile("%s_feats.txt" % i), t.getmember("%s_feats.txt" % i)))
    out_tar = tarfile.open("%s" % out_f, "w:gz")
    out_tar.addfile(pt2acc_tarinfo, pt2acc_s)
    for m in members:
        out_tar.addfile(m[1], m[0])
    out_tar.close()

def extend(archive_f, model_dir, pf2acc_desc_f, pts, pfam_v): 
    """extend an existing model archive by additional models"""
    #validate new phenotype models 
    validate(model_dir, pt2acc_desc_f)
    ptc = PhenotypeCollection.PhenotypeCollection(archive_f)
    pt2acc = ptc.get_pt2acc()
    #check if phenotype models not already exist
    #add phenotype models to the archive

def new(models_dir, pf2acc_f, pt2desc_f, hmm_name,  archive_name, is_standardized, pt_table, tree):
    """create new archive with phenotype models"""
    #read in pf and pt accessions 
    if pt2desc_f is not None:
        pts = pd.read_csv(pt2desc_f, sep = "\t", index_col = 0)
        #check if models for all phenotypes exist
        validate(models_dir, pts)
    create_tar(models_dir, pf2acc_f, pt2desc_f, hmm_name,  archive_name, is_standardized, pt_table, tree)

def create_tar(models_dir, pf2acc_f, pt2desc_f, hmm_name,  archive_name, is_standardized, pt_table, tree):
    #create tar archive
    config = {}
    t = tarfile.open("%s.tar.gz" % archive_name, "w:gz")
    t.add(pf2acc_f, arcname = "pf2acc_desc.txt")
    if not pt2desc_f is None:
        pt2desc = pd.read_csv(pt2desc_f, sep = "\t", index_col = 0) 
        t.add(pt2desc_f, arcname = "pt2acc.txt")
        for i in pt2desc.index:
            #check if standardization parameters are available and add to model archive
            if is_standardized:
                mfs.append("%s_mean.txt")
                mfs.append("%s_scale.txt")
            for j in mfs:
                t.add(os.path.join(models_dir, j % i), arcname = os.path.basename(os.path.join(models_dir, j % i)))
    if not pt_table is None:
        t.add(pt_table, arcname = "phenotype_table.txt")
    if not models_dir is None:
        cv_acc_f = os.path.join(models_dir, "cv_acc.txt")
        if os.path.exists(cv_acc_f): 
            t.add(os.path.join(models_dir, "cv_acc.txt"), arcname = "cv_acc.txt")
        if tree is not None:
            t.add(pt_table, arcname = "tree.txt")
        with open(os.path.join(models_dir, "config.json")) as jf:    
                config = json.load(jf)
    config["archive_name"] = archive_name.split("/")[-1]
    config["annot_name"]  = hmm_name
    if is_standardized:
        config["is_standardized"] = "True"
    else:
        config["is_standardized"] = "False"
    config_s = StringIO.StringIO(json.dumps(config, indent=4, separators=(',', ': ')))
    config_tarinfo = tarfile.TarInfo("config.txt")
    config_tarinfo.size = len(config_s.buf)
    t.addfile(config_tarinfo, config_s)
    t.close()
         
    
    
    
    
