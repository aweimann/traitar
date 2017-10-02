#!/usr/bin/env python
import pandas as ps
import os
from . import modify
import subprocess
import sys
import shutil
from . import get_external_data
from . import __file__ 
from . import hmm2gff
import json
import os
import sys
import tarfile
from .PhenotypeCollection import PhenotypeCollection
from ._version import __version__
import re
from . import evaluation
env = os.environ.copy()
primary_default_models = "%(phenolyzer)s/data/models/phypat.tar.gz" %{"phenolyzer" : os.path.abspath(os.path.dirname(__file__))}
secondary_default_models = "%(phenolyzer)s/data/models/phypat+PGL.tar.gz" %{"phenolyzer" : os.path.abspath(os.path.dirname(__file__))}

def annotate(args):
    """annotate input"""
    p = Traitar(args.input_dir, args.output_dir, args.sample2file, args.cpus, gene_gff_type=args.gene_gff_type, primary_models = args.primary_models, secondary_models = args.secondary_models, primary_hmm_db = args.primary_hmm_db, secondary_hmm_db = args.secondary_hmm_db, annotation_summary = args.annotation_summary)
    p.annotate(args.mode)

def evaluate(args):
    evaluation.evaluate.evaluate(args.out, args.gold_standard_f, args.traitar_pred_f, args.min_samples, args.are_pt_ids, args.phenotype_archive)

def phenolyze(args):
    """annotate and then run phenotyping"""
    p = Traitar(args.input_dir, args.output_dir, args.sample2file, args.cpus, args.rearrange_heatmap, args.heatmap_format, args.no_heatmap_phenotype_clustering, args.no_heatmap_sample_clustering, args.gene_gff_type, args.primary_models, args.secondary_models, args.primary_hmm_db, args.secondary_hmm_db, args.annotation_summary)
    #check if user wants to supply pre-computed annotation summary
    if not args.mode == "from_annotation_summary":
        p.annotate(args.mode)
    else:
        if not "annotation_summary" in args:
            sys.stderr.write("Provide an annotation file with -a / --annotation_summary along with from_annotation_summary")
            sys.exit(1)
    p.phenolyze(args.mode)

def new(args):
    """create new phenotype model archive"""
    modify.new(args.models_dir, args.feat2desc, args.pt2desc, args.annot_name,  args.archive_name, args.is_standardized, args.pt_table, args.tree)

def show(args):
    """show features for the given phenotype"""
    if args.models_f is not None:
        pc = [PhenotypeCollection(args.models_f)]
    else:
        if args.predictor == "phypat":
            pc = [PhenotypeCollection(primary_default_models)]
        else:
            pc = [PhenotypeCollection(secondary_default_models)]
    for i in pc:
        i.get_selected_features(args.phenotype, args.strategy, args.include_negative).to_csv(sys.stdout, sep = "\t", float_format='%.3f')

def remove(args):
    """remove phenotypes from phenotype tar archive"""
    modify.remove(args.archive_f, args.phenotypes, args.out_f, args.keep)

def extend(args):
    pass

class Traitar:

    def __init__(self, input_dir, output_dir, sample2file, cpu = 1, heatmap_out = None, heatmap_format = "pdf", no_heatmap_phenotype_clustering = False, no_heatmap_sample_clustering = False ,gene_gff_type = None, primary_models = None, secondary_models = None, primary_hmm_db = None, secondary_hmm_db = None, annotation_summary = None):
        self.user_message = "output dir %s already exists; press 1 to continue with data from a previous run; press 2 to remove this directory; press 3 to abort followed by [ENTER]"
        self.error_message =  "directory %s already exists; delete directory or run in interactive mode if this step is done and you want to continue from there"
        self.heatmap_format = heatmap_format 
        self.no_heatmap_phenotype_clustering = no_heatmap_phenotype_clustering
        self.no_heatmap_sample_clustering = no_heatmap_sample_clustering
        self.sample2file = sample2file
        self.input_dir = input_dir
        self.gene_gff_type = gene_gff_type
        self.s2f = self.parse_sample_f()
        self.cpu = cpu
        self.output_dir = output_dir
        self.phenolyzer_dir = os.path.abspath(os.path.dirname(__file__)) 
        self.annotation_summary = annotation_summary
        #load  configuration file
        with open(os.path.join(self.phenolyzer_dir, "config.json" ), "r") as cf:
            self.config = json.load(cf)
        #set primary models either to user specified model or to default phypat model
        if primary_models is not None:
            self.primary_models = PhenotypeCollection(primary_models)
            self.primary_models.hmm_f = primary_hmm_db
        else:
            self.primary_models = PhenotypeCollection(primary_default_models)
            self.primary_models.hmm_f =  config["hmm_f"]
        #set secondary models either to user specified model or to default phypat model
        if secondary_models is not None and primary_models is not None:
            self.secondary_models = PhenotypeCollection(secondary_models)
            self.secondary_models.hmm_f = secondary_hmm_db
        elif self.primary_models.get_name() != "phypat":
            self.secondary_models = None
        else:
            self.secondary_models = PhenotypeCollection(secondary_default_models)
            self.secondary_models.hmm_f =  config["hmm_f"]
        self.is_gnu_parallel_available = self.is_exe("parallel")
        self.heatmap_out = heatmap_out
        #check if GNU parallel is available when running with more than one cpu
        if cpu != 1 and not self.is_gnu_parallel_available:
            sys.stderr.write("GNU parallel is not available on the command line; make sure you installed it properly or decrease number of cpus to 1\n")
            sys.exit(1)
        #check if config exists
        if not os.path.isfile(os.path.join(self.phenolyzer_dir, "config.json")):
            sys.stderr.write("config.json does not exists; make sure that you have run traitar pfam")
            sys.exit(1)
            #check if Pfam hmm specified in config.json exists
        #create output dir
        self.check_dir(output_dir)
        #pred dir
        self.pred_dir = os.path.join(self.output_dir, "phenotype_prediction")
        self.phypat_dir = os.path.join(self.pred_dir, self.primary_models.get_name())
        if not self.secondary_models is None:
            self.phypat_pgl_dir = os.path.join(self.pred_dir, self.secondary_models.get_name())

    def _special_match(self, strg, search = re.compile(r'[^A-Za-z0-9.\-_]').search):
        return not bool(search(strg))
   
    def parse_sample_f(self):
        """read sample file with pandas and make sure the content is appropriate"""
        #check if sample files exist
        if not os.path.exists(self.sample2file):
            sys.exit("sample file %s does not exist" % self.sample2file)
        s2f = ps.read_csv(self.sample2file, dtype = 'string', sep = "\t")
        col_check = dict((i, False) for i in ["sample_file_name", "sample_name", "category", "gene_gff"])
        for i in s2f.columns:
            if i not in ["sample_file_name", "sample_name", "category", "gene_gff"]:
                sys.exit("%s is not a valid column identifier" % i)
            col_check[i] = True
        if not col_check["sample_file_name"]:
            sys.exit("sample_file_name column in input sample_file missing")
        if not col_check["sample_name"]:
            sys.exit("sample_name colun in input sample_file missing")
        for i in s2f.loc[:, "sample_file_name"]:
            if not os.path.exists(os.path.join(self.input_dir,i)):
                #sys.exit("sample file %s does not exist in the output directory %s" % (self.input_dir, i))
                pass
        for i in s2f.loc[:, "sample_name"]:
            if not self._special_match(i):
                sys.exit("invalid character in sample name %s; only [a-zA-Z0-9.-_] allowed" % i)
            if len(i) > 41:
                sys.exit("sample names may only have 40 characters %s" %i)
        if col_check["category"]:
            uq = s2f.loc[:, "category"].unique()
            if len(uq) > 12:
                sys.exit("reduce the number of sample categories to less than 15")
            for i in uq:
                if len(i) > 30:
                    sys.exit("sample categories may not be longer than 30 characters %s" %i)
        if col_check["gene_gff"]:
            for i in s2f.loc[:, "gene_gff"]:
                if not os.path.exists(os.path.join(self.input_dir,i)):
                    sys.stderr.write("WARNING: sample file %s does not exist in the output directory %s; program will fail if not in 'from_summary_annotation' mode" % (self.input_dir, i))
            if self.gene_gff_type is None:
                sys.exit("gene gff type needs to be specified with -g / --gene_gff_type <gene_gff_type> if sample file contains gene_gff column")
        return s2f

    def is_exe(self, program):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return True 
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return True 
                else:
                    False 

    def check_dir(self, out_dir):
        if os.path.exists(out_dir):
            #check if the process is run in background
            try: 
                if os.getpgrp() == os.tcgetpgrp(sys.stdout.fileno()):
                    print  self.user_message % out_dir
                    user_input = raw_input()
                    while user_input not in ["1", "2", "3"]:
                        user_input = raw_input().strip()
                    if user_input == "1":
                        return False 
                    if user_input == "2":
                        shutil.rmtree(out_dir)
                        os.mkdir(out_dir)
                        return True 
                    if user_input == "3":
                        sys.exit(1)
                else:
                    sys.stderr.write(self.error_message % out_dir)
            except OSError:
                    sys.stderr.write(self.error_message % out_dir)
        else: 
            os.mkdir(out_dir)
            return True

    def annotate(self, mode):
        pfam_msg = "running annotation with hmmer. This step can take a while. A rough estimate for sequential Pfam annotation of genome samples of ~3 Mbs is 10 min per genome."
        #check if executables are available 
        if not self.is_exe("hmmsearch"):
            sys.stderr.write("hmmsearch not available on command line; please make sure you have properly installed it\n")
            sys.exit(1)
        if mode == "from_nucleotides":
            if not self.is_exe("prodigal"):
                sys.stderr.write("prodigal not available on command line; please make sure you have properly installed it\n")
                sys.exit(1)
            print "running gene prediction with Prodigal"
            sys.stdout.flush()
            self.run_gene_prediction(self.s2f.loc[:,"sample_file_name"], self.s2f.loc[:,"sample_name"])
            print pfam_msg
            sys.stdout.flush()
            self.run_hmmer_annotation(self.s2f.loc[:,"sample_name"], self.s2f.loc[:,"sample_name"], mode)
        if mode == "from_genes": 
            print pfam_msg
            sys.stdout.flush()
            self.run_hmmer_annotation(self.s2f.loc[:,"sample_file_name"], self.s2f.loc[:,"sample_name"], mode)

    def phenolyze(self, mode):
        print "running phenotype prediction"
        sys.stdout.flush()
        is_recompute = self.check_dir(self.pred_dir)
        if is_recompute:
            self.run_phenotype_prediction()
        print "running feature track generation"
        if not mode == "from_annotation_summary":
            if not "gene_gff" in self.s2f.columns:
                self.run_feature_track_generation(self.s2f.loc[:,"sample_name"], mode)
            else:
                self.run_feature_track_generation(self.s2f.loc[:, "sample_name"],  mode, self.s2f.loc[: , "gene_gff"], self.gene_gff_type)
        sys.stdout.flush()
        print "running heatmap generation"
        self.run_heatmap_generation(self.s2f.loc[:,"sample_name"])
        sys.stdout.flush()
    
    def execute_commands(self, commands, joblog = None):
        devnull = open('/dev/null', 'w')
        from subprocess import Popen, PIPE
        if len(commands) == 0:
            #nothing to do
            return
        if self.cpu > 1:
            #run with parallel
            #ps.DataFrame(commands).to_csv(tf, index = False, header = False) 
            p = Popen("parallel --will-cite %s -j %s" %  ("--joblog %s" % joblog if joblog is not None else "", self.cpu),  stdout = devnull, shell = True,  executable = "/bin/bash", stdin = PIPE, env = env)
            p.communicate(input = "\n".join(commands))
            if p.returncode != 0:
                if not joblog is None:
                    sys.stderr.write("Non-zero exit value; commands can be found in %s\n" % joblog) 
                else:
                    sys.stderr.write("Non-zero exit value; command(s) %s failed\n" % "\n".join(commands)) 
                sys.exit()
        else:
            #run in sequential order
            for i in commands:
                p = Popen(i,  stdout = devnull, shell = True,  executable = "/bin/bash", stdin = PIPE, env = env)
                p.communicate(input = i)
                if p.returncode != 0:
                    sys.stderr.write("Non-zero exit value; command %s failed\n" % i) 
                    sys.exit()

    def run_gene_prediction(self, in_samples, out_samples):
        #create output directory for the gene prediction 
        gp_dir = os.path.join(self.output_dir, "gene_prediction")
        is_recompute = self.check_dir(gp_dir)
        prodigal = "prodigal < \"%(in_dir)s/%(in_sample)s\" > %(gp_dir)s/%(out_sample)s.gff  -a %(gp_dir)s/%(out_sample)s.faa  -f gff"
        prodigal_commands = []
        for i in range(len(in_samples)):
            prodigal_commands.append(prodigal % {"in_dir": self.input_dir, "in_sample": in_samples[i], "out_sample":out_samples[i], "gp_dir":gp_dir})
        if is_recompute:    
            self.execute_commands(prodigal_commands, os.path.join(gp_dir, "joblog.txt"))

    def run_hmmer_annotation(self, in_samples, out_samples, mode):
        param_dict = {
                "out_dir" : self.output_dir, 
                "phenolyzer_dir" : self.phenolyzer_dir}
        if mode == "from_nucleotides":
            param_dict["file_extension"] = ".faa" 
            param_dict["in_dir"] = os.path.join(self.output_dir, "gene_prediction")
        else:
            param_dict["file_extension"] = "" 
            param_dict["in_dir"] = self.input_dir
        #create output directory for the pfam annotation 
        a_dir_base = os.path.join(self.output_dir, "annotation")
        #check if output directory already exists and trigger user input if in interactive mode
        is_recompute = self.check_dir(a_dir_base)
        #run hmmer annotation
        hmmer = "hmmsearch --cpu 1 --cut_ga  --domtblout %(a_dir)s/%(out_sample)s_domtblout.dat  %(hmm_f)s > /dev/null \"%(in_dir)s/%(in_sample)s%(file_extension)s\""
        filter_and_aggregate = "hmmer2filtered_best %(a_dir)s/%(out_sample)s_domtblout.dat   %(a_dir)s/%(out_sample)s_filtered_best.dat %(hmm_name)s"

        if self.secondary_models is not None and self.primary_models.get_hmm_name() != self.secondary_models.get_hmm_name():
            models = [self.primary_models, self.secondary_models]
        else:
            models = [self.primary_models]
        for pt_models in models:
            a_dir = os.path.join(a_dir_base, pt_models.get_hmm_name())
            param_dict["a_dir"] = a_dir
            #param_dict["hmms"] =  self.config["hmms"] 
            param_dict["hmm_f"] =  pt_models.get_hmm_f()
            param_dict["hmm_name"] =  pt_models.get_hmm_name()
            param_dict["archive_f"] =  pt_models.get_archive_f()
            if not os.path.exists(a_dir):
                os.mkdir(a_dir)
            hmmer_commands = []
            fae_commands = []
            for i in range(len(in_samples)):
                fname_hmm = os.path.join(a_dir, out_samples[i] + "_domtblout.dat")
                fname_filtered = os.path.join(a_dir, out_samples[i] + "_filtered_best.dat")
                param_dict["in_sample"] = in_samples[i]
                param_dict["out_sample"] = out_samples[i]
                if not os.path.isfile(fname_hmm):
                    hmmer_commands.append(hmmer % param_dict) 
                    if is_recompute:
                        print "hmmer output %s is missing; recomputing ..." %fname_hmm  
                #run filtering and best domain hit aggregation
                if not os.path.isfile(fname_filtered):
                    if is_recompute:
                        print "filtered hmmer output %s is missing; recomputing ..." %fname_filtered
                    fae_commands.append(filter_and_aggregate % param_dict)
            self.execute_commands(hmmer_commands, joblog = os.path.join(a_dir, "joblog_hmmer.txt"))
            self.execute_commands(fae_commands, joblog = os.path.join(a_dir, "joblog_filter_and_aggregate.txt"))
            if not os.path.exists(os.path.join(a_dir, "summary.dat")):
                if is_recompute:
                    print "annotation summary matrix is missing; recomputing ..." 
                #run summary matrix computation
                domtblout2gene_generic = "domtblout2gene_generic %(a_dir)s/summary.dat  <(ls %(a_dir)s/*_filtered_best.dat) %(archive_f)s" % param_dict
                self.execute_commands([domtblout2gene_generic], joblog = os.path.join(a_dir, "joblog_generate_annotation_matrix.txt"))


    def run_phenotype_prediction(self):
        #create output directory for the phenotype prediction 
        if not os.path.exists(self.phypat_dir):
            os.mkdir(self.phypat_dir)
        #run phenotype prediction for primary and secondary models 
        param_dict = {
                "out_dir" : self.output_dir, 
                "pred_dir" : self.phypat_dir, 
                "primary_models" : self.primary_models.get_archive_f(), 
                "hmm_f" : self.primary_models.get_hmm_name(),
                "phenolyzer_dir" : self.phenolyzer_dir}
        if self.annotation_summary is None:
            param_dict["annotation_summary"] = "%(out_dir)s/annotation/%(hmm_f)s/summary.dat" % param_dict
        else: 
            param_dict["annotation_summary"] = self.annotation_summary 
        predict_phypat = "predict %(primary_models)s %(pred_dir)s  %(annotation_summary)s  -k 5" % param_dict 
        #check if secondary models are available
        if not self.secondary_models is None:
            if not os.path.exists(self.phypat_pgl_dir):
                os.mkdir(self.phypat_pgl_dir)
            param_dict["secondary_models"] = self.secondary_models.get_archive_f()
            param_dict["pred_dir"] = self.phypat_pgl_dir
            predict_phypat_pgl = "predict %(secondary_models)s %(pred_dir)s %(annotation_summary)s -k 5" % param_dict 
            param_dict["primary_name"] = self.primary_models.get_name()
            param_dict["secondary_name"] = self.secondary_models.get_name()
            param_dict["phypat_dir"] = self.phypat_dir
            param_dict["phypat_pgl_dir"] = self.phypat_pgl_dir 
            #combine phypat and phypat+PGL predictions
            merge_preds = "merge_preds %(out_dir)s/phenotype_prediction %(phypat_dir)s %(phypat_pgl_dir)s %(primary_name)s %(secondary_name)s -k 5" % param_dict
            self.execute_commands([predict_phypat, predict_phypat_pgl])
            self.execute_commands([merge_preds])
        else:
            self.execute_commands([predict_phypat])
    
    def run_feature_track_generation(self, in_samples, mode, gene_gffs = None, gene_gff_type = "prodigal"):
        """map the phenotype relevant protein families and write mapping to disk"""
        #create output directory for the pfam annotation 
        #hmm2gff command for the full pfam annotation
        hmm2gff = "hmm2gff %(out_dir)s/annotation/%(hmm_name)s/%(sample)s_filtered_best.dat  %(out_gff_dir)s %(sample)s  %(model_tar)s %(predicted_pts)s " 
        gene_gff_extd = "--gene_gff %(gene_gff)s --gene_gff_type " + gene_gff_type 
        #read in phypat predictions
        phypat_preds = ps.read_csv(os.path.join(self.phypat_dir, "predictions_majority-vote.txt"), index_col = 0, sep = "\t", encoding = "utf-8")
        phypat_preds.index = phypat_preds.index.values.astype('string')
        #read pfam phenotype id mapping file from model tar
        pt2desc_phypat = self.primary_models.get_acc2pt()
        param_dict = {
                "out_dir" : self.output_dir, 
                "pred_dir" : self.phypat_dir, 
                "phenolyzer" : self.phenolyzer_dir}
        #secondary models
        if not self.secondary_models is None:
            phypat_pgl_preds = ps.read_csv(os.path.join(self.phypat_pgl_dir, "predictions_majority-vote.txt"), index_col = 0, sep = "\t", encoding = "utf-8") 
            phypat_pgl_preds.index = phypat_preds.index.values.astype('string')
            pt2desc_phypat_pgl = self.secondary_models.get_acc2pt() 
        #collect predictions and compose hmm2gff command for each samples
        h2gff_commands = []
        for i in range(len(in_samples)):
            param_dict["sample"] = in_samples[i]
            param_dict["model_tar"] = self.primary_models.get_archive_f()
            param_dict["out_gff_dir"] = "%s/feat_gffs/" % self.phypat_dir
            param_dict["hmm_name"] = self.primary_models.get_hmm_name() 
            predicted_pts_phypat = []
            predicted_pts_phypat_pgl = []
            for j in phypat_preds.columns:
                if phypat_preds.loc[in_samples[i], j] == 1:
                    predicted_pts_phypat.append(str(pt2desc_phypat.loc[j,][0]))
            if not self.secondary_models is None:
                for j in phypat_pgl_preds.columns:
                    if not self.secondary_models is None and phypat_pgl_preds.loc[in_samples[i], j] == 1:
                        predicted_pts_phypat_pgl.append(str(pt2desc_phypat_pgl.loc[j,][0]))
            if not len(predicted_pts_phypat) == 0:
                param_dict["predicted_pts"] =  ",".join(predicted_pts_phypat)
                cmd = hmm2gff % param_dict 
                if gene_gffs is not None:
                    cmd += gene_gff_extd % {"gene_gff" : ("%s/gene_prediction/%s.gff" % (self.output_dir, in_samples[i])) if gene_gffs is None else os.path.join(self.input_dir, gene_gffs[i])}
                h2gff_commands.append(cmd)
            if not len(predicted_pts_phypat_pgl) == 0:
                param_dict["model_tar"] = self.secondary_models.get_archive_f()
                param_dict["out_gff_dir"] = "%s/feat_gffs/" % self.phypat_pgl_dir
                param_dict["predicted_pts"] = ",".join(predicted_pts_phypat_pgl)
                param_dict["hmm_name"] = self.secondary_models.get_hmm_name() 
                if not self.secondary_models is None:
                    cmd = hmm2gff % param_dict 
                    if gene_gffs is not None:
                        cmd += gene_gff_extd % {"gene_gff" : ("%s/gene_prediction/%s.gff" % (self.output_dir, in_samples[i])) if gene_gffs is None else os.path.join(self.input_dir, gene_gffs[i])}
                    h2gff_commands.append(cmd)
        #create output dirs for feature tracks
        is_recompute = self.check_dir(os.path.join(self.phypat_dir, "feat_gffs"))
        if is_recompute:
            if not self.secondary_models is None:
                os.mkdir(os.path.join(self.phypat_pgl_dir, "feat_gffs"))
            self.execute_commands(h2gff_commands) 
        if mode != "from_nucleotides" and gene_gffs is None:
            ftco = os.path.join(self.pred_dir, "feature_track_commands.txt")
            sys.stderr.write("tracks with Pfams relevant for the predictions cannot be ad-hoc generated because the input is amino acids and no gene prediction GFF files have been generated\n commands are saved to %s and can be modified and run manually\n" % ftco)
            with open(ftco, 'w') as f:
                f.write("\n".join(h2gff_commands))

    def run_heatmap_generation(self, in_samples, compatible = False):
        """generate a heatmap from the results"""
        #TODO introduce nans via the merge prediction routine; mark the nans in the heatmap with a different color
        if self.heatmap_out is None:
            self.heatmap_out = self.pred_dir
        hm_cmd = "heatmap %(pred_dir)s/%(pred_f)s %(out)s/heatmap_%(predictor)s.%(heatmap_format)s %(secondary_models)s --sample_f %(sample_file)s %(model_archive)s %(phenolyzer)s/data/colors.txt %(phenotype_clustering)s %(sample_clustering)s"
        hm_dict = {"phenolyzer" : self.phenolyzer_dir,  "out": self.heatmap_out, "sample_file" : self.sample2file  , "heatmap_format" : self.heatmap_format, "phenotype_clustering" : "--column_method None" if self.no_heatmap_phenotype_clustering else "", "sample_clustering" : "--row_method None" if self.no_heatmap_sample_clustering else "" }
        hm_dict_phypat = {"pred_dir" : self.phypat_dir, "pred_f" : "predictions_majority-vote.txt","predictor": self.primary_models.get_name(),  "model_archive" : self.primary_models.get_archive_f(), "secondary_models" : ""}
        cmds = []
        if self.secondary_models is None:
            hm_dict_phypat.update(hm_dict)
            cmds.append(hm_cmd % hm_dict_phypat)
        else:
            hm_dict_phypat_pgl = {"pred_dir" : self.phypat_pgl_dir, "pred_f" : "predictions_majority-vote.txt", "predictor": self.secondary_models.get_name(),"model_archive" : self.secondary_models.get_archive_f(), "secondary_models" : ""}
            hm_dict_comb = {"pred_dir" : self.pred_dir ,"pred_f" : "predictions_majority-vote_combined.txt","predictor": "combined", "model_archive" : self.primary_models.get_archive_f(), "secondary_models" : "--secondary_model_tar %s" % self.secondary_models.get_archive_f()}
            for i in [hm_dict_phypat, hm_dict_phypat_pgl, hm_dict_comb]:
                i.update(hm_dict)
                cmds.append(hm_cmd % i)
        self.execute_commands(cmds)
