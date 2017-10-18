import tarfile
import zipfile
import pandas as pd
import re
import os
import os.path
from .traitar import phenolyze
from shutil import copyfile


def get_sample_names(namelist):
    """parse sample names"""
    sufs = [".fna", ".faa", ".fasta"]
    out_sample_file_names = []
    out_sample_names = []
    for f in namelist:
        #replace any non standard characters with underlines
        repl = re.sub("[^A-Za-z0-9.-]", "_", f)
        out_sample_file_names.append(repl)
        #check if there is a file ending that needs to be replace for the sample name
        sn = f
        for suf in sufs:
            if repl.endswith(suf):
                sn = repl.replace(suf, "")
            else:
                sn = repl
        out_sample_names.append(sn)
    return out_sample_file_names, out_sample_names
    
    

def read_archive(input_archive, archive_type, mode, sample2cat, input_dir, input_names):
    """read archive"""
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)

    if archive_type == "zip" or archive_type == "tar.gz":
        if archive_type == "zip":
            archive = zipfile.open(input_archive)
            namelist = archive.namelist()
        if archive_type == "tar.gz":
            archive = tarfile.open(input_archive, "r")
            namelist = archive.getnames()
        sample_file_names, sample_names = get_sample_names(namelist)
        for tf, sfn in zip(namelist, sample_file_names):
                extracted = archive.extractfile(tf) 
                with open("%s/%s" % (input_dir, sfn), 'w') as sample_file_out:
                    for line in extracted:
                        sample_file_out.write(line) 
                extracted.close()
    elif archive_type == "directory":
        sample_names = input_names.split(',')
        sample_file_names = []
        for input_part in input_archive.split(','):
            input_dir_part=os.path.basename(input_part)
            sample_file_names.append(input_dir_part)
            os.symlink(input_part, input_dir+"/"+input_dir_part)

            
    #create sample table
    if sample2cat is not None:
        sample_cat = pd.read_csv(sample2cat, index_col = 0, sep = "\t")
        #replace index with cleaned file names
        if archive_type != "directory":
            sample_cat.index.rename(str, dict([(tf, sfn) for sfn, tf in zip(sample_file_names, namelist)]))
            sample_table = pd.DataFrame(sample_names)
            categories = pd.Series(sample_cat.loc[sample_file_names, ]['category'].tolist())
        else:
            sample_table = pd.DataFrame(sample_file_names)
            categories = pd.Series(sample_cat.loc[sample_names, ]['category'].tolist())
        sample_table['category'] = categories          
        sample_table.columns = ["sample_file_name", "category"]
    else:
        sample_table = pd.DataFrame(sample_file_names)
        sample_table.columns = ["sample_file_name"]
    sample_table.index = sample_names
    sample_table.index.name = "sample_name"
    sample_table.to_csv("%s/sample_table.txt" % input_dir, sep = "\t")  
    
         

def call_traitar(args):
    args.rearrange_heatmap = args.no_heatmap_phenotype_clustering = args.no_heatmap_sample_clustering = args.gene_gff_type = args.primary_models = args.secondary_models = None
    args.sample2file = "%s/sample_table.txt" % args.input_dir 
    phenolyze(args)
    #compress output
    
    if args.generate_galaxy_html is not None:
        (html_file, html_dir) = args.generate_galaxy_html. split(':')
        os.makedirs(html_dir)
        image_name = args.output_dir+"/phenotype_prediction/heatmap_combined.%s" % args.heatmap_format
        copyfile(image_name, html_dir)
        with tarfile.open(html_dir+"/archive.tar.gz", "w:gz") as tar:
            tar.add(args.output_dir, arcname=os.path.basename(args.output_dir))
        copyfile('html/sample.html', html_file)
    else:
        with tarfile.open(args.out_archive, "w:gz") as tar:
            tar.add(args.output_dir, arcname=os.path.basename(args.output_dir))

        if args.output_image is not None:
            image_source = args.output_dir+"/phenotype_prediction/heatmap_combined.%s" % args.heatmap_format
            if args.output_image[0:1] == '/':
                output_image = args.output_image
            else:
                output_image = os.path.dirname(args.out_archive)+'/'+args.output_image
        
            copyfile(image_source, output_image)
