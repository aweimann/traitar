import sys
import pandas as ps
import tarfile
import StringIO
import warnings
from .PhenotypeCollection import PhenotypeCollection
warnings.filterwarnings("ignore", category=FutureWarning) 
#ignore the following warning; this piece of code needs to be adjusted in a future version of traitar
#/home/aaron/traitar/traitar/hmm2gff.py:226: FutureWarning: sort(columns=....) is deprecated, use sort_values(by=.....)
#  out_table_df.sort(columns = ["Phenotype", "cor"], ascending = [True, False]).to_csv("%s/%s.dat" % (out_gff_dir, sample), sep = "\t")
""" produce gff file from a hmmer hit file """

def get_protein_acc(fasta_id):
    """extract from a fasta id like gi|17233466|ref|NP_490586.1| the protein accession NP_490586.1"""
    return fasta_id.split("|")[3]
    #return "_".join(fasta_id.split("|")[1].split("_")[3:5])


def read_rel_feats(pt_models, pts):
    """read in pfams relevant for all predicted phenotypes"""
    #pfam to phenotype mapping
    pf2pt = {}
    #store relevant pfam lists
    rel_feat_list = []
    #extract pfam lists
    for pt in pts:
        rel_feat_list.append(pt_models.get_selected_features(pt, strategy = "majority", include_negative = False))
    #process pfam lists
    for pfs, pt in zip(rel_feat_list, pts): 
        for pf in pfs.index:
            if not pf in pf2pt:
                if not type(pfs.index.values) == type(ps.Series()):
                    pf2pt[pf] = [(pt, pfs.loc[pf, "cor"])]
                else: 
                    print pf, pt
                    print pfs.loc[pf, "Pfam_acc"]
            else:
                if not pt in pf2pt[pf]:
                    pf2pt[pf].append((pt, pfs.loc[pf, "cor"]))
    return pf2pt


def read_gff(gff_file, mode):
    """read a gene calling gff"""
    f = open(gff_file)
    gene_dict = {}
    i = 0
    fasta = False
    for l in f:
        i += 1
        # print i
        if l.startswith("#") or l.startswith("\n"):
            continue
        else:
            if mode == "metagenemark":
                read_genemark_entry(l, gene_dict)
            elif mode == "prodigal":
                read_prodigal_entry(l, gene_dict)
            elif mode == "ncbi":
                read_ncbi_entry(l, gene_dict)
            elif mode == "refseq":
                read_refseq_entry(l, gene_dict)
            elif mode == "img":
                read_img_entry(l, gene_dict)
            elif mode == "genbank":
                #in genbank gene gffs the original nucleotide sequence can be appended to the gff"
                if l.startswith(">"):
                    fasta = True 
                    continue
                if fasta and l.startswith(('a','g', 'c', 't')):
                    continue
                else:
                    fasta = False
                read_genbank_entry(l, gene_dict)
    f.close()
    return gene_dict


def read_genemark_entry(l, gene_dict):
    """read and parse one line from a genemark gff"""
    elems = l.strip().split("\t")
    gene_dict[
        elems[8].split(" ")[0].strip(",").replace("=", "_")] = (
        elems[0].split(" ")[0], int(
            elems[3]), int(
            elems[4]), elems[6])

def read_prodigal_entry(l, gene_dict):
    """read and parse one line from a prodigal gff"""
    elems = l.strip().split("\t")
    #print elems[8].split(";"), len(elems)
    attrs = dict(
        [(i.split("=")[0], i.split("=")[1])
        for i in elems[8].strip(";").split(";")])
    gene_dict["%s_%s"% (elems[0], attrs["ID"].split("_")[1])] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])

def read_img_entry(l, gene_dict):
    """read and parse one line from an IMG gff"""
    elems = l.strip().split("\t")
    if elems[2] == "CDS":
        try:
            attrs = dict(
                [(i.split("=")[0], i.split("=")[1]) for i in elems[8].strip(";").split(";")])
        except IndexError:
            sys.stderr.write("something went wrong in line; skipping \n%s\n" % l)
            return 
        gene_dict[attrs["ID"]] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])


def read_ncbi_entry(l, gene_dict):
    """read and parse one line from a ncbi gff"""
    elems = l.strip().split("\t")
    if elems[2] == "CDS":
        attrs = dict(
            [(i.split("=")[0], i.split("=")[1])
             for i in elems[8].split(";")])
        gene_dict[
            attrs["Name"]] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])

def read_refseq_entry(l, gene_dict):
    """read and parse one line from a refseq gff"""
    elems = l.strip().split("\t")
    if elems[2] == "CDS":
        attrs = dict(
            [(i.split("=")[0], i.split("=")[1])
             for i in elems[8].split(";")])
        gene_dict[
            attrs["ID"]] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])

def read_genbank_entry(l, gene_dict):
    """read and parse one line from a genbank gff"""
    elems = l.strip().split("\t")
    if elems[2] == "CDS":
        attrs = dict(
            [(i.split("=")[0], i.split("=")[1])
             for i in elems[8].split(";")])
        gene_dict[
            attrs["locus_tag"]] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])



def get_coords(gene_start, gene_end, ali_from, ali_to, strand):
    n_start = ali_from * 3 - 2
    n_end = ali_to * 3 - 2
    # if strand == "+":
    return (gene_start + n_start, gene_start + n_end)
    # else:
    #    return (gene_end - n_end, gene_end  - n_start)


def write_hmm_gff(hmmer_file, out_gff_dir, gene_dict, sample, skip_genes, mode,  pt_models, predicted_phenotypes, description = False):
    #read in pfam accession to description mapping
    if description:
        acc2desc = pt_models.get_pf2desc()
        id2desc = pt_models.get_pt2acc()
        #change dtype of index to string
        id2desc.index = id2desc.index.values.astype('string')
    out_table = []
    with open(hmmer_file, 'r') as hmmer_fo:
        # skip param line and header
        hmmer_fo.readline()
        hmmer_fo.readline()
        hmmer_fo.readline()
        rel_feats_dict = None
        # read in relevant which shall be highlighted
        pts = predicted_phenotypes.split(",") 
        if pts is not None:
            pf2pt = read_rel_feats(pt_models, pts)
        if not gene_dict is None:
            # open global output table
            # open out gffs for reading
            out_gffs = dict([(i, open("%s/%s_%s_important_features.gff"%(out_gff_dir, sample, id2desc.loc[i,].iloc[0].replace(" ", "_").replace("(", "_").replace("(", "_")), 'w')) for i in pts])
            out_gffs["Pfams"] = open("%s/%s_complete_annotation.gff" %(out_gff_dir, sample), 'w')
            for i in out_gffs:
                out_gffs[i].write("""##gff-version 3\n""")
        for l in hmmer_fo:
            elems = l.strip().split("\t")
            # change pfam accession PF00001.3 into PF00001
            gid, acc, name, ali_from, ali_to, ieval = [
                get_protein_acc(elems[0]) if mode == "ncbi" or mode == "refseq" else elems[0], 
                    elems[4].split(".")[0] if pt_models.get_hmm_name() == "pfam" else elems[3].split(".")[0],
                    elems[3], int(elems[17]), int(elems[18]), elems[11]]
            if gene_dict is not None and gid not in gene_dict: 
                if not skip_genes:
                    print >> sys.stderr, gid, "not found in input orf gff file"
                    if mode == "ncbi":
                        print >> sys.stderr, "entry might be outdated, skipping"
                        continue
                    if mode == "img":
                        print >> sys.stderr, "might be a semicolon ; delimiter within an attribute; sample is %s\n" % sample
                        continue
                else:
                    continue
            colour = ""
            # check if this feature is among the highest ranking features and add
            # colour and rank if so 
            #print acc, rel_f
            description_string = ""
            if gene_dict is not None and description and acc in acc2desc.index:
                contig_name, gene_start, gene_end, strand = gene_dict[gid]
                hmm_coords = get_coords(gene_start, gene_end, ali_from, ali_to, strand)
                description_string = ";Description=%s" % (acc2desc.loc[acc,].iloc[0])
                gff_line = "%s\tHMMER\t%s\t%s\t%s\t%s\t%s\t.\tParent=%s;Name=%s;Note=%s%s%s\n" %(contig_name, pt_models.get_name(), hmm_coords[0], hmm_coords[1], ieval, strand, gid, acc, name, description_string, colour)
                out_gffs["Pfams"].write(gff_line)
            if pf2pt is None or acc in pf2pt:
                #if not rel_feats_dict is None:
                #    colour = ";colour=3;rank=%s" % (rel_feats_dict[acc])
                for i in pf2pt[acc]:
                    if not gene_dict is None:
                        out_gffs[i[0]].write(gff_line)
                    out_table.append("%s\t%s\t%s\t%s\t%s\n" %(id2desc.loc[i[0],"accession"], gid, acc, acc2desc.loc[acc, "description"], i[1]))
    if not len(out_table) == 0:
        out_table_df = ps.read_csv(StringIO.StringIO("".join(out_table)), sep = "\t", header = None)
        out_table_df.columns = ["Phenotype", "gene_id", "acc", "description", "cor"]
        out_table_df.sort_values(by = ["Phenotype", "cor"], ascending = [True, False]).to_csv("%s/%s_important_features.dat" % (out_gff_dir, sample), sep = "\t", index = None)
    #close out gffs
    if not gene_dict is None:
        for i in out_gffs:
            out_gffs[i].close()

def run(in_file, out_gff_f, gene_gff_f, sample, gene_gff_mode, pt_models, predicted_pts = None):
    """run feature mapping"""
    gene_dict = None
    if not gene_gff_f is None:
        gene_dict = read_gff(gene_gff_f, gene_gff_mode)
    write_hmm_gff(in_file, out_gff_f, gene_dict, sample, False, gene_gff_mode, pt_models, predicted_pts, True)
