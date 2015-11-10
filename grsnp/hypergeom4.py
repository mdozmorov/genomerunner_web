#!/usr/bin/env python2
from __future__ import division
import argparse
import collections
import math
import sys
import logging
from logging import FileHandler,StreamHandler
#from bx.intervals.intersection import IntervalTree
from scipy.stats import hypergeom
import numpy as np
import scipy
import pdb
import os
import json
import rpy2.robjects as robjects
import gzip
import tarfile
import traceback
import StringIO
import dbcreator_ucsc as bedfilecreator
import textwrap
import subprocess
import sys
import commands
import mako
import simplejson
import zipfile
import inspect
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import grsnp.dbcreator_util as grsnp_util
import random
import string
import collections

# Logging configuration
logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
logger.propagate = 0

# This line outputs logging info to the console
matrix_outpath = None
detailed_outpath = None
progress_outpath = None
run_files_dir = None
console_output = False 
print_progress = False




def get_overlap_statistics(gf,fois):
    """Returns a dictionary with indicating how many hits exist for each foi against the gf
    gf: filepath for GF
    fois: list of FOI filepaths
    """
    results = []
    out = ""
    # use temporary files instead of piping out to console because large amounts of output to console can cause deadlock
    # this creates unique random file names
    tmp_path = get_tmp_file('grsnptmp')        
    tmp_error_path = get_tmp_file('grsnperrortmp')
    tmp_file = open(tmp_path,'wb')
    tmp_error_file = open(tmp_error_path,'wb')

    try:
        # Runs overlapStatistics with preprocessed background stats if they exist
        out = subprocess.Popen(["overlapStatistics"] + [gf] + fois,stdout=tmp_file,stderr=tmp_error_file)
        out.wait()
        tmp_file.close()
        tmp_error_file.close()
        tmp = open(tmp_path).read()
        tmp_er = open(tmp_error_path).read()
        if tmp_er != "": logger.error(tmp_er)
        if tmp[:6] == "ERROR:": 
            logger.error(tmp[7:])
            raise Exception(tmp)

        for x in tmp.split("\n")[1:]:
            if x != "":
                tmp = x.split("\t")
                foi_name,n,hit_count = os.path.split(tmp[0])[-1],tmp[2],tmp[3]
                results.append({"queryfile": foi_name,"queryregions": int(n),"intersectregions": int(hit_count),"indexregions": int(tmp[1])})
        # remove the temporary output files
        if os.path.exists(tmp_path): os.remove(tmp_path)
        if os.path.exists(tmp_error_path): os.remove(tmp_error_path)
    except Exception, e:       
        if not tmp_file.closed: tmp_file.close()
        if not tmp_error_file.closed: tmp_error_file.close()
        # remove the temporary output files
        if os.path.exists(tmp_path): os.remove(tmp_path)
        if os.path.exists(tmp_error_path): os.remove(tmp_error_path)
        logger.error(traceback.format_exc())
        raise e
    return results

def get_tmp_file(prefix):
    tmp_path = prefix + "_" + ''.join(random.choice(string.lowercase+string.digits) for _ in range(32))+'.tmp'
    while (os.path.exists(tmp_path)):
        tmp_path = prefix + "_" + ''.join(random.choice(string.lowercase+string.digits) for _ in range(32))+'.tmp'
    return tmp_path

def get_bgobs(bg,gf,root_data_dir,organism): 
    ''' Check if pre-calculated GF and background overlap data exist.
    If they do not, it manually calculates them.
    '''
    # get the grsnp_db_[filt] folder
    filt_grsnp_db = gf.replace(root_data_dir,"").lstrip("/").split("/")[0]
    bkg_overlap_path = os.path.join(root_data_dir,filt_grsnp_db,organism,'bkg_overlaps.gr')

    # See if pre-calculated values exist
    if os.path.exists(bkg_overlap_path):       
        data = open(bkg_overlap_path).read().split("\n")
        data = [x.split("\t") for x in data if x != ""]
        d_gf = [x[1] for x in data if os.path.join(root_data_dir,x[0]) == gf and x[1]  != ""]

        if len(d_gf) != 0:
            bg_obs = [x.split(":")[1] for x in d_gf[0].split(",") if x.split(":")[0] == os.path.basename(bg)]
            if len(bg_obs) != 0:
                logger.info("Pre-calculated values found for background and {} ".format(base_name(gf)))
                return bg_obs[0]
    # manually get overlap values
    logger.info("Calculating overlap stats on background and {}".format(base_name(gf)))
    _write_progress("Calculating overlap stats on background and {}".format(base_name(gf)))
    result = get_overlap_statistics(gf,[bg])
    try:
        result = int(result[0]["intersectregions"])
    except Exception, e:
        result = None
        logger.error(traceback.format_exc())
    return result

def output_p_value(foi_obs,n_fois,bg_obs,n_bgs,foi_path,gf_path,background_path,run_randomization_test=False):    
    """Return the shrunken odds-ratio and signed p-value of all FOIs against the GF so they can be written to
    matrix files. Outputs stats to the detailed results file.
    """
    global logger_path
    foi_name = base_name(foi_path)
    gf_name = base_name(gf_path)
    sign,pval,odds_ratio,shrunken_or,ci_lower,ci_upper = calculate_p_value_odds_ratio(foi_obs,n_fois,bg_obs,n_bgs,foi_name,gf_path)

    if sign == 1 or str(odds_ratio) == "inf":
        direction  = "overrepresented" 
    else: direction =  "underrepresented"

    # calculate the p_rand
    prnd = 1 # default prnd for non-significant results
    if pval > 0.05:
        direction = "nonsignificant"
    else:
        if run_randomization_test: 
            _write_progress("Running randomization test on {}".format(foi_name))
            prnd = p_rand(foi_path,n_fois,background_path,bg_obs,n_bgs,gf_path)  
    
    pval_unmod = pval
    pval = np.power(10,-(np.log10(prnd)- np.log10(pval))) # adjust p_value using randomization test
    # write out to the detailed results file
    strpval,strprnd = "","" 
    if run_randomization_test:
        strpval = "%.2e" % pval if type(pval) != type("") else pval 
        strprnd = "%.2e" % prnd if type(prnd) != type("") else prnd 

   
    write_output("\t".join(map(str, [foi_name.rpartition('/')[-1], foi_obs, n_fois, bg_obs, n_bgs, 
                _format_type(odds_ratio), 
                _format_type(ci_lower), 
                _format_type(ci_upper), 
                _format_type(shrunken_or),                 
                "%.2e" % pval_unmod if type(pval_unmod) != type("") else pval_unmod,
                strprnd,strpval])) + "\n",detailed_outpath)

    if pval < 1E-307:
        # set to value obtained from sys.float_info.min_10_exp
        pval = 1E-306   
    return [sign * pval,shrunken_or]

def _format_type(num):
    ''' Sets format to be either scientific or float depending on num value
    '''
    if type(num) != type(""):
        if num > 100 or num < 0.01:
            return "%.2e" % num
        else:
            return "%.2f" % num
    else:
        return num

def p_rand(foi_path,n_fois,background_path,bg_obs,n_bgs,gf_path):
    ''' Calculated by generating 'num' random feature files and running them against gf_path.
    Calculates the mean of the p_values for the overrepresented and underrepresented random features separately.
    '''
    num = 10
    rnds_paths = generate_randomsnps(foi_path,background_path,n_fois,gf_path,num)
    rnd_stats = get_overlap_statistics(gf_path,rnds_paths)
    p_rand = [1]
    for r in rnd_stats:
        sign,pval,odds_ratio,_,_,_ = calculate_p_value_odds_ratio(r["intersectregions"],r["queryregions"],bg_obs,n_bgs,base_name(foi_path),gf_path)
        p_rand.append(pval)
    return np.min(p_rand)

def calculate_p_value_odds_ratio(foi_obs,n_fois,bg_obs,n_bgs,foi_name,gf_path):
    """Calculates the p-value,confidence intervals and the shrunken odds ratio.
    Returns [sign,pval,odds_ratio,shrunken_or,ci_lower,ci_upper]
    """
    bg_obs,n_bgs = int(bg_obs),int(n_bgs)
    ctable = [[foi_obs, n_fois-foi_obs],
              [bg_obs-foi_obs,n_bgs-n_fois-(bg_obs-foi_obs)]]
    # Ensure there are no negative values in the ctable
    do_chi_square = True
    for i in ctable:
        for k in i:
            if k < 0:
                logger.warning("Cannot calculate p-value for {} and {}. Is the background too small? foi_obs {}, n_fois {}, bg_obs {}, n_bgs {}".format(base_name(gf_path),foi_name,foi_obs,n_fois,bg_obs,n_bgs))
                return [1,1,1,1,1,1]
            if k < 5:
                do_chi_square = False
    # check for zeros and add 0.5 if one of the cells is 0
    if ctable[0][0] == 0 or ctable[0][1] == 0 or ctable[1][0] == 0 or ctable[1][1] == 0:
        ctable[0][0] += 0.5
        ctable[0][1] += 0.5
        ctable[1][0] += 0.5
        ctable[1][1] += 0.5
    # if bg_obs < foi_obs:
    #     odds_ratio, pval = "nan", 1
    #     logger.warning("P-value cannot be calculated for {} and {} (pvalue = 1.0, odds_ratio = 'nan'). Number of SNPs overlapping with GF > number of background SNPs overlapping with GF. foi_obs {}, n_fois {}, bg_obs {}, n_bgs {}".format(gf_name,foi_name,foi_obs,n_fois,bg_obs,n_bgs))
    # else:  
    if do_chi_square:
        chi_result = scipy.stats.chi2_contingency(ctable)
        pval = chi_result[1]
        odds_ratio = (ctable[0][0]*ctable[1][1])/(ctable[0][1]*ctable[1][0])
    else:  
        odds_ratio, pval = scipy.stats.fisher_exact(ctable)
    # Adjustments of outliers
    if pval == 1.0:
        odds_ratio = 1
    if odds_ratio == 0.0:
        odds_ratio = sys.float_info.min
    if np.isinf(odds_ratio):
        odds_ratio = sys.float_info.max

    # calculate the shrunken odds ratio
    log_or = scipy.log(odds_ratio)
    conf_coe = 1.96 # the confidence coefficient of a standard norm dist
    # calculate the standard error
    se = math.sqrt(1/ctable[0][0] + 1/ctable[1][0] + 1/ctable[0][1] + 1/ctable[1][1])
    # calculate the upper and lower confidence interval
    ci_upper = scipy.exp(log_or + conf_coe * se)
    ci_lower = scipy.exp(log_or - conf_coe * se)
    # Precaution against CI overflow
    if np.isinf(ci_upper):
        ci_upper = sys.float_info.max
    if ci_lower == 0.0:
        ci_lower = sys.float_info.min
    # shrunken_or is the ci (either upper or lower) that is closest to 1
    if ci_lower<1 and ci_upper>1:
        shrunken_or,odds_ratio = 1,1
    else:
        # find which value is closer to 1
        ci_index = scipy.array([[abs(math.log(ci_lower)),abs(math.log(ci_upper))]]).argmin()
        shrunken_or = [ci_lower,ci_upper][ci_index]

    sign = -1 if (odds_ratio < 1) else 1
    return [sign,pval,odds_ratio,shrunken_or,ci_lower,ci_upper]



def generate_randomsnps(foi_path,background,n_fois,gf_path,num):
    paths = []
    out_dir = os.path.join(os.path.dirname(foi_path),"random")
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    for n in range(num):        
        rnd_snp_path = os.path.join(out_dir,"random{}_".format(n)+base_name(foi_path)+".bed")
        # generate random snps from background
        if background.endswith('.gz'):
            command = "zcat {} | shuf -n {} | cut -f 1-3 > {}".format(background,str(n_fois),rnd_snp_path)
            out = commands.getstatusoutput(command)
        else:        
            command = "shuf -n {} {}  | cut -f 1-3 > {}".format(str(n_fois),str(background),rnd_snp_path)
            out = commands.getstatusoutput(command)
        paths.append(rnd_snp_path)
    return paths

def _chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]

def get_annotation(foi,gfs):
    """
    fois: list of FOI filepath
    gfs: filepaths for GF
    """
    results = []
    out = ""
    # use temporary files instead of piping out to console because large amounts of output to console can cause deadlock
    # this creates unique random file names
    tmp_path = get_tmp_file('grsnptmp')
    tmp_error_path = get_tmp_file('grsnperrortmp')
    
    tmp_file = open(tmp_path,'wb')
    tmp_error_file = open(tmp_error_path,'wb')
    try:
        out = subprocess.Popen(["annotationAnalysis"] + [foi] + gfs,stdout=tmp_file,stderr=tmp_error_file) # TODO enable ["--print-region-name"]
        out.wait()
        tmp_file.close()
        tmp_error_file.close()
        tmp = open(tmp_path).read()
        tmp_er = open(tmp_error_path).read()
        if tmp_er != "": logger.error(tmp_er)
        if tmp[:6] == "ERROR:": 
            logger.error(tmp[7:])
            raise Exception(tmp)
        # remove the temporary output files
        if os.path.exists(tmp_path): os.remove(tmp_path)
        if os.path.exists(tmp_error_path): os.remove(tmp_error_path)

    except Exception, e:
        if not tmp_file.closed: tmp_file.close()
        if not tmp_error_file.closed: tmp_error_file.close()
        # remove the temporary output files
        if os.path.exists(tmp_path): os.remove(tmp_path)
        if os.path.exists(tmp_error_path): os.remove(tmp_error_path)
        logger.error(traceback.format_exc())
        raise e
    return tmp

# Writes the output to the file specified.  Also prints to console if console_output is set to true
def write_output(content,outpath=None):
    if outpath:
        with open(outpath,"a") as out:
            out.write(content)
    if console_output:
        print >> sys.stderr, content

# Collect lines from a file into a list
def read_lines(path):
    elems = []
    with open(path) as h:
        for line in h:
            if line.strip():
                elems.append(line.strip())
    return elems

def base_name(k):
    return os.path.basename(k).split(".")[0]

def _write_progress(line):
    """Saves the current progress to the progress file
    """
    global progress_outpath
    if progress_outpath:
        global curprog, progmax
        progress = {"status": line, "curprog": curprog,"progmax": progmax}
        with open(progress_outpath,"wb") as progfile:
            progfile.write(json.dumps(progress))
    if print_progress:
        print line


def _write_head(content,outpath):
    f = front_appender(outpath)
    f.write(content)
    f.close()


def check_background_foi_overlap(bg,fois):
    """ Calculates the overlap of the FOIs with the background.
    Removes FOIs that are poorly formed with the background.
    """
    _write_progress("Validating FOIs against background")
    good_fois = []
    if len(fois) == 0:
        return [[],[]]
    # Runs overlapStatistics on background and FOIs
    foi_bg_stats =  get_overlap_statistics(bg,fois)
    for f in foi_bg_stats:
        isgood = True
        foi_name,n_bgs,n_fois,foi_in = f["queryfile"],f["indexregions"],f["queryregions"],f["intersectregions"]
        if n_fois < 5:
            isgood = False
            logger.warning("Number of SNPs in {} < 5. Removing it from analysis.".format(foi_name))
        elif n_bgs < n_fois:
            isgood = False
            logger.warning("Number of SNPs in {} > than in background. Removing it from analysis.".format(foi_name))
        if isgood:
            # ensure that overlapStatistics output filename with extension for queryFile field
            good_fois.append([x for x in fois if os.path.basename(x) == f["queryfile"]][0])
        if foi_in < n_fois:
            logger.warning("{} out of {} {} SNPs are not a part of the background. P-value are unreliable. Please, include all SNPs in the background and re-run analysis.".format(n_fois-foi_in,n_fois,foi_name))
    return [foi_bg_stats, good_fois]
                                                                                                                       


def get_description(gf,track_descriptions):
    ''' Get the GF feature description
    '''
    desc = [x[1] for x in track_descriptions if x[0] == gf and len(x[0]) > 1]
    if len(desc) is not 0: 
        return desc[0]
    else: 
        return "No Description"




def _zip_run_files(fois,gfs,bg_path,outdir,id=""):
    '''
    File paths of FOIs and GFs as a list. Gathers all the files together in one zipped file
    '''    

    # zip annotation result folder if it exists
    anno_dir = os.path.join(outdir,'annotations')
    if os.path.exists(anno_dir):
        z_ano = zipfile.ZipFile(os.path.join(outdir,'annotations.zip'),'a')
        for f in os.listdir(anno_dir):
            z_ano.write(os.path.join(anno_dir,f),f)
        z_ano.close()

    # zip processed FOIs directory
    proc_dir = os.path.join(outdir,'processed_fois')
    if os.path.exists(proc_dir):
        z_ano = zipfile.ZipFile(os.path.join(outdir,'processed_fois.zip'),'a')
        for f in os.listdir(proc_dir):
            z_ano.write(os.path.join(proc_dir,f),f)
        z_ano.close()

    f = open(os.path.join(outdir,"gr_log.txt"))
    f_log = f.read()
    f.close()
    path_settings,f_sett =os.path.join(outdir,".settings"),""
    if os.path.exists(path_settings):
        f = open(path_settings)
        f_sett = f.read() + "\n###LOG###\n"
        f.close()
    new_log_path = os.path.join(outdir,"gr_log.txt")
    new_log = open(new_log_path,'wb')
    new_log.write(f_sett+f_log)
    new_log.close()
    tar_path = os.path.join(outdir,'GR_{}.tar'.format(id))
    tar = tarfile.TarFile(tar_path,"a")    
    output_files =  [os.path.join(outdir,x) for x in os.listdir(outdir) if x.endswith(".txt") or x.endswith(".pdf") or x.endswith('.zip')]
    fls = output_files
    for f in fls:
        tar.add(f,os.path.basename(f))
    tar.close()
    tar_file = open(tar_path,'rb')
    with gzip.open(tar_path+".gz","wb") as gz:
        gz.writelines(tar_file)
    tar_file.close()
    if os.path.exists(tar_path): os.remove(tar_path)

def validate_filenames(file_paths):
    ''' Checks if there are spaces before or after the file extension.
    EX. 'dir1/dir2/test .bed' is not valid. 'dir1/dir2/test.bed is valid.
    'dir1/dir2/test.bed .gz' is not valid.
    '''
    invalid = []
    for file in file_paths:
        for t in os.path.basename(file).split("."):
            if len(t.strip()) != len(t):
                # there are spaces before or after the '.'. Add the file to the list of invalids.
                logger.error("Cannot have space before/in file extension: {}".format(file))
                invalid.append(os.path.basename(file))
    files_wo_ext = [base_name(x) for x in file_paths]
    file_duplicates = [x for x, y in collections.Counter(files_wo_ext).items() if y > 1]
    for f in file_duplicates:
        logger.error("{} exists multiple times. (i.e. 'foi.txt.gz' has same basename as 'foi.txt')".format(f))
        invalid.append(f)
    return invalid


def get_score_strand_settings(gf_path):
    ''' Parses the gf_path and determines if gf is filtered by score and/or strand.
    '''
    str_strand,str_scorethresh = "Strand: Both","Score threshold: NA"
    gfsplit = gf_path.split("/grsnp_db_")
    if len(gfsplit) == 2:
        str_score_strand = gfsplit[-1].split("/")[0].split("_")
        for s in str_score_strand:
            if s.isdigit():
                str_scorethresh = "Score threshold: " + s
            else:
                str_strand = "Strand: " + s
    return str_strand + "\t" + str_scorethresh

def validate_rsids(foi_path):
    ''' Checks if the first line is contains an rsIDs. If it does, the file is sorted
    '''   
    if foi_path.endswith('.gz'):
        infile = gzip.open(foi_path)
    else:
        infile = open(foi_path)
    line = infile.readline().rstrip('\n')
    infile.close()
    if line[:2] == 'rs':
        # sort the rsIDs in place
        script =  "sort -k1,1 " + '-o ' +foi_path+' '+ foi_path
        out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out.wait()
        tmp = out.stdout.read()
        tmp_er = out.stderr.read()
        if tmp_er != "":
            logger.error(tmp_er)
            raise Exception(tmp_er) 
        return True
    else:
        return False

def preprocess_fois(fois,run_files_dir,root_data_dir,organism):
    processed_fois = []
    output_dir = os.path.join(run_files_dir,'processed_fois')
    # Sort the fois 
    out = ""
    data_dirs = [os.path.join(root_data_dir,'grsnp_db',organism),
                os.path.join(root_data_dir,'custom_data')]
    try: 
        for f in fois:
            if len([x for x in data_dirs if x in f]) == 0:
                # extract file if it is gzipped.
                unzipped_f = f
                if f.endswith('.gz'):
                    out = subprocess.Popen(["gunzip {}".format(f)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    out.wait()
                    # filepath without the .gz extension
                    unzipped_f = ".".join(unzipped_f.split(".")[:-1])
                # copy the FOI to the output_dir
                out_fname = os.path.split(unzipped_f)[1]
                if not out_fname.endswith('.bed'):
                    out_fname = base_name(out_fname)+".bed"
                out_f = os.path.join(output_dir,out_fname) # the processed foi file
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                out = subprocess.Popen(['cp {} {}'.format(unzipped_f,out_f)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                out.wait()
                # remove the header from the files
                grsnp_util.remove_headers(out_f)

                # Check if items are rsIDs. Convert to bed coordinates if they are
                if validate_rsids(out_f):
                    # check if a file exists in the database for rsID conversion and construct the path to it
                    rsid_path = os.path.join(root_data_dir,'custom_data','rsid_conversion',organism)
                    if not os.path.exists(rsid_path):
                        logger.error('rsID conversion not available for this organism. Feature set {} removed'.format(f)) 
                        continue
                    files = [x for x in os.listdir(rsid_path) if os.path.isfile(os.path.join(rsid_path,x)) and x.endswith('.bed')]
                    # if conversion files found, perform conversion
                    if len(files) > 0:
                        rsid_path = os.path.join(rsid_path,files[0])
                        script = """join {} {} -1 1 -2 4 -o 2.1 -o 2.2 -o 2.3 -o 2.4 -o 2.5 -o 2.6 | sed 's/\ /\t/g' > {}.temp""".format(out_f,rsid_path,out_f)
                        out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        out.wait()
                        tmp_er = out.stderr.read()
                        if tmp_er != "":
                            logger.error(tmp_er)
                            raise Exception(tmp_er)
                        # we remove the original out_f FOI file and replace with the out_f.temp created with the join command above
                        os.remove(out_f)
                        out = subprocess.Popen(['cp {} {}'.format(out_f+'.temp',out_f)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        out.wait()  
                        os.remove(out_f+'.temp')
                    else:
                        logger.error('rsID conversion not available for this organism. Analysis terminated.') 
                        return [] 

                # sort the FOI and bgzip
                grsnp_util.sort_convert_to_bgzip(out_f,out_f+".gz")
                script = "tabix -f " + out_f+".gz"
                out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE)
                out.wait()
                # check if processed FOI file is empty
                if os.stat(out_f+".gz")[6]!=0:                   
                    processed_fois.append(out_f + ".gz")
                else:
                    logger.error("{} is empty. Removing.")                    
            else:
                processed_fois.append(f)
                # print "{} is part of the database".format(f)
    except Exception, e:
        logger.error("Error while processing the FOIs")
        logger.error(traceback.format_exc())
        raise Exception(e)
    return processed_fois

def preprocess_gf_files(file_paths,root_data_dir,organism):
    ''' Used to preprocess the GF files and the background file.
    Returns gzipped file paths
    '''
    processed_files = []
    data_dirs = [os.path.join(root_data_dir,'grsnp_db',organism),
                os.path.join(root_data_dir,'custom_data')]
    try:
        for out_f in file_paths:
            # if the file is not uploaded by the user, do not preprocess
            if len([x for x in data_dirs if x in out_f]) == 0:
                # unzip the file
                unzipped_f = out_f
                if out_f.endswith('.gz'):
                    out = subprocess.Popen(["gunzip {}".format(out_f)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    out.wait()
                    # filepath without the .gz extension
                    unzipped_f = ".".join(unzipped_f.split(".")[:-1])

                out_fname = os.path.split(unzipped_f)[1]
                if not out_fname.endswith('.bed'):
                    out_fname = base_name(out_fname)+".bed"
                # replace the ".txt" extension with ".bed"
                output_dir = os.path.split(unzipped_f)[0]
                out_bed_f = os.path.join(output_dir,out_fname) # the processed file

                grsnp_util.remove_headers(unzipped_f)
                # perform rsid conversion
                if validate_rsids(unzipped_f):
                    # check if a file exists in the database for rsID conversion and construct the path to it
                    rsid_path = os.path.join(root_data_dir,'custom_data','rsid_conversion',organism)
                    if not os.path.exists(rsid_path):
                        logger.error('rsID conversion not available for this organism. Feature set {} removed'.format(unzipped_f)) 
                        continue
                    files = [x for x in os.listdir(rsid_path) if os.path.isfile(os.path.join(rsid_path,x)) and x.endswith('.bed')]
                    # if conversion files found, perform conversion
                    if len(files) > 0:
                        rsid_path = os.path.join(rsid_path,files[0])
                        # sort the RSID
                        script = """sort -k1,1 -o {} {}""".format(unzipped_f,unzipped_f)
                        out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        out.wait()   
                        # join the RSID with the SNP data in custom_data
                        script = """join {} {} -1 1 -2 4 -o 2.1 -o 2.2 -o 2.3 -o 2.4 -o 2.5 -o 2.6 | sed 's/\ /\t/g' > {}.temp""".format(unzipped_f,rsid_path,unzipped_f)
                        out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        out.wait()             
                        tmp_er = out.stderr.read()
                        if tmp_er != "":
                            logger.error(tmp_er)
                            raise Exception(tmp_er)
                        # we remove the original unzipped_f FOI file and replace with the unzipped_f.temp created with the join command above
                        out = subprocess.Popen(['cp {} {}'.format(unzipped_f+'.temp',out_bed_f)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        out.wait()  
                        os.remove(unzipped_f+'.temp')
                    else:
                        logger.error('rsID conversion not available for this organism. Analysis terminated.') 
                        return []

                # sort the FOI
                grsnp_util.sort_convert_to_bgzip(out_bed_f,out_bed_f+".gz")
                # check if processed FOI file is empty
                if os.stat(out_bed_f+".gz")[6]!=0:                   
                    # add the processed file (the ".bed" file) to the list.
                    processed_files.append(out_bed_f + ".gz")
                else:
                    logger.error("{} is empty. Removing.")
                    
            else:
                processed_files.append(out_f)
                # print "{} is part of the database".format(out_f)
    except Exception, e:
       logger.error("Error while processing the GF/background files")
       logger.error(traceback.format_exc())
       raise Exception(e)
    return processed_files




def run_hypergeom(fois, gfs, bg_path,outdir,job_name="",zip_run_files=False,bkg_overlaps_path="",root_data_dir = "" ,run_annotation=True,run_randomization_test=False,pct_score="",organism = ""):
    global formatter
    global detailed_outpath,matrix_outpath, progress_outpath, curprog, progmax,run_files_dir
    if not os.path.exists(os.path.normpath(outdir)): os.mkdir(os.path.normpath(outdir))
    fh = logging.FileHandler(os.path.join(outdir,'gr_log.txt'))
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    run_files_dir = outdir
    curprog,progmax = 0,1
    logger.propagate = False

    try:
        track_descriptions = []
        logger.info("Enrichment analysis started")
        decriptions_path = os.path.join(root_data_dir,"grsnp_db",organism,"gf_descriptions.txt")
        if os.path.exists(decriptions_path):
            track_descriptions = [x.split("\t") for x in open(decriptions_path).read().split("\n") if x != ""]
        # set output settings
        detailed_outpath =  os.path.join(outdir, "detailed.txt") 
        matrix_outpath = os.path.join(outdir,"matrix_PVAL.txt")
        matrix_sor_outpath = os.path.join(outdir,"matrix_OR.txt")
        progress_outpath = os.path.join(outdir,".prog")
        _write_progress("Starting analysis.")
        f = open(matrix_outpath,'wb') 
        f.close()
        f = open(matrix_sor_outpath,'wb')
        f.close()
        f = open(detailed_outpath,'wb')
        f.close()

        # Read in the paths
        fois = [line for line in read_lines(fois) if not line.endswith(".tbi")]
        gfs = [line for line in read_lines(gfs) if not line.endswith(".tbi")]
        # check if there are spaces in invalid parts of the file name
        invalid_names = validate_filenames(fois + gfs + [bg_path])
        if len(invalid_names) != 0:
            logger.error("The following file(s) have invalid file names:\n" + "\n".join(invalid_names))
            _write_progress("ERROR: Files have invalid filenames. See log file. Terminating run. See Analysis Log.")
            print "ERROR: Files have invalid filenames. See log file. Terminating run. See Analysis Log."
            return         
        if bg_path.endswith(".tbi"):
            logger.error("Background has invalid extension (.tbi). Terminating run.")
            _write_progress("ERROR: Background has invalid extension (.tbi). Terminating run. See Analysis Log.")
            print "ERROR: Background has invalid extension (.tbi). Terminating run. See Analysis Log."
            return

        # pre-process the FOIs
        fois = preprocess_fois(fois,run_files_dir,root_data_dir,organism)
        if len(fois) == 0:
            logger.error('No valid FOIs to supplied')
            _write_progress("ERROR: No valid FOI files supplied. Terminating run. See Analysis Log.")
            return
        # pre-process the GFs and the background
        bg_path = preprocess_gf_files([bg_path],root_data_dir,organism)[0]
        gfs = preprocess_gf_files(gfs,root_data_dir,organism)
        # Validate FOIs against background. Also get the size of the background (n_bgs)
        foi_bg,good_fois  = check_background_foi_overlap(bg_path,fois)   
        write_output("\t".join(map(base_name,good_fois))+"\n", matrix_outpath)
        write_output("\t".join(map(base_name,good_fois))+"\n", matrix_sor_outpath)
        write_output("\t".join(['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio',"ci_lower","ci_upper", 'shrunken_odds_ratio', 'p_val','p_rand' if run_randomization_test else "",'p_mod' if run_randomization_test else ""]) + "\n",detailed_outpath)
        curprog,progmax = 0,len(gfs)
        # check if any good fois exist after background filtering
        if len(good_fois) == 0:
            logger.error('No valid FOIs to supplied')
            _write_progress("ERROR: No valid FOI files supplied. Terminating run. See Analysis Log.")
            return
        # remove old detailed enrichment result files if they exit
        enr_path =  os.path.join(run_files_dir,"enrichment")
        for f in good_fois:
            f_path = os.path.join(enr_path, base_name(f)+'.txt')
            if os.path.exists(f_path): os.remove(f_path)
        _write_progress("Performing calculations on the background.")
        for gf in gfs: 
            current_gf = base_name(gf)    
            _write_progress("Performing Hypergeometric analysis for {}".format(base_name(gf)))
            write_output("###"+base_name(gf)+"\t"+get_score_strand_settings(gf)+"\t"+get_description(base_name(gf),track_descriptions)+"###"+"\n",detailed_outpath)
            res = get_overlap_statistics(gf,good_fois) 

            # calculate bg_obs
            bg_obs = get_bgobs(bg_path,gf,root_data_dir,organism)
            if bg_obs == None: 
                logger.error("Skipping {}".format(gf))
                continue

            n_bgs = foi_bg[0]["indexregions"]  

            # calculate the pvalues and output the matrix line for the current gf
            pvals,sors = [],[] # sors = shrunken odds-ratios
            for i in range(len(good_fois)):
                [pvalue,shrunken_or] = output_p_value(res[i]["intersectregions"],res[i]["queryregions"],bg_obs,n_bgs ,good_fois[i],gf,bg_path,run_randomization_test)
                pvals.append(str(pvalue))
                sors.append(str(shrunken_or))

            # output the matrices file lines
            write_output("\t".join([base_name(gf)] + pvals)+"\n",matrix_outpath)
            write_output("\t".join([base_name(gf)] + sors)+"\n",matrix_sor_outpath)

            curprog += 1
        if run_annotation:
            logger.info("Annotation started")
            annot_outdir = os.path.join(outdir,"annotations")
            if not os.path.exists(annot_outdir): os.mkdir(annot_outdir)
            curprog,progmax = 0,len(fois)
            for f in fois:                
                _write_progress("Running Annotation Analysis for {}.".format(base_name(f)))
                logger.info("Running annotation analysis for {}".format(base_name(f)))
                for i, g in enumerate(_chunks(gfs,100)):
                    with open(os.path.join(annot_outdir,base_name(f) + str(i) + ".txt"),"wb") as wr:
                        anot = get_annotation(f,g).split("\n")
                        anot[0] = anot[0].replace("Region\t\t","Region\t")
                        wr.write("Region"+"\t"+"\t".join(base_name(x) for x in reversed(anot[0].split("\t")[1:])) + "\tTotal") # annotationAnalysis column order is reverse of input order
                        for ind, a in enumerate(anot[1:]):
                            if a.strip() != "":
                                cur_row = a.split("\t")
                                wr.write("\n" + str(ind) + "|"+"\t".join(cur_row + [str(sum([int(x) for x in cur_row[1:] if x != ""]))]))                            
                curprog += 1
            logger.info("Annotation finished")
        if zip_run_files:
            _write_progress("Preparing run files for download")
            _zip_run_files(fois,gfs,bg_path,outdir,job_name)
        curprog,progmax = 1,1
        _write_progress("Analysis Completed")
        logger.info("Analysis Completed")       
    except Exception, e: 
        logger.error( traceback.print_exc())
        _write_progress("Run crashed. See end of log for details.")
        raise Exception(e)

class front_appender:
    '''
    Appends content to start of file.
    '''
    def __init__(self, fname, mode='a'):
        self.__write_queue = []
        self.__old_content = ""
        if mode == 'a': 
            self.__old_content = open(fname).read()
        self.__f = open(fname, 'w')


    def write(self, s):
        self.__write_queue.append(s)

    def close(self):
        self.__f.writelines(self.__write_queue + [self.__old_content])
        self.__f.close()


def _load_minmax(path):
    data = {}
    if not os.path.exists(path):
        return data
    score = [x for x in open(path).read().split("\n") if x != ""]
    for s in score:
        name,min_max = s.split('\t')
        data[name] = min_max
    return data

def main():
        global matrix_outpath, detailed_outpath, progress_outpath, run_files_dir, console_output, print_progress, print_progress
        print_progress = True     
        parser = argparse.ArgumentParser(description="Enrichment analysis of several sets of SNPs (FOIs) files against several genomic features (GFs). Example: python hypergeom4.py foi_full_names.txt gf_full_names.txt /path_to_background/snp137.bed.gz")
        parser.add_argument("fois", nargs=1, help="Text file with paths to FOI files (unless -p used). Required") 
        parser.add_argument("gfs" ,nargs=1, help="Text file with pathrs to GF files (unless -p used). GF files may be gzipped. Required")
        parser.add_argument("bg_path", nargs=1, help="Path to background, or population of all SNPs. Required")
        parser.add_argument("--run_annotation" , "-a", help="Run annotation analysis", action="store_true" )
        parser.add_argument("--run_files_dir" , "-r", nargs="?", help="Set the directory where the results should be saved. Use absolute path. Example: /home/username/run_files/.", default=os.getcwd())
        parser.add_argument("--pass_paths", "-p", help="Pass fois and gfs as comma separated paths. Paths are saved in .fois and .gfs file.", action="store_true")
        parser.add_argument("--data_dir" , "-d", nargs="?",type=str, help="Set the directory containing the database. Required for rsID conversion. Use absolute path. Example: /home/username/db_#.##_#.##.####/.", default="")
        parser.add_argument('--organism','-g', nargs="?", help="The UCSC code of the organism to use. Required for rsID conversion. Default: hg19 (human).", default="hg19")
        args = vars(parser.parse_args())
        if args['organism'] is None:
            print "--organism cannot be blank"
            return None
        if args['run_files_dir'] is None:
            print "--run_files_dir cannot be blank"
            return None
        if args["pass_paths"]: 
            gf = args["gfs"][0].split(",")      
            foi = args["fois"][0].split(",")  
            if not os.path.exists(args["run_files_dir"]) and args['run_files_dir'] != "":
                os.mkdir(args['run_files_dir'])
            # write out the passed gf and foi paths into .gfs and .fois files.
            args["gfs"][0],args["fois"][0] = os.path.join(args["run_files_dir"],".gfs"),os.path.join(args["run_files_dir"],".fois")       
            with open(args["gfs"][0],'wb') as writer:       
                writer.write("\n".join(gf))     
            with open(args["fois"][0],"wb") as writer:      
                writer.write("\n".join(foi))
        run_hypergeom(args["fois"][0],args["gfs"][0],args["bg_path"][0],args["run_files_dir"],"",False,"",args['data_dir'],args["run_annotation"],run_randomization_test=False,organism=args['organism'])




if __name__ == "__main__":
    main()
