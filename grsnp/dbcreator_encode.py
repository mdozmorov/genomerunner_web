#!/usr/bin/env python2
import sys
import errno
import logging
from logging import FileHandler,StreamHandler
import os
import ftplib
import sqlite3
import string
from contextlib import closing
import subprocess
import argparse
import gzip
import re
import collections
import copy
import traceback  as trace
import pdb
from xml.sax.saxutils import quoteattr as xml_quoteattr
from bs4 import BeautifulSoup
import urllib2
from grsnp.dbcreator_util import *
from time import sleep
from collections import namedtuple
import numpy as np

logger = logging.getLogger('genomerunner.dbcreator')

# connection information for the ucsc ftp server

username = 'anonymous'
password = ''

download_dir = ""
DEBUG = True

class GF_ALREADY_EXISTS(Exception):
    pass

gf_description_colnames = ["file_name","full_path", "URL", "full_description", "category", "category_desc", "cell", "cell_desc",
	"factor", "factor_desc", "source", "source_desc"]
GF_description = namedtuple("gf_description",gf_description_colnames)

gf_descriptions = {}
organism = ""


# downloads the specified file from ucsc.  Saves it with a .temp extension untill the download is complete.
def download_encode_file(organism,gf_group,gf_file):
	''' Downloads the gf_file from the UCSC ftp server and saves it
	in a folder with the same name as the organism.
	'''
	global download_dir
	outputpath = ''
	try:
		if os.path.exists(download_dir) == False and download_dir != '':
			logger.info( "creating directory {}".format(download_dir))
			os.makedirs(download_dir)
	except Exception, e:
		logger.warning( e)
		logger.warning("Could not create folder at {} for {}".format(download_dir,gf_file))
		return '' 
	
	try:
		outputpath = os.path.join(download_dir,gf_file)
		if not os.path.exists(outputpath):
			ftp = ftplib.FTP(gf_grp_sett[gf_group]['ftp_server'], timeout=1800) # Connection timeout 0.5h
			ftp.login(username,password)			
			with open(outputpath + ".temp",'wb') as fhandle:  
				global ftp
				logger.info( 'Downloading {} from {}'.format(gf_file,
				gf_grp_sett[gf_group]['ftp_server']+gf_grp_sett[gf_group]['directory'].format(organism)))			
				ftp.cwd(gf_grp_sett[gf_group]['directory'].format(organism))
				ftp.retrbinary('RETR ' + "{}".format(gf_file),fhandle.write)
				os.rename(outputpath+".temp",outputpath)
			ftp.quit()
			# remove header lines
			remove_headers(outputpath)
		else:
			logger.info( '{} already exists, skipping download'.format(outputpath))
	except IOError as e:
		if e.errno == errno.EPIPE: # Broken pipe error handling
			logger.error('FTP download error. Restart the dbcreator. Exiting now...')
			ftp.quit()
			sys.exit(2)
	except Exception, e:
		logger.warning(e)
		logger.warning("Could not download the {} file. Names ARE case sensitive.".format(gf_file))
		return '' 
	return outputpath

def download_roadmap_file(gf_group,gf_file):
	''' Downloads the gf_file from the http server and saves it
	in a folder with the same name as the organism.
	'''
	global download_dir
	outputpath = ''
	
	try:
		if os.path.exists(download_dir) == False and download_dir != '':
			logger.info( "creating directory {}".format(download_dir))
			os.makedirs(download_dir)
	except Exception, e:
		logger.warning( e)
		logger.warning("Could not create folder at {} for {}".format(download_dir,gf_file))
		return ''	
	try:
		outputpath = os.path.join(download_dir,gf_file)
		outputpath = outputpath.replace(".Z.","Z.") # special case: H2A.Z needs to be H2AZ
		if not os.path.exists(outputpath):
			url = "".join([gf_grp_sett[gf_group]['html_server'],gf_grp_sett[gf_group]['directory'],"/",gf_file])
			req = urllib2.urlopen(url)
			CHUNK = 16 * 1024
			logger.info( 'Downloading {} from {}'.format(gf_file,url))
			with open(outputpath + ".temp",'wb') as fp:
				while True:
					chunk = req.read(CHUNK)
					if not chunk: break
					fp.write(chunk)
			os.rename(outputpath+".temp",outputpath)
			# remove header lines
			remove_headers(outputpath)
		else:
			logger.info('{} already exists, skipping download'.format(outputpath))	
	except Exception, e:
		logger.warning(e)
		logger.warning("Could not download the {} file. Names ARE case sensitive.".format(gf_file))
		return '' 
	return outputpath


def get_gf_filepaths(organism,gf_group):
	'''Returns the list of file in the gf_group folder on the ucsc server.
	'''
	ftp = ftplib.FTP(gf_grp_sett[gf_group]['ftp_server'], timeout=1800) # Connection timeout 0.5h
	ftp.login(username,password)

	root_dir = gf_grp_sett[gf_group]['directory'].format(organism)
	ftp.cwd(root_dir)
	gf_paths = []
	# get all folders in the root directory
	ftp.dir(gf_paths.append)
	gf_paths = [x.split()[-1] for x in gf_paths if x[0] == '-']
	ftp.quit()
	# for chromstats in Roadmap we only want certain files
	if gf_group.startswith('chromStates_'):
		gf_paths = [x for x in gf_paths if x.endswith("_dense.bed.gz")]
	# Do some filtering as to which files are included 
	if gf_group in ["wgEncodeCaltechHist","wgEncodeCaltechTfbs","wgEncodeSydhTfbs","wgEncodePsuDnase"]:
		gf_paths = [x for x in gf_paths if x.endswith(".narrowPeak.gz")]
	if gf_group in ["wgEncodeLicrHistone", "wgEncodePsuHistone", "wgEncodeLicrTfbs", "wgEncodePsuTfbs", "wgEncodeUwDnase", "wgEncodeUwDgf"]:
		gf_paths = [x for x in gf_paths  if x.endswith(".broadPeak.gz")]
	if gf_group == "wgEncodePsuHistone":
		gf_paths = [x for x in gf_paths if "InputPk" not in x]
	return gf_paths
	

def get_road_gf_filepaths(gf_group):
	'''Returns the list of the files on the roadmap server for the 'gf_group'.
	Does this by parsing the returned html file for hrefs'''
	url = "".join([gf_grp_sett[gf_group]['html_server'],gf_grp_sett[gf_group]['directory']])
	url_file = gf_group +'.html'
	if os.path.exists(url_file) and DEBUG:
		remotefile = open(url_file).read()
	else:
		remotefile = urllib2.urlopen(url).read().decode('utf-8')
		if DEBUG:
			with open(url_file,'wb') as writer:
				writer.write(remotefile)
	soup = BeautifulSoup(remotefile)
	links = soup.body.find_all('a', href=True)
	# exclude header line links (they lack a '.')
	file_names = [x['href'] for x in links if '.' in x.contents[0]]
	if gf_group in ["chromStates15", "chromStates18", "chromStates25"]:
		file_names = [x for x in file_names if x.endswith('_dense.bed.gz')]
	if gf_group.startswith('Histone_'):
		# filter out the DNase files as these will be processes separately
		file_names = [x for x in file_names if '-DNase' not in x]
		# filter out all non compressed .bed file names
		file_names = [x for x in file_names if not x.endswith('.bed')]
	if gf_group.startswith('DNase_'):
		file_names = [x for x in file_names if '-DNase' in x and "DNase.imputed.narrowPeak.bed.gz" not in x]
	return file_names



def _format_bed(line):
	line = line.split("\t")
	# removes special characters and the '.' used when a field is blank
	if len(line) >= 6:
		line[3] = ''.join(e for e in line[3] if e.isalnum())
		line[4] = line[4] if line[4] != "." else "0"
		line[5] = line[5] if line[5] in ["+","-"] else ""
		line = line[0:6]
	elif len(line) == 5:
		line[3] = ''.join(e for e in line[3] if e.isalnum())
		line[4] = line[4] if line[4] != "." else "0"
		line = line[0:5]
	elif len(line) == 4:
		line[3] = ''.join(e for e in line[3] if e.isalnum())
		line = line[0:4]
	return line

def _format_peak(line):
	'''Handles broadPeak and narrowPeak. Returns line of formated bed as a list of fields
	'''
	line = line.split("\t")
	# removes special characters and the '.' used when a field is blank
	line[3] = ''.join(e for e in line[3] if e.isalnum())
	line[4] = line[6] if line[6] != "." else "0" # for peak data we use the SignalValue column
	line[5] = line[5] if line[5] in ["+","-"] else ""
	line = line[0:6]
	return line

def _format_gapped_peak(line):
	'''Handles gappedPeak. Returns line of formated bed as a list of fields
	'''
	line = line.split("\t")
	# removes special characters and the '.' used when a field is blank
	line[3] = ''.join(e for e in line[3] if e.isalnum())
	line[4] = line[12] if line[12] != "." else "0" # for peak data we use the SignalValue column
	line[5] = line[5] if line[5] in ["+","-"] else ""
	line = line[0:6]
	return line

def _format_bedRNA(line):
	offset = 1
	line = line.split("\t")
	# removes special characters and the '.' used when a field is blank
	if len(line) >= 7:
		line[4] = ''.join(e for e in line[4] if e.isalnum())
		line[5] = line[5] if line[5] != "." else "0"
		line[6] = line[6] if line[6] in ["+","-"] else ""
		line = line[1:7]
	elif len(line) == 6:
		line[4] = ''.join(e for e in line[4] if e.isalnum())
		line[5] = line[5] if line[5] != "." else "0"
		line = line[1:6]
	elif len(line) == 5:
		line[4] = ''.join(e for e in line[4] if e.isalnum())
		line = line[1:5]
	return line




def _get_tier(line,outputdir): # Generating paths for the ENCODE data tables using tiers, and cell types
	CELLS1 = re.compile('Gm12878|K562|H1hesc')
	CELLS2 = re.compile('A549|Cd20ro01778|Cd20ro01794|Cd20|H1neurons|Helas3|Hepg2|Huvec|Imr90|Lhcnm2|Mcf7|Monocd14ro1746|Sknsh|Monocytescd14ro01746')
	CELLS3 = re.compile('Ag04449|Ag04450|Ag09309|Ag09319|Ag10803|Aoaf|Aosmc|Be2c|Bj|Caco2|Cmk|Dnd41|Ecc1|Gm06990|Gm12801|Gm12864|Gm12865|Gm12872|Gm12873|Gm12875|Gm12891|Gm12892|Gm19239|H7es|Hac|Hae|Hah|Hasp|Hbmec|Hcfaa|Hcf|Hcm|Hcpe|Hct116|Hee|Hek293|Hffmyc|Hff|Hgf|Hipe|Hl60|Hmec|Hmf|Hmvecdblad|Hnpce|Hpae|Hpaf|Hpdlf|Hpf|Hrce|Hre|Hrpe|Hsmmfshd|Hsmmtubefshd|Hsmmt|Hsmm|Htr8|Hvmf|Jurkat|Lncap|M059j|Mcf10aes|Nb4|Nha|Nhbe|Nhdfad|Nhdfneo|Nhek|Nhlf|Nt2d1|Osteobl|Osteo|Ovcar3|Panc1|Panislets|Pfsk1|Prec|Progfib|Rpmi7951|Rptec|Saec|Skmc|Sknmc|Sknshra|T47d|Th1|Th2|U87|Werirb1|Wi38|8988t|Cd34mobilized|Chorion|Cll|Fibrobl|Fibrop|Gliobla|Gm08714|Gm10847|Gm12874|Gm15510|Gm18505|Gm18507|Gm18526|Gm18951|Gm19099|Gm19193|Gm19238|Gm19240|H7hesc|H9es|Hconf|Hepatocytes|Hmvecdad|Hmvecdblneo|Hmvecdlyad|Hmvecdlyneo|Hmvecdneo|Hmveclbl|Hmveclly|Hpde6e6e7|Hrgec|Huh7|Huh75|Ips|Ishikawaestradiol|Ishikawatamoxifen|Medullo|Melano|Myometr|Panisletd|Pbde|Pbdefetal|Phte|Raji|Rwpe1|Shsy5y|Stellate|Th0|U2os|Urothelia|Urotheliaut189')
	m1 = CELLS1.search(line)
	m2 = CELLS2.search(line)
	m3 = CELLS3.search(line)
	if m1:
		Tier = 'Tier1'
	elif m2:
		Tier = 'Tier2'
	elif m3:
		Tier = 'Tier3'
	else:
		Tier = 'Tier3'
		with open(os.path.join(outputdir,"missing_tier.log"),'a') as writer:
			writer.write(line+"\n")

	return Tier

def preparebed_splitby(gf_outputdir,organism,gf_group, gf_file):
	''' A function that creates separate bed files for each value in field with index 'splitby'
	gf_outputdir: the directory in which the gf_folder should be created in which the split GFs should be outputted to
	EX: /[root]/grsnp_db/[organism]/[tier]/[celltype]/

	gf_file: file name with extension i.e gfname.bed.gz
	'''
	# TODO handle the case of partially finished database
	added_features = []
	# download the GF file	
	if "ftp_server" in gf_grp_sett[gf_group].keys():
		dwnl_file = download_encode_file(organism, gf_group, gf_file)
	elif "html_server" in gf_grp_sett[gf_group].keys():
		dwnl_file = download_roadmap_file(gf_group, gf_file)
	full_gf_paths = []
	# convert it to bed format
	logger.info( "Converting into proper bed format: {}".format(gf_file))
	if dwnl_file.endswith('.gz'):
		infile = gzip.open(dwnl_file)
	else:
		infile = open(dwnl_file)

	gf_file_ext = '.'.join(gf_file.split('.')[-2:]) # get file extension i.e. 'bed.gz'
	z_gf_file = gf_file.replace(".Z.","Z.") # special case: H2A.Z needs to be H2AZ
	output_gf_file = '_'.join(z_gf_file.split('.')[:-2]) # replace all other '.' with '_' for rest of filename
	out_gf_file = '.'.join([output_gf_file,gf_file_ext])

	file_writers = {} # {'cur_split_value': file_writer_object} A writer is created for each name field
	if not os.path.exists(gf_outputdir):
		os.makedirs(gf_outputdir)
	while True:
		line = infile.readline().rstrip('\n')
		if line == "":
			break
		cur_gf = preparebed_line[out_gf_file.replace(".gz",'').split(".")[-1]](line)
		cur_split_value = cur_gf[3]
		# check if current TFBS already has a file writer
		if cur_split_value not in file_writers:
			o_dir = gf_outputdir
			# get formated file name [cell]-[factor]-[source]
			eid =  _get_EID(gf_file,gf_group)
			tmp_name = cur_split_value
			if eid != "": # will be '' if working with ENCODE data
				tmp_name = eid + "-" + cur_split_value
			else:
				tmp_name = _get_celltype(gf_file,gf_group)+"-"+cur_split_value #add cell type to front
			# create formated file name
			form_dwnl_file = tmp_name + '-' + _get_source(gf_file, gf_group)

			out_path = os.path.join(o_dir,''.join(e for e in base_name(form_dwnl_file) if e.isalnum() or e=='.' or e=='_' or e=='-')) + ".bed.gz.temp"
			new_path = os.path.join(o_dir,''.join(e for e in base_name(form_dwnl_file) if e.isalnum() or e=='.' or e=='_' or e=='-')) + ".bed.gz"
			if os.path.exists(new_path):
				raise GF_ALREADY_EXISTS("{} already exists. Not going to overwrite it.".format(out_path))
			file_writers[cur_split_value] = open(out_path,'wb')
			full_gf_paths.append(out_path)
		file_writers[cur_split_value].write("\t".join(cur_gf)+"\n")
	#close all open files
	for k in file_writers.keys():
		file_writers[k].close()
	infile.close()
	# convert all created files into gzip
	converted_paths = []
	for f_path in full_gf_paths:
		o_dir = os.path.dirname(f_path)
		new_path = os.path.join(o_dir,''.join(e for e in base_name(f_path) if e.isalnum() or e=='.' or e=='_' or e=='-')) + ".bed.gz"
		sort_convert_to_bgzip(f_path,new_path)
		converted_paths.append(new_path)
	return converted_paths

def preparebed(gf_outputdir, organism, gf_group, gf_file):
	''' Converts the file to the correct bed format, sorts it, and gzips it.
	gf_outputdir: the directory in which the gf_folder should be created in which the split GFs should be outputted to
	EX: /[root]/grsnp_db/[organism]/[tier]/[celltype]/

	gf_file: file name with extension i.e gfname.bed.gz
	'''
	added_features = []
	if not os.path.exists(gf_outputdir): os.makedirs(gf_outputdir)
	# download the GF file
	if "ftp_server" in gf_grp_sett[gf_group].keys():
		dwnl_file = download_encode_file(organism, gf_group, gf_file)
	elif "html_server" in gf_grp_sett[gf_group].keys():
		dwnl_file = download_roadmap_file(gf_group,	gf_file)
	f_path = os.path.join(gf_outputdir,base_name(dwnl_file)+".bed.temp")
	o_dir = os.path.dirname(f_path)
	tmp_gf_file = dwnl_file.replace(".imputed.gappedPeak.bed.gPk.gz",".gPk.gz").replace(".imputed.narrowPeak.bed.nPk.gz",".nPk.gz")
	# replace all '.' with '_' except for the file extension portion
	gf_file_ext = '.'.join(tmp_gf_file.split('.')[-2:]) # get file extension i.e. 'bed.gz'
	z_gf_file = tmp_gf_file.replace(".Z.","Z.") # special case: H2A.Z needs to be H2AZ
	output_gf_file = '_'.join(z_gf_file.split('.')[:-2]) # replace all other '.' with '_' for rest of filename
	out_gf_file = '.'.join([output_gf_file,gf_file_ext])
	# get formated file name [cell]-[factor]-[source]
	form_dwnl_file = _get_formated_file_name(gf_group,os.path.split(dwnl_file)[1])

	new_path = os.path.join(o_dir,''.join(e for e in base_name(form_dwnl_file) if e.isalnum() or e=='.' or e=='_' or e=='-')) + ".bed.gz"
	if os.path.exists(new_path) == True:
		raise GF_ALREADY_EXISTS("{} already exists. Not going to overwrite it.".format(new_path))
	logger.info( "Converting into proper bed format: {}".format(gf_file))
	if gf_file.endswith('.gz'):
		infile = gzip.open(dwnl_file)
	else:
		infile = open(dwnl_file)
	# convert it to bed format	
	with open(f_path,'wb') as writer:
		while True:
			line = infile.readline().rstrip('\n')
			if line == "":
				break
			cur_gf = preparebed_line[gf_file.replace(".gz",'').split(".")[-1]](line)
			writer.write("\t".join(cur_gf)+"\n")
	# convert all created files into gzip
	sort_convert_to_bgzip(f_path,new_path)
	infile.close()
	return [new_path]


	

# Dictates how much of the filename to strip off before searching for cell type etc.
# Leave one word in front of cell name.  I.E. for  wgEncodeLicrHistoneBatH3k04me1MAdult24wksC57bl6StdAlnRep1
# padding should be 'wgEncodeLicr' in order to extract 'Bat' as cell type
padding = {
 "wgEncodeAwgTfbsUniform": "wgEncodeAwgTfbs",
 "wgEncodeRegTfbsClustered": "wgEncode",
 "wgEncodeBroadHmm": "wgEncodeBroad",
 "wgEncodeBroadHistone": "wgEncode",
 "wgEncodeUwHistone": "wgEncode",
 "wgEncodeSydhHistone": "wgEncode",
 "wgEncodeAwgDnaseUniform": "wgEncodeAwgDnase",
 "Encode_chromeStates": "wgEncodeAwgSegmentation",
 # mouse
 "wgEncodeCaltechHist": 'wgEncodeCaltech',
 "wgEncodeLicrHistone": 'wgEncodeLicr',
 "wgEncodePsuHistone": 'wgEncodePsu',
 "wgEncodeCaltechTfbs": 'wgEncodeCaltech',
 "wgEncodeLicrTfbs": 'wgEncodeLicr',
 "wgEncodePsuTfbs": 'wgEncodePsu',
 "wgEncodeSydhTfbs": 'wgEncodeSydh',
 "wgEncodeUwDnase": 'wgEncodeUw',
 "wgEncodeUwDgf": 'wgEncode',
 "wgEncodePsuDnase": 'wgEncodePsu' 
}

source = {
 "wgEncodeAwgTfbsUniform": "AwgTfbsUniform",
 'wgEncodeBroadHmm': "BroadHmm",
 "wgEncodeRegTfbsClustered": "TfbsV3",
 "wgEncodeBroadHistone": "BroadHistone",
 "wgEncodeUwHistone": "UwHistone",
 "wgEncodeSydhHistone": "SydhHistone",
 "wgEncodeAwgDnaseUniform": "AwgDnaseUniform",
 "chromStates15": "chromStates15",
 "chromStates18": "chromStates18",
 "chromStates25": "chromStates25",
 "Histone_processed_broadPeak": "processed",
 "Histone_processed_narrowPeak": "processed",
 "Histone_processed_gappedPeak": "processed",
 "Histone_imputed_narrowPeak": "imputed",
 "Histone_imputed_gappedPeak": "imputed",
 "DNase_processed_broadPeak": "processed",
 "DNase_processed_narrowPeak": "processed",
 "DNase_processed_gappedPeak": "processed",
 "DNase_imputed_narrowPeak": "imputed",
 "DNase_imputed_gappedPeak": "imputed",
 # mouse
 "wgEncodeCaltechHist": 'Caltech',
 "wgEncodeLicrHistone": 'Licr',
 "wgEncodePsuHistone": 'Psu',
 "wgEncodeCaltechTfbs": 'Caltech',
 "wgEncodeLicrTfbs": 'Licr',
 "wgEncodePsuTfbs": 'Psu',
 "wgEncodeSydhTfbs": 'Sydh',
 "wgEncodeUwDnase": 'Uwdnase',
 "wgEncodeUwDgf": 'Uwdgf',
 "wgEncodePsuDnase": 'Psu'
}


root_folder = {
"wgEncodeAwgTfbsUniform": "ENCODE/TFBS_cellspecific",
 "wgEncodeRegTfbsClustered": "ENCODE/TFBS_clustered",
 "wgEncodeBroadHmm": "ENCODE/chromStates/BroadHmm",
 "wgEncodeBroadHistone": "ENCODE/Histone",
 "wgEncodeUwHistone": "ENCODE/Histone",
 "wgEncodeSydhHistone": "ENCODE/Histone",
 "wgEncodeAwgDnaseUniform": "ENCODE/DNase",
 "chromStates15": "ROADMAP/chromStates15",
 'chromStates18':"ROADMAP/chromStates18",
 "chromStates25": "ROADMAP/chromStates25",
 "Histone_processed_broadPeak": "ROADMAP/Histone_bPk-processed",
 "Histone_processed_narrowPeak": "ROADMAP/Histone_nPk-processed",
 "Histone_processed_gappedPeak": "ROADMAP/Histone_gPk-processed",
 "Histone_imputed_narrowPeak": "ROADMAP/Histone_nPk-imputed",
 "Histone_imputed_gappedPeak": "ROADMAP/Histone_gPk-imputed",
 "DNase_processed_broadPeak": "ROADMAP/DNase_bPk-processed",
 "DNase_processed_narrowPeak": "ROADMAP/DNase_nPk-processed",
 "DNase_processed_gappedPeak": "ROADMAP/DNase_gPk-processed",
 "DNase_imputed_narrowPeak": "ROADMAP/DNase_nPk-imputed",
 "DNase_imputed_gappedPeak": "ROADMAP/DNase_gPk-imputed",
 "Encode_chromeStates": "ENCODE/chromStates",
 # mouse
 "wgEncodeCaltechHist": "ENCODE/Histone",
 "wgEncodeLicrHistone": "ENCODE/Histone",
 "wgEncodePsuHistone": "ENCODE/Histone",
 "wgEncodeCaltechTfbs": "ENCODE/TFBS",
 "wgEncodeLicrTfbs": "ENCODE/TFBS",
 "wgEncodePsuTfbs": "ENCODE/TFBS",
 "wgEncodeSydhTfbs": "ENCODE/TFBS" ,
 "wgEncodeUwDnase": "ENCODE/DNase",
 "wgEncodeUwDgf": "ENCODE/DNase",
 "wgEncodePsuDnase": "ENCODE/DNase" 
}

def _get_celltype(f_name, gf_group):
	''' Extracts the cell type from the ENCODE genomic feature file name
	Note that for the filepath, this function is NOT used (see _get_formated_file_name function)
	'''
	if f_name == "wgEncodeAwgDnaseDuke8988tUniPk":
		return "8988t" # special case since cell name starts with a number
	if gf_group == "wgEncodeRegTfbsClustered":
		return "Clustered"
	if f_name.startswith(padding[gf_group]):
		f_name = f_name[len(padding[gf_group]):]	

	# # these are special cases
	# if "K562b" in f_name:
	# 		cell_type = "K562"
	# 		f_name = f_name.replace("K562b","K562")
	# 		categories = re.findall('[A-Z][^A-Z]*', f_name.split('.')[0])
	# elif "K562E" in f_name and "K562Ezh2" not in f_name and "HaibK562" not in f_name and "SydhK562" not in f_name:
	# 	cell_type = "K562"
	# 	f_name = f_name.replace("K562Efos","K562Fos") # Capitalize factor letters
	# 	f_name = f_name.replace("K562Egata","K562Gata")
	# 	f_name = f_name.replace("K562Ehdac8","K562Hdac8")
	# 	f_name = f_name.replace("K562Ejun","K562Jun")
	# 	pdb.set_trace()
	# 	categories = re.findall('[A-Z][^A-Z]*', f_name.split('.')[0])		
	categories = re.findall('[A-Z][^A-Z]*', f_name.split('.')[0])
	return categories[1]

def _get_road_tissuegrp(f_name):
	'''Looks up tissue group using the E*** roadmap cell name'''
	# roadmapCellData is located in dbcreator_util.py
	road_conversion = [x.split("\t") for x in roadmapCellData.split("\n")[1:]]
	f_name_eid =  f_name.split("_")[0].split('-')[0] # need to split by both '_' and '-' since DNase uses '-' to separate the EID
	# [3] is the GR column in roadmapCellData
	tissue_group = [x for x in road_conversion if x[0] == f_name_eid][0][3]
	return tissue_group

def _get_EID(f_name, gf_group):
	'''Returns the roadmap EID of the current file'''
	if root_folder[gf_group].split('/')[0] == 'ENCODE':
		return ''
	eid = f_name.split("_")[0].split('-')[0] # need to split by both '_' and '-' since DNase uses '-' to separate the EID
	return eid

def _get_source(gf_name,gf_group):
	if gf_group == 'Encode_chromeStates':
		if gf_name.startswith(padding[gf_group]):
			gf_name = gf_name[len(padding[gf_group]):]
		categories = re.findall('[A-Z][^A-Z]*', gf_name)
		return categories[0]
	else:
		return source[gf_group]

def _get_formated_file_name(gf_group,gf_name):
	global organism
	rfold = root_folder[gf_group].split("/")[0] # 'ENCODE' or 'ROADMAP'
	peak_to_p = {'broadPeak': 'bPk', 'narrowPeak': 'nPk', 'gappedPeak': 'gPk'}
	if rfold == "ROADMAP":
			
		eid = _get_EID(gf_name,gf_group)
		factor = gf_name.replace(eid,'')[1:].split('.')[0]

		# append abbreviated peak to factor if it exists. i.e. 'narrowPeak' is 'nPk'
		peak = gf_group.split("_")[-1]
		if peak in peak_to_p.keys():
			factor = factor + "_" + peak_to_p[peak]
		return "-".join([eid, factor, _get_source(gf_name,gf_group)])
	elif rfold == "ENCODE":			
		if gf_name.startswith(padding[gf_group]):
			gf_name = gf_name[len(padding[gf_group]):]
		# Histone is special case in hg19 only, we append this to the gf_source later
		if "Histone" in gf_name and organism == "hg19":
			gf_name = gf_name.replace("Histone","")
			gf_source = "Histone"
		else: gf_source = ""
		# split filename by capital letters
		categories = re.findall('[A-Z][^A-Z]*', gf_name.split('.')[0])
		# get the cell type
		cell_type = categories[1]
		# special cell types
		if gf_group in ['wgEncodeUwDnase', 'wgEncodeUwDgf']:
			cell_type =  gf_name.split('Hotspot')[0].replace("Dnase",'').replace("Dgf","").replace("Uw",'')
		if gf_group in ['wgEncodePsuDnase']:
			cell_type =  gf_name.split('Pk')[0].replace('Dnase','').replace("Dgf","")
		if "8988t" in gf_name: # cell name that doesn't start with capital letter
			cell_type = "8988t"
			gf_name = gf_name.replace("8988t","Abc") # replace with dummy cell_type so splitting by letters works
			categories = re.findall('[A-Z][^A-Z]*', gf_name.split('.')[0])
		elif "K562b" in gf_name:
			cell_type = "K562"
			gf_name = gf_name.replace("K562b","K562")
			categories = re.findall('[A-Z][^A-Z]*', gf_name.split('.')[0])
		elif "K562E" in gf_name and "K562Ezh2" not in gf_name:
			cell_type = "K562"
			gf_name = gf_name.replace("K562E","K562")
			categories = re.findall('[A-Z][^A-Z]*', gf_name.split('.')[0])
		
		# get gf_source
		gf_source = categories[0] + gf_source
		if organism == 'mm9':
			gf_source = source[gf_group]
		# get factor
		if gf_group in ["wgEncodeUwDnase","wgEncodecAwgDnaseUniform","wgEncodePsuDnase"]:
			factor = 'DNase'
			if gf_group in ["wgEncodeUwDnase"]:
				# appends the 'Rep' portion of the file name to the factor
				if "Rep" in gf_name:
					factor = 'DNase' + '_Rep'+gf_name.split(".")[0].split("Rep")[1]
		elif gf_group in ["wgEncodeUwDgf"]:
			factor = "Dgf"
		else:
			if gf_group in ["wgEncodePsuHistone","wgEncodeCaltechTfbs" ,"wgEncodePsuTfbs"]:
				if 'PeaksRep1' in gf_name or 'PkRep1' in gf_name:
					factor = categories[2] + "_rep1"
				elif 'PeaksRep2' in gf_name or 'PkRep2' in gf_name:
					factor = categories[2] + '_rep2'
				else:
					factor = categories[2]			
			else:
				factor = categories[2]
		return "-".join([cell_type, factor, gf_source])



def _get_gf_directory(outputdir,gf_group,gf_name):
	''' Returns the output_dir of the gf.
	'''

	# Dictates the structure of the directory
	dirstruture = {
	 "wgEncodeAwgTfbsUniform": ['tier','cell'],
	 'wgEncodeBroadHmm': ['cell'],
	 "wgEncodeRegTfbsClustered": [],
	 "wgEncodeBroadHistone": ['tier','cell'],
	 "wgEncodeUwHistone": ['tier','cell'],
	 "wgEncodeSydhHistone": ['tier','cell'],
	 "wgEncodeAwgDnaseUniform": ['tier','cell'],
	 "chromStates15": ['roadmap_tissue','EID'],
	 "chromStates18": ['roadmap_tissue','EID'],
	 "chromStates25": ['roadmap_tissue','EID'],
	 "Histone_processed_broadPeak": ['roadmap_tissue','EID'],
	 "Histone_processed_narrowPeak": ['roadmap_tissue','EID'],
	 "Histone_processed_gappedPeak": ['roadmap_tissue','EID'],
	 "Histone_imputed_narrowPeak": ['roadmap_tissue','EID'],
	 "Histone_imputed_gappedPeak": ['roadmap_tissue','EID'],
	 "DNase_processed_broadPeak": ['roadmap_tissue','EID'],
	 "DNase_processed_narrowPeak": ['roadmap_tissue','EID'],
	 "DNase_processed_gappedPeak": ['roadmap_tissue','EID'],
	 "DNase_imputed_narrowPeak": ['roadmap_tissue','EID'],
	 "DNase_imputed_gappedPeak": ['roadmap_tissue','EID'],
	 "Encode_chromeStates": ['source','cell'],
	 # mouse
	 "wgEncodeCaltechHist": ['source'],
	 "wgEncodeLicrHistone": ['source'],
	 "wgEncodePsuHistone": ['source'],
	 "wgEncodeCaltechTfbs": ['source'],
	 "wgEncodeLicrTfbs": ['source'],
	 "wgEncodePsuTfbs": ['source'],
	 "wgEncodeSydhTfbs": ['source'],
	 "wgEncodeUwDnase": ['source'],
	 "wgEncodeUwDgf": ['source'],
	 "wgEncodePsuDnase": ['source']
	}

	gf_name = gf_name.replace("Histone","")
	dir_structure = dirstruture[gf_group]
	gf_directory = [root_folder[gf_group]]
	for folder in dirstruture[gf_group]:
		if folder == 'tier':
			gf_directory.append(_get_tier(gf_name,outputdir))
		elif folder == 'cell':
			gf_directory.append(_get_celltype(gf_name,gf_group))
		elif folder == 'roadmap_tissue':
			gf_directory.append(_get_road_tissuegrp(gf_name))
		elif folder == 'EID':
			gf_directory.append(_get_EID(gf_name,gf_group))
		elif folder == 'source':
			gf_directory.append(_get_source(gf_name,gf_group))
	gf_directory = "/".join(gf_directory)
	return os.path.join(outputdir,gf_directory)

def create_feature_set(data_dir,organism,gf_group,max_install = None):
	outputdir = os.path.join(data_dir,organism)	
	added_features = [] 
	outpath = ""
	summary_path = os.path.join(outputdir,"summary.log")
	if not os.path.exists(outputdir): os.makedirs(outputdir)
	open(summary_path,'wb')


	grp_count,gf_file_paths = 0,[]
	if "gf_files" in gf_grp_sett[gf_group].keys():
		gf_file_paths = gf_grp_sett[gf_group]["gf_files"]
		url = gf_grp_sett[gf_group]['ftp_server'] +  gf_grp_sett[gf_group]['directory'].format(organism)
	else:
		url = ""
		if "ftp_server" in gf_grp_sett[gf_group].keys():
			gf_file_paths = get_gf_filepaths(organism,gf_group)
			url = gf_grp_sett[gf_group]['ftp_server'] +  gf_grp_sett[gf_group]['directory'].format(organism)
		elif "html_server" in gf_grp_sett[gf_group].keys():
			gf_file_paths = get_road_gf_filepaths(gf_group)
			url = gf_grp_sett[gf_group]['html_server'] +  gf_grp_sett[gf_group]['directory'].format(organism)
	for gf_file in gf_file_paths:		
		# check if gf_type is supported
		if gf_file.replace(".gz",'').split(".")[-1] not in preparebed_line.keys():
			continue
		# limit the number of GFs to install per group	
		if max_install != None and grp_count >= max_install:
			break
		try:
			grp_count += 1			
			gf_type = ""
			gf_outputdir = _get_gf_directory(outputdir,gf_group,base_name(gf_file))
			if not os.path.exists(gf_outputdir):
				os.makedirs(gf_outputdir)
			# removes the .temp file, to prevent duplicate data from being written
			if os.path.exists(outpath+".temp"):
				os.remove(outpath+".temp")
			try:
				# converts the ucsc data into proper bed format		
				gf_paths = gf_grp_sett[gf_group]["prep_method"](gf_outputdir,organism,gf_group,gf_file)
				# output description
				for f in gf_paths:
					# get the relative full path i.e. grsnp/ENCODE/....
					relative_outputdir = "/grsnp_db/" + os.path.split(f)[0].split("/grsnp_db")[-1]
					rel_f = os.path.join(relative_outputdir, os.path.split(f)[-1]) 
					if not rel_f in gf_descriptions.keys():
						f_parts = base_name(f).split("-")
						gf_descriptions[rel_f] = GF_description(**{"file_name": base_name(f), "full_path": rel_f, "URL": url+"/" + gf_file,
							"full_description": "", "category": gf_group, "category_desc": "", "cell": f_parts[0], "cell_desc": "",	
							"factor": f_parts[1], "factor_desc": "", "source": f_parts[2], "source_desc": ""})
				_write_description_file(data_dir,organism)
			except GF_ALREADY_EXISTS:
				logger.info("{} already exists, skipping extraction".format(outpath.replace(".gz","")))
				continue
			
			# cannot detect type, skip
			if gf_type == "failed":
				write_line("\t".join([gf_file,"Not supported","None"]),summary_path)
				logger.warning( "Unable to convert {} into bed".format(gf_file))
				continue
			added_features.append(outpath)
			write_line("\t".join([gf_file,gf_type,"None"]),summary_path)

		except Exception, e:
			write_line("\t".join([gf_file,"Failed",str(e)]),summary_path)
			exc = trace.format_exc()
			logger.warning( "Unable to convert {} into bed".format(gf_file))
			logger.warning(exc)
			continue	

	# cleanup the temporary files
	if os.path.exists(outpath + ".temp"): os.remove(outpath+".temp")

	# logger.info( "The following types are not supported (includes all 'big' file types):\n " + str(notsuptypes))
	# logger.info("The following features were added to the database: \n{}".format(added_features))
	return "Created UCSC database"

def _write_description_file(data_dir,organism):
	# write gf_description to temp file which is renamed after writing is complete
	outpath = os.path.join(data_dir,organism,"gf_descriptions.txt")
	with open(outpath + ".tmp",'wb') as writer:
		writer.write("\t".join(gf_description_colnames)+"\n")
		for full_path,gf_desc in gf_descriptions.iteritems():
			writer.write("\t".join([str(i) for i in gf_desc])+"\n")
	if os.path.exists(outpath):
		os.remove(outpath)
	os.rename(outpath+".tmp",outpath)

def _read_description_file(data_dir,organism):
	descriptions = {}
	outpath = os.path.join(data_dir,organism,"gf_descriptions.txt")
	if not os.path.exists(outpath):
		return descriptions
	with open(outpath) as reader:
		while True:			
			line = reader.readline().rstrip('\n')
			if line == "":
				break
			row = GF_description(**dict(zip(gf_description_colnames, line.split("\t"))))
			descriptions[row.full_path] = row
	return descriptions

def update_progress(progress):
    print '\r[{0}] {1}%'.format('#'*(progress/10), progress),


# each genomic feature group must have an entry in this dictionary
# f_names is optional and can be left out.  It limits the download of the gf group to only the files listed
# 'ftp_server' and 'directory' are used to download files from the database. 
# '{}' in 'directory' is filled in with the organism name when present
gf_grp_sett = {
	# human
	"wgEncodeAwgTfbsUniform": {'prep_method': preparebed,
				'ftp_server': 'hgdownload.cse.ucsc.edu', 'directory': '/goldenPath/{}/encodeDCC/wgEncodeAwgTfbsUniform'},
	"wgEncodeRegTfbsClustered": {"prep_method":  preparebed_splitby, "gf_files": ["wgEncodeRegTfbsClusteredV3.bed.gz"],
				"ftp_server": 'hgdownload.cse.ucsc.edu', "directory": '/goldenPath/{}/encodeDCC/wgEncodeRegTfbsClustered'},
	"wgEncodeBroadHmm": {"prep_method":  preparebed_splitby,
				"ftp_server": 'hgdownload.cse.ucsc.edu', "directory": '/goldenPath/{}/encodeDCC/wgEncodeBroadHmm'},
	"wgEncodeBroadHistone": {"prep_method":  preparebed,
				"ftp_server": 'hgdownload.cse.ucsc.edu', "directory": '/goldenPath/{}/encodeDCC/wgEncodeBroadHistone'},
	"wgEncodeUwHistone": {"prep_method":  preparebed,
				"ftp_server": 'hgdownload.cse.ucsc.edu', "directory": '/goldenPath/{}/encodeDCC/wgEncodeUwHistone'},
	"wgEncodeSydhHistone": {"prep_method":  preparebed,
				"ftp_server": 'hgdownload.cse.ucsc.edu', "directory": '/goldenPath/{}/encodeDCC/wgEncodeSydhHistone'},
	"wgEncodeAwgDnaseUniform": {"prep_method":  preparebed,
				"ftp_server": 'hgdownload.cse.ucsc.edu', "directory": '/goldenPath/{}/encodeDCC/wgEncodeAwgDnaseUniform'},
	"chromStates15": {'prep_method': preparebed_splitby,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final"},
	"chromStates18": {'prep_method': preparebed_splitby,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final"},
	"chromStates25": {'prep_method': preparebed_splitby,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final"},
	"Histone_processed_broadPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidated/broadPeak"},
	"Histone_processed_narrowPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidated/narrowPeak"},
	"Histone_processed_gappedPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidated/gappedPeak"},
	"Histone_imputed_narrowPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak"},
	"Histone_imputed_gappedPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidatedImputed/gappedPeak/"},
	"DNase_processed_broadPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidated/broadPeak"},
	"DNase_processed_narrowPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidated/narrowPeak"},
	"DNase_processed_gappedPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidated/gappedPeak"},
	"DNase_imputed_narrowPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak"},
	"DNase_imputed_gappedPeak": {'prep_method': preparebed,
	 			"html_server": 'http://egg2.wustl.edu', 'directory': "/roadmap/data/byFileType/peaks/consolidatedImputed/gappedPeak"},
	"Encode_chromeStates": {'prep_method': preparebed_splitby,
				'ftp_server': 'hgdownload.cse.ucsc.edu', 'directory': "/goldenPath/{}/encodeDCC/wgEncodeAwgSegmentation"},
	# mouse
	"wgEncodeCaltechHist": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodeCaltechHist"},
	"wgEncodeLicrHistone": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone"},
	"wgEncodePsuHistone": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodePsuHistone"},
	"wgEncodeCaltechTfbs": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodeCaltechTfbs"},
	"wgEncodeLicrTfbs": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodeLicrTfbs"},
	"wgEncodePsuTfbs": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodePsuTfbs"},
	"wgEncodeSydhTfbs": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodeSydhTfbs"},
	"wgEncodeUwDnase": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodeUwDnase"},
	"wgEncodeUwDgf": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodeUwDgf"},
	"wgEncodePsuDnase": {'prep_method': preparebed,
				'ftp_server': "hgdownload.cse.ucsc.edu", 'directory': "/goldenPath/mm9/encodeDCC/wgEncodePsuDnase"},



}

org_gfgroup = {
	'hg19': [
		"wgEncodeAwgTfbsUniform",
		"wgEncodeRegTfbsClustered",
		"wgEncodeBroadHmm",
		"wgEncodeBroadHistone",
		"wgEncodeUwHistone",
		"wgEncodeSydhHistone",
		"wgEncodeAwgDnaseUniform",
		"chromStates15",
		"chromStates18",
		"chromStates25",
		"Histone_processed_broadPeak",
		"Histone_processed_narrowPeak",
		"Histone_processed_gappedPeak",
		"Histone_imputed_narrowPeak",
		"Histone_imputed_gappedPeak",
		"DNase_processed_broadPeak",
		"DNase_processed_narrowPeak",
		"DNase_processed_gappedPeak",
		"DNase_imputed_narrowPeak",
		"DNase_imputed_gappedPeak",
		"Encode_chromeStates",
	],
	'mm9':[
		"wgEncodeCaltechHist",
		"wgEncodeLicrHistone",
		"wgEncodePsuHistone",
		"wgEncodeCaltechTfbs",
		"wgEncodeLicrTfbs",
		"wgEncodePsuTfbs",
		"wgEncodeSydhTfbs",
		"wgEncodeUwDnase",
		"wgEncodeUwDgf",
		"wgEncodePsuDnase"
	]
}

# Dictionary of supported file type mapped to function used to convert to proper bed format
preparebed_line = {
	"bed":_format_bed,
	"bed9":_format_bed,
	"bedRrbs":_format_bed,
	"peptideMapping":_format_bed,
	"broadPeak":_format_peak,
	"narrowPeak":_format_peak,
	"gappedPeak": _format_gapped_peak,
	"bedRnaElements":_format_bed,
	"nPk": _format_bed,
	"gPk": _format_bed
}



def main():
	global username,password, download_dir, gfs, organism
	parser = argparse.ArgumentParser(prog="python -m grsnp.dbcreator", description='Creates the GenomeRunner SNP Database. Example: python -m grsnp.dbcreator -d /home/username/grsnp_db/ -g mm9', epilog='IMPORTANT: Execute DBCreator from the database folder, e.g., /home/username/grsnp_db/. Downloaded files from UCSC are placed in ./downloads database created in ./grsnp_db.')
	parser.add_argument("--data_dir" , "-d", nargs="?", help="Set the directory where the database to be created. Use absolute path. Example: /home/username/grsnp_db/. Required", required=True)
	parser.add_argument('--organism','-g', nargs="?", help="The UCSC code of the organism to use for the database creation. Default: hg19 (human). Required", default="hg19")
	parser.add_argument('--featuregroups','-f', nargs="?", help='The names of the specific genomic feature groups to download.  List available for hg19 at ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/', default="")
	parser.add_argument('--galaxy', help="Create the xml files needed for Galaxy. Outputted to the current working directory.", action="store_true")
	parser.add_argument('--quantiles', '-s', help="Commas separated list of score quantiles.", nargs='?',default="")
	#parser.add_argument('--filteronly','-o', help="Only filter by score and strand. Skips downloading and installing new GFs.", action="store_true")
	parser.add_argument('--max','-m',help="Limit the number of features to be created within each group.",type=int, default=None)



	args = vars(parser.parse_args())

	hdlr = logging.FileHandler(os.path.join(args['data_dir'], 'genomerunner_dbcreator.log'))
	hdlr_std = StreamHandler()
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.addHandler(hdlr_std)
	logger.setLevel(logging.INFO)

	if not args["data_dir"]:
		print "ERROR: --data_dir is required"
		sys.exit()
	if args['quantiles'] == "":
		quantiles = [25,50,75]
	else:
		quantiles = set(map(int,args['quantiles'].split(','))) # remove duplicate scores
		quantiles = sorted(list(quantiles))
	data_dir=os.path.join(args["data_dir"],'grsnp_db')

	if args['organism'] is not None: # Only organism is specified. Download all organism-specific features
		global download_dir, gf_grp_sett
		organism = args['organism']
		download_dir = os.path.join(args["data_dir"],"downloads",args['organism'])
		gfs = args["featuregroups"].split(",")
		gf_descriptions = _read_description_file(data_dir,args["organism"])
		for grp in org_gfgroup[args["organism"]]:
			create_feature_set(data_dir,args['organism'],grp,args['max'])
	else:
		print "ERROR: Requires UCSC organism code.  Use --help for more information"
		sys.exit()

	### Second Step: Create subdirectories for score and filter data by score quantiles
	# create sub directories for score percentiles and populate with score-filtered GF data
	# gather all directories (groups) in the database
	print "Filtering GFs by strand and score..."
	orgdir = os.path.join(data_dir,args['organism'])
	dirs = [name for name in os.listdir(orgdir)
		if os.path.isdir(os.path.join(orgdir, name))]
	for d in dirs:
		# gather all paths
		gfs = []
		for base, tmp, files in os.walk(os.path.join(orgdir,d)):
				gfs += [os.path.join(base,f) for f 	 
					in files if f.endswith(('.gz'))]
		for gf_path in gfs:			
			print "Filtering {} ...".format(gf_path)
			# filter the original GF by strand. filter_by_strand checks if this step has been done
			strand_paths = filter_by_strand(data_dir,gf_path)
			# filter original and strand filtered files by score
			for p in [gf_path] + strand_paths:
				# check if score column exists
				with gzip.open(p) as dr:
					line = dr.readline().rstrip('\n')
					if len(line.split('\t')) < 5:
						print "{} lacks score column..."
						continue
				# check if the file has already been processed.
				# Done before loading into numpy array since that takes time
				scores_processed = True
				for i in range(len(quantiles)):
					gf_scorepath_out = p.replace('/grsnp_db','/grsnp_db_{}'.format(quantiles[i]))
					if os.path.exists(gf_scorepath_out):
						scores_processed = False
				if not scores_processed:
					continue
				# read in score column into numpy array
				dt = np.loadtxt(p, usecols=(4,), delimiter  = '\t') 
				score_min,score_max = np.min(dt),np.max(dt)
				if score_max == score_min:
					print "Min_score == Max_score == {} ...skipping".format(np.round(score_min,2))
					continue
				# calculate the quantiles using numpy
				quantile_thresh = np.percentile(dt,quantiles)
				# filter by the quantile thresholds
				for i in range(len(quantile_thresh)):
					gf_scorepath_out = p.replace('/grsnp_db','/grsnp_db_{}'.format(quantiles[i]))					
					if not os.path.exists(os.path.split(gf_scorepath_out)[0]):
						os.makedirs(os.path.split(gf_scorepath_out)[0])
					# filepath without extension is needed by the filter function
					gf_path_out_woext = os.path.join(os.path.split(gf_scorepath_out)[0],base_name(gf_scorepath_out))
					# filter by score
					filter_by_score(p, gf_path_out_woext,quantile_thresh[i])
				logger.info("MinMax stats for {}: Min={}, Max={}, quantile values={}".format(base_name(p), score_min,score_max,np.round(quantile_thresh,2)))

	root_dir = os.path.dirname(os.path.realpath(__file__))
	readme = open(os.path.join(root_dir,"grsnp_db_readme.txt")).read()
	with open("grsnp_db_readme.txt","wb") as writer:
		writer.write(readme)
	print "FINISHED: Downloaded files from UCSC are placed in {}.  Database created in {}".format(os.path.join(args["data_dir"],"downloads"),os.path.join(args["data_dir"],"grsnp_db`"))



if __name__ == "__main__":
	main()


