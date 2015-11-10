#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os
import subprocess
import pdb
import gzip


def remove_headers(datapath):
	''' Removes all lines that start with special characters, such as header lines.
	Replaces the datapath file with a file that does not contain header lines.
	Handles both gziped files and raw text (based on file extension)
	'''
	outputpath = datapath + ".tmp"
	if datapath.endswith('.gz'):
		infile,outfile = gzip.open(datapath),gzip.open(outputpath,'wb')
	else:
		infile,outfile = open(datapath),open(outputpath,'wb')
	while True:
		line = infile.readline().rstrip('\n')
		if line == "":
			break
		if line[:5] == 'track':
			continue
		if line[0].isalnum():
			outfile.write(line+"\n")
	infile.close(),outfile.close()
	os.remove(datapath)
	os.rename(outputpath,datapath)

def load_minmax(path):
	data = {}
	if not os.path.exists(path):
		return data
	score = [x for x in open(path).read().split("\n") if x != ""]
	for s in score:
		name,min_max = s.split('\t')
		data[name] = min_max
	return data


def save_minmax(data,path):
	''' Saves the dictionary of key value pairs of minmax data to a text file.
	'''
	with open(path,'wb') as writer:
		for k,v in data.items():
			writer.write("{}\t{}\n".format(k,v))


def base_name(k):
	return os.path.basename(k).split(".")[0]

def sort_convert_to_bgzip(path,outpath):
	script = "sort -k1,1 -k2,2n -k3,3n " + path +" | bgzip -c > " + outpath + ".gz.temp"
	out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE)
	out.wait()
	os.remove(path)# remove the .temp file extension to activate the GF		
	os.rename(outpath+".gz.temp",outpath)
	script = "tabix -f " + outpath
	out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE)
	out.wait()



def filter_by_score(gf_path_input,gf_path_output,thresh_score):
	''' Read in the gf data from gf_path_input and filter out each GF that does not
	have a score greater than the thresh_score threshold.
	gf_path_output should be WITHOUT file extension
	'''

	tmp_path = gf_path_output+'.temp'
	# check if score column exists
	with gzip.open(gf_path_input) as dr:
		line = dr.readline()
		if len(line.split("\t")) < 5:
			return []
	with open(tmp_path,"wb") as bed:
		with gzip.open(gf_path_input) as dr:
			while True:
				line = dr.readline().rstrip('\n')
				cur_gf = line.split('\t')
				if line == "":
					break
				score  = cur_gf[4]
				# if the score is >= to the threshold, output that GF
				if float(score) >= float(thresh_score):
					bed.write(line+"\n")
	sort_convert_to_bgzip(tmp_path,gf_path_output+'.bed.gz')
	return gf_path_output + '.bed.gz'

def filter_by_strand(data_dir,gf_path):
	''' Read in the gf data from gf_path and create new GFs for each strand.
	gf_path should be the full path to the gf WITHOUT file extension.
	data_dir should be the directory containing the database + grsnp_db_[score]
	'''
	sub_data_path = gf_path.replace(data_dir+"/","")
	plus_path, minus_path = os.path.join(data_dir+"_plus/",sub_data_path+'.temp'),os.path.join(data_dir+"_minus/",sub_data_path+'.temp')
	if not os.path.exists(os.path.split(plus_path)[0]): os.makedirs(os.path.split(plus_path)[0])
	if not os.path.exists(os.path.split(minus_path)[0]): os.makedirs(os.path.split(minus_path)[0])
	final_plus_path = os.path.join(os.path.split(plus_path)[0],base_name(plus_path)+'.bed.gz')
	final_minus_path = os.path.join(os.path.split(minus_path)[0],base_name(minus_path)+'.bed.gz')
	if os.path.exists(final_plus_path) and os.path.exists(final_minus_path):
		pdb.set_trace()
		return [final_plus_path,final_minus_path]
	plus_file,minus_file = open(plus_path,'wb'),open(minus_path,'wb')
	strand_file = {"+":plus_file,'-':minus_file}
	# check if strand column exists
	with gzip.open(gf_path) as dr:
		line = dr.readline()
		if len(line.split("\t")) < 6:
			return []
	with gzip.open(gf_path) as dr:		
		while True:
			line = dr.readline().rstrip('\n')
			cur_gf =  line.split('\t')			
			if line == "":
				break
			strand  = cur_gf[5]
			# check if a valid strand exists and output to the appropriate file.
			if strand in strand_file:
				strand_file[strand].write(line+"\n")
	plus_file.close()
	minus_file.close()
	# remove the  strand filtered gf file if empty
	strand_paths = []
	if os.stat(plus_path).st_size==0:
		os.remove(plus_path)
	else:
		sort_convert_to_bgzip(plus_path,final_plus_path)
		strand_paths.append(final_plus_path)
	if os.stat(minus_path).st_size==0:
		os.remove(minus_path)
	else:
		sort_convert_to_bgzip(minus_path,final_minus_path)
		strand_paths.append(final_minus_path)
	return strand_paths


class MinMax:
	"""Used to keep track of the min and max score. 
	By Default the min and max score is None.  It is sufficient to create this class for GFs that lack a score field.
	The class will return an NA string if the score is never updated"""
	def __init__(self):
		self.max,self.min = None, None

	def update_minmax(self,n):
		
		try:
			n = float(n)		
		except ValueError:
		    return
		# Assign the first value found to min and max
		if self.max == None:
			self.max = n
		if self.min == None:
			self.min = n
		# Check if we have found a new min or max
		if n < self.min:
			self.min = n
		elif n > self.max:
			self.max = n
	
	def str_minmax(self):
		n_max, n_min = 'NA','NA'
		if self.max != None:
			n_max = self.max
		if self.min != None:
			n_min = self.min
		if n_min == n_max:
			n_max, n_min = 'NA','NA'
		return "{},{}".format(n_min,n_max)

def write_line(line,path):
	with open(path, 'a') as writer:
		writer.write(line+"\n")





roadmapCellData ="""EID	GROUP	ANATOMY	GR	Epigenome Mnemonic	Standardized Epigenome name (from Wouter)	Epigenome name (from EDACC Release 9 directory)	TYPE
E033	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD3.CPC	Primary T cells from cord blood	CD3_Primary_Cells_Cord_BI	PrimaryCell
E034	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD3.PPC	Primary T cells from peripheral blood	CD3_Primary_Cells_Peripheral_UW	PrimaryCell
E037	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.MPC	Primary T helper memory cells from peripheral blood 2	CD4_Memory_Primary_Cells	PrimaryCell
E038	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.NPC	Primary T helper naive cells from peripheral blood	CD4_Naive_Primary_Cells	PrimaryCell
E039	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.CD25M.CD45RA.NPC	Primary T helper naive cells from peripheral blood	CD4+_CD25-_CD45RA+_Naive_Primary_Cells	PrimaryCell
E040	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.CD25M.CD45RO.MPC	Primary T helper memory cells from peripheral blood 1	CD4+_CD25-_CD45RO+_Memory_Primary_Cells	PrimaryCell
E041	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.CD25M.IL17M.PL.TPC	Primary T helper cells PMA-I stimulated	CD4+_CD25-_IL17-_PMA-Ionomycin_stimulated_MACS_purified_Th_Primary_Cells	PrimaryCell
E042	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.CD25M.IL17P.PL.TPC	Primary T helper 17 cells PMA-I stimulated	CD4+_CD25-_IL17+_PMA-Ionomcyin_stimulated_Th17_Primary_Cells	PrimaryCell
E043	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.CD25M.TPC	Primary T helper cells from peripheral blood	CD4+_CD25-_Th_Primary_Cells	PrimaryCell
E044	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.CD25.CD127M.TREGPC	Primary T regulatory cells from peripheral blood	CD4+_CD25+_CD127-_Treg_Primary_Cells	PrimaryCell
E045	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD4.CD25I.CD127.TMEMPC	Primary T cells effector/memory enriched from peripheral blood	CD4+_CD25int_CD127+_Tmem_Primary_Cells	PrimaryCell
E047	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD8.NPC	Primary T CD8+ naive cells from peripheral blood	CD8_Naive_Primary_Cells	PrimaryCell
E048	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.CD8.MPC	Primary T CD8+ memory cells from peripheral blood	CD8_Memory_Primary_Cells	PrimaryCell
E062	Blood & T-cell	BLOOD	Blood_and_Tcell	BLD.PER.MONUC.PC	Primary mononuclear cells from peripheral blood	Peripheral_Blood_Mononuclear_Primary_Cells	PrimaryCell
E053	Neurosph	BRAIN	Brain	BRN.CRTX.DR.NRSPHR	Cortex derived primary cultured neurospheres	Neurosphere_Cultured_Cells_Cortex_Derived	PrimaryCulture
E054	Neurosph	BRAIN	Brain	BRN.GANGEM.DR.NRSPHR	Ganglion Eminence derived primary cultured neurospheres	Neurosphere_Cultured_Cells_Ganglionic_Eminence_Derived	PrimaryCulture
E067	Brain	BRAIN	Brain	BRN.ANG.GYR	Brain Angular Gyrus	Brain_Angular_Gyrus	PrimaryTissue
E068	Brain	BRAIN	Brain	BRN.ANT.CAUD	Brain Anterior Caudate	Brain_Anterior_Caudate	PrimaryTissue
E069	Brain	BRAIN	Brain	BRN.CING.GYR	Brain Cingulate Gyrus	Brain_Cingulate_Gyrus	PrimaryTissue
E070	Brain	BRAIN	Brain	BRN.GRM.MTRX	Brain Germinal Matrix	Brain_Germinal_Matrix	PrimaryTissue
E071	Brain	BRAIN	Brain	BRN.HIPP.MID	Brain Hippocampus Middle	Brain_Hippocampus_Middle	PrimaryTissue
E072	Brain	BRAIN	Brain	BRN.INF.TMP	Brain Inferior Temporal Lobe	Brain_Inferior_Temporal_Lobe	PrimaryTissue
E073	Brain	BRAIN	Brain	BRN.DL.PRFRNTL.CRTX	Brain_Dorsolateral_Prefrontal_Cortex	Brain_Mid_Frontal_Lobe	PrimaryTissue
E074	Brain	BRAIN	Brain	BRN.SUB.NIG	Brain Substantia Nigra	Brain_Substantia_Nigra	PrimaryTissue
E081	Brain	BRAIN	Brain	BRN.FET.M	Fetal Brain Male	Fetal_Brain_Male	PrimaryTissue
E082	Brain	BRAIN	Brain	BRN.FET.F	Fetal Brain Female	Fetal_Brain_Female	PrimaryTissue
E075	Digestive	GI_COLON	Digestive	GI.CLN.MUC	Colonic Mucosa	Colonic_Mucosa	PrimaryTissue
E077	Digestive	GI_DUODENUM	Digestive	GI.DUO.MUC	Duodenum Mucosa	Duodenum_Mucosa	PrimaryTissue
E079	Digestive	GI_ESOPHAGUS	Digestive	GI.ESO	Esophagus	Esophagus	PrimaryTissue
E084	Digestive	GI_INTESTINE	Digestive	GI.L.INT.FET	Fetal Intestine Large	Fetal_Intestine_Large	PrimaryTissue
E085	Digestive	GI_INTESTINE	Digestive	GI.S.INT.FET	Fetal Intestine Small	Fetal_Intestine_Small	PrimaryTissue
E092	Digestive	GI_STOMACH	Digestive	GI.STMC.FET	Fetal Stomach	Fetal_Stomach	PrimaryTissue
E094	Digestive	GI_STOMACH	Digestive	GI.STMC.GAST	Gastric	Gastric	PrimaryTissue
E101	Digestive	GI_RECTUM	Digestive	GI.RECT.MUC.29	Rectal Mucosa Donor 29	Rectal_Mucosa.Donor_29	PrimaryTissue
E102	Digestive	GI_RECTUM	Digestive	GI.RECT.MUC.31	Rectal Mucosa Donor 31	Rectal_Mucosa.Donor_31	PrimaryTissue
E106	Digestive	GI_COLON	Digestive	GI.CLN.SIG	Sigmoid Colon	Sigmoid_Colon	PrimaryTissue
E109	Digestive	GI_INTESTINE	Digestive	GI.S.INT	Small Intestine	Small_Intestine	PrimaryTissue
E110	Digestive	GI_STOMACH	Digestive	GI.STMC.MUC	Stomach Mucosa	Stomach_Mucosa	PrimaryTissue
E114	ENCODE2012	LUNG	ENCODE2012	LNG.A549.ETOH002.CNCR	A549 EtOH 0.02pct Lung Carcinoma Cell Line	A549_EtOH_0.02pct_Lung_Carcinoma	CellLine
E115	ENCODE2012	BLOOD	ENCODE2012	BLD.DND41.CNCR	Dnd41 TCell Leukemia Cell Line	Dnd41_TCell_Leukemia	CellLine
E116	ENCODE2012	BLOOD	ENCODE2012	BLD.GM12878	GM12878 Lymphoblastoid Cells	GM12878_Lymphoblastoid	PrimaryCulture
E117	ENCODE2012	CERVIX	ENCODE2012	CRVX.HELAS3.CNCR	HeLa-S3 Cervical Carcinoma Cell Line	HeLa-S3_Cervical_Carcinoma	CellLine
E118	ENCODE2012	LIVER	ENCODE2012	LIV.HEPG2.CNCR	HepG2 Hepatocellular Carcinoma Cell Line	HepG2_Hepatocellular_Carcinoma	CellLine
E119	ENCODE2012	BREAST	ENCODE2012	BRST.HMEC	HMEC Mammary Epithelial Primary Cells	HMEC_Mammary_Epithelial	PrimaryCulture
E120	ENCODE2012	MUSCLE	ENCODE2012	MUS.HSMM	HSMM Skeletal Muscle Myoblasts Cells	HSMM_Skeletal_Muscle_Myoblasts	PrimaryCulture
E121	ENCODE2012	MUSCLE	ENCODE2012	MUS.HSMMT	HSMM cell derived Skeletal Muscle Myotubes Cells	HSMMtube_Skeletal_Muscle_Myotubes_Derived_from_HSMM	PrimaryCulture
E122	ENCODE2012	VASCULAR	ENCODE2012	VAS.HUVEC	HUVEC Umbilical Vein Endothelial Primary Cells	HUVEC_Umbilical_Vein_Endothelial_Cells	PrimaryCulture
E123	ENCODE2012	BLOOD	ENCODE2012	BLD.K562.CNCR	K562 Leukemia Cells	K562_Leukemia	PrimaryCulture
E124	ENCODE2012	BLOOD	ENCODE2012	BLD.CD14.MONO	Monocytes-CD14+ RO01746 Primary Cells	Monocytes-CD14+_RO01746	PrimaryCell
E125	ENCODE2012	BRAIN	ENCODE2012	BRN.NHA	NH-A Astrocytes Primary Cells	NH-A_Astrocytes	PrimaryCulture
E126	ENCODE2012	SKIN	ENCODE2012	SKIN.NHDFAD	NHDF-Ad Adult Dermal Fibroblast Primary Cells	NHDF-Ad_Adult_Dermal_Fibroblasts	PrimaryCulture
E127	ENCODE2012	SKIN	ENCODE2012	SKIN.NHEK	NHEK-Epidermal Keratinocyte Primary Cells	NHEK-Epidermal_Keratinocytes	PrimaryCulture
E128	ENCODE2012	LUNG	ENCODE2012	LNG.NHLF	NHLF Lung Fibroblast Primary Cells	NHLF_Lung_Fibroblasts	PrimaryCulture
E129	ENCODE2012	BONE	ENCODE2012	BONE.OSTEO	Osteoblast Primary Cells	Osteoblasts	PrimaryCulture
E027	Epithelial	BREAST	Epithelial	BRST.MYO	Breast Myoepithelial Primary Cells	Breast_Myoepithelial_Cells	PrimaryCell
E028	Epithelial	BREAST	Epithelial	BRST.HMEC.35	Breast variant Human Mammary Epithelial Cells (vHMEC)	Breast_vHMEC	PrimaryCulture
E055	Epithelial	SKIN	Epithelial	SKIN.PEN.FRSK.FIB.01	Foreskin Fibroblast Primary Cells skin01	Penis_Foreskin_Fibroblast_Primary_Cells_skin01	PrimaryCulture
E056	Epithelial	SKIN	Epithelial	SKIN.PEN.FRSK.FIB.02	Foreskin Fibroblast Primary Cells skin02	Penis_Foreskin_Fibroblast_Primary_Cells_skin02	PrimaryCulture
E057	Epithelial	SKIN	Epithelial	SKIN.PEN.FRSK.KER.02	Foreskin Keratinocyte Primary Cells skin02	Penis_Foreskin_Keratinocyte_Primary_Cells_skin02	PrimaryCulture
E058	Epithelial	SKIN	Epithelial	SKIN.PEN.FRSK.KER.03	Foreskin Keratinocyte Primary Cells skin03	Penis_Foreskin_Keratinocyte_Primary_Cells_skin03	PrimaryCulture
E059	Epithelial	SKIN	Epithelial	SKIN.PEN.FRSK.MEL.01	Foreskin Melanocyte Primary Cells skin01	Penis_Foreskin_Melanocyte_Primary_Cells_skin01	PrimaryCulture
E061	Epithelial	SKIN	Epithelial	SKIN.PEN.FRSK.MEL.03	Foreskin Melanocyte Primary Cells skin03	Penis_Foreskin_Melanocyte_Primary_Cells_skin03	PrimaryCulture
E004	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.H1.BMP4.MESO	H1 BMP4 Derived Mesendoderm Cultured Cells	H1_BMP4_Derived_Mesendoderm_Cultured_Cells	ESCDerived
E005	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.H1.BMP4.TROP	H1 BMP4 Derived Trophoblast Cultured Cells	H1_BMP4_Derived_Trophoblast_Cultured_Cells	ESCDerived
E006	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.H1.MSC	H1 Derived Mesenchymal Stem Cells	H1_Derived_Mesenchymal_Stem_Cells	ESCDerived
E007	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.H1.NEUR.PROG	H1 Derived Neuronal Progenitor Cultured Cells	H1_Derived_Neuronal_Progenitor_Cultured_Cells	ESCDerived
E009	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.H9.NEUR.PROG	H9 Derived Neuronal Progenitor Cultured Cells	H9_Derived_Neuronal_Progenitor_Cultured_Cells	ESCDerived
E010	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.H9.NEUR	H9 Derived Neuron Cultured Cells	H9_Derived_Neuron_Cultured_Cells	ESCDerived
E011	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.CD184.ENDO	hESC Derived CD184+ Endoderm Cultured Cells	hESC_Derived_CD184+_Endoderm_Cultured_Cells	ESCDerived
E012	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.CD56.ECTO	hESC Derived CD56+ Ectoderm Cultured Cells	hESC_Derived_CD56+_Ectoderm_Cultured_Cells	ESCDerived
E013	ES-deriv	ESC_DERIVED	ES-deriv	ESDR.CD56.MESO	hESC Derived CD56+ Mesoderm Cultured Cells	hESC_Derived_CD56+_Mesoderm_Cultured_Cells	ESCDerived
E001	ESC	ESC	ESC	ESC.I3	ES-I3 Cells	ES-I3_Cell_Line	PrimaryCulture
E002	ESC	ESC	ESC	ESC.WA7	ES-WA7 Cells	ES-WA7_Cell_Line	PrimaryCulture
E003	ESC	ESC	ESC	ESC.H1	H1 Cells	H1_Cell_Line	PrimaryCulture
E008	ESC	ESC	ESC	ESC.H9	H9 Cells	H9_Cell_Line	PrimaryCulture
E014	ESC	ESC	ESC	ESC.HUES48	HUES48 Cells	HUES48_Cell_Line	PrimaryCulture
E015	ESC	ESC	ESC	ESC.HUES6	HUES6 Cells	HUES6_Cell_Line	PrimaryCulture
E016	ESC	ESC	ESC	ESC.HUES64	HUES64 Cells	HUES64_Cell_Line	PrimaryCulture
E024	ESC	ESC	ESC	ESC.4STAR	ES-UCSF4  Cells	4star	PrimaryCulture
E065	Heart	VASCULAR	Heart	VAS.AOR	Aorta	Aorta	PrimaryTissue
E083	Heart	HEART	Heart	HRT.FET	Fetal Heart	Fetal_Heart	PrimaryTissue
E095	Heart	HEART	Heart	HRT.VENT.L	Left Ventricle	Left_Ventricle	PrimaryTissue
E104	Heart	HEART	Heart	HRT.ATR.R	Right Atrium	Right_Atrium	PrimaryTissue
E105	Heart	HEART	Heart	HRT.VNT.R	Right Ventricle	Right_Ventricle	PrimaryTissue
E029	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.CD14.PC	Primary monocytes from peripheral blood	CD14_Primary_Cells	PrimaryCell
E030	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.CD15.PC	Primary neutrophils from peripheral blood	CD15_Primary_Cells	PrimaryCell
E031	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.CD19.CPC	Primary B cells from cord blood	CD19_Primary_Cells_Cord_BI	PrimaryCell
E032	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.CD19.PPC	Primary B cells from peripheral blood	CD19_Primary_Cells_Peripheral_UW	PrimaryCell
E035	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.CD34.PC	Primary hematopoietic stem cells	CD34_Primary_Cells	PrimaryCell
E036	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.CD34.CC	Primary hematopoietic stem cells short term culture	CD34_Cultured_Cells	PrimaryCell
E046	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.CD56.PC	Primary Natural Killer cells from peripheral blood	CD56_Primary_Cells	PrimaryCell
E050	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.MOB.CD34.PC.F	Primary hematopoietic stem cells G-CSF-mobilized Female	Mobilized_CD34_Primary_Cells_Female	PrimaryCell
E051	HSC & B-cell	BLOOD	HSC_and_Bcell	BLD.MOB.CD34.PC.M	Primary hematopoietic stem cells G-CSF-mobilized Male	Mobilized_CD34_Primary_Cells_Male	PrimaryCell
E018	iPSC	IPSC	iPSC	IPSC.15b	iPS-15b Cells	iPS-15b_Cell_Line	PrimaryCulture
E019	iPSC	IPSC	iPSC	IPSC.18	iPS-18 Cells	iPS-18_Cell_Line	PrimaryCulture
E020	iPSC	IPSC	iPSC	IPSC.20B	iPS-20b Cells	iPS-20b_Cell_Line	PrimaryCulture
E021	iPSC	IPSC	iPSC	IPSC.DF.6.9	iPS DF 6.9 Cells	iPS_DF_6.9_Cell_Line	PrimaryCulture
E022	iPSC	IPSC	iPSC	IPSC.DF.19.11	iPS DF 19.11 Cells	iPS_DF_19.11_Cell_Line	PrimaryCulture
E023	Mesench	FAT	Mesench	FAT.MSC.DR.ADIP	Mesenchymal Stem Cell Derived Adipocyte Cultured Cells	Mesenchymal_Stem_Cell_Derived_Adipocyte_Cultured_Cells	PrimaryCulture
E025	Mesench	FAT	Mesench	FAT.ADIP.DR.MSC	Adipose Derived Mesenchymal Stem Cell Cultured Cells	Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells	PrimaryCulture
E026	Mesench	STROMAL_CONNECTIVE	Mesench	STRM.MRW.MSC	Bone Marrow Derived Cultured Mesenchymal Stem Cells	Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells	PrimaryCulture
E049	Mesench	STROMAL_CONNECTIVE	Mesench	STRM.CHON.MRW.DR.MSC	Mesenchymal Stem Cell Derived Chondrocyte Cultured Cells	Chondrocytes_from_Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells	PrimaryCulture
E052	Myosat	MUSCLE	Muscle	MUS.SAT	Muscle Satellite Cultured Cells	Muscle_Satellite_Cultured_Cells	PrimaryCulture
E089	Muscle	MUSCLE	Muscle	MUS.TRNK.FET	Fetal Muscle Trunk	Fetal_Muscle_Trunk	PrimaryTissue
E090	Muscle	MUSCLE_LEG	Muscle	MUS.LEG.FET	Fetal Muscle Leg	Fetal_Muscle_Leg	PrimaryTissue
E100	Muscle	MUSCLE	Muscle	MUS.PSOAS	Psoas Muscle	Psoas_Muscle	PrimaryTissue
E107	Muscle	MUSCLE	Muscle	MUS.SKLT.M	Skeletal Muscle Male	Skeletal_Muscle_Male	PrimaryTissue
E108	Muscle	MUSCLE	Muscle	MUS.SKLT.F	Skeletal Muscle Female	Skeletal_Muscle_Female	PrimaryTissue
E017	IMR90	LUNG	Other	LNG.IMR90	IMR90 fetal lung fibroblasts Cell Line	IMR90_Cell_Line	CellLine
E063	Adipose	FAT	Other	FAT.ADIP.NUC	Adipose Nuclei	Adipose_Nuclei	PrimaryTissue
E066	Other	LIVER	Other	LIV.ADLT	Liver	Adult_Liver	PrimaryTissue
E076	Sm. Muscle	GI_COLON	Other	GI.CLN.SM.MUS	Colon Smooth Muscle	Colon_Smooth_Muscle	PrimaryTissue
E078	Sm. Muscle	GI_DUODENUM	Other	GI.DUO.SM.MUS	Duodenum Smooth Muscle	Duodenum_Smooth_Muscle	PrimaryTissue
E080	Other	ADRENAL	Other	ADRL.GLND.FET	Fetal Adrenal Gland	Fetal_Adrenal_Gland	PrimaryTissue
E086	Other	KIDNEY	Other	KID.FET	Fetal Kidney	Fetal_Kidney	PrimaryTissue
E087	Other	PANCREAS	Other	PANC.ISLT	Pancreatic Islets	Pancreatic_Islets	PrimaryTissue
E088	Other	LUNG	Other	LNG.FET	Fetal Lung	Fetal_Lung	PrimaryTissue
E091	Other	PLACENTA	Other	PLCNT.FET	Placenta	Fetal_Placenta	PrimaryTissue
E093	Thymus	THYMUS	Other	THYM.FET	Fetal Thymus	Fetal_Thymus	PrimaryTissue
E096	Other	LUNG	Other	LNG	Lung	Lung	PrimaryTissue
E097	Other	OVARY	Other	OVRY	Ovary	Ovary	PrimaryTissue
E098	Other	PANCREAS	Other	PANC	Pancreas	Pancreas	PrimaryTissue
E099	Other	PLACENTA	Other	PLCNT.AMN	Placenta Amnion	Placenta_Amnion	PrimaryTissue
E103	Sm. Muscle	GI_RECTUM	Other	GI.RECT.SM.MUS	Rectal Smooth Muscle	Rectal_Smooth_Muscle	PrimaryTissue
E111	Sm. Muscle	GI_STOMACH	Other	GI.STMC.MUS	Stomach Smooth Muscle	Stomach_Smooth_Muscle	PrimaryTissue
E112	Thymus	THYMUS	Other	THYM	Thymus	Thymus	PrimaryTissue
E113	Other	SPLEEN	Other	SPLN	Spleen	Spleen	PrimaryTissue
"""