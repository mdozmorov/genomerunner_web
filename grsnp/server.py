import cherrypy, os, cgi, tempfile, sys, itertools
from mako.template import Template
from mako.lookup import TemplateLookup
from mako.exceptions import RichTraceback as MakoTraceback
from contextlib import closing
import sqlite3
import re
from operator import attrgetter
from multiprocessing import Process
import cPickle
from path import PathNode
import path as grsnp_path
from operator import itemgetter
from path import base_name as basename
import logging
from logging import FileHandler,StreamHandler
import json
import pdb
import grsnp.worker_hypergeom4
from time import gmtime, strftime
import simplejson
import string
import random
import traceback
import dbcreator_ucsc as uscsreader
import argparse
import shutil
import subprocess
import celeryconfiguration
from celery import Celery
import grp
import pwd

os.environ['GR_COMPATIBILITY_MODE'] = 'y'

root_dir = os.path.dirname(os.path.realpath(__file__))
lookup = TemplateLookup(directories=[os.path.join(root_dir,"frontend/templates")])


sett = {}
DEBUG_MODE = True
logger = logging.getLogger('genomerunner.server')





# Each function in this class is a web page 
class WebUI(object):
	def __init__(self):		
		# go through each database directory and create custom_data if it does not exist.
		for db_ver,db_dir in sett["data_dir"].items():
			# create all directories in the custom_data dir if they do not already exist
			for org in self.get_org(db_ver):
				custom_dir = os.path.join(os.path.split(db_dir)[0],"custom_data")
				logger.info("Processing genomic features for {}".format(org))
				if not os.path.exists(custom_dir): os.mkdir(custom_dir)
				cust_sub_dir = ["backgrounds","gfs","fois","rsid_conversion"]
				for c in cust_sub_dir:
					tmp = os.path.join(custom_dir,c)
					if not os.path.exists(tmp): os.mkdir(tmp)
					c_dir = os.path.join(custom_dir,c,org)
					if not os.path.exists(c_dir): os.mkdir(c_dir)
				# Read the genomic feature files and generate html files
				paths = PathNode()
				paths.name = "Root"
				paths.organisms = self.get_org(db_ver) 
				paths.traverse(os.path.join(db_dir,org))
				grsnp_path.write_treeview_json(os.path.join(db_dir,org))
		self._index_html = {}

	@cherrypy.expose
	def index(self,organism=None,db_version=None):
		global results_dir, uploads_dir, sett
		if not organism: organism = sett["default_organism"]

		if DEBUG_MODE or not organism in self._index_html:		
			tmpl = lookup.get_template("master.mako")
			paths = PathNode()
			paths.name = "Root"
			[html_dbversion, db_version] =  grsnp_path.get_database_versions_html(sett["data_dir"],db_version)
			paths.organisms = self.get_org(db_version)
			# Check if the organism actually exists in the current database.
			# If it is not, select the default organism
			if not organism in paths.organisms:
				organism = sett["default_organism"]
			custom_dir = os.path.join(os.path.split(sett["data_dir"][db_version])[0],"custom_data")
			# Use mako to render index.html
			body = lookup.get_template("index.mako").render(paths=paths,default_background=paths.get_backgrounds_combo(organism,custom_dir),
									custom_gfs=paths.get_custom_gfs(organism,custom_dir),demo_snps=paths.get_custom_fois(organism,custom_dir),
									data_dir=os.path.join(sett["data_dir"][db_version],organism),default_organism=organism,
									database_versions=html_dbversion,pct_scores=paths.get_scores(os.path.split(sett["data_dir"][db_version])[0]))
			script = lookup.get_template("index.js").render(default_organism=organism)
			self._index_html[organism] = tmpl.render(body=body,script=script)

		return self._index_html[organism]

	def get_org(self,db_version):
		organisms = []
		files = os.listdir(sett["data_dir"][db_version])
		for f in files:
			if f.find(".") == -1:
				organisms.append(f)
		return organisms	

	@cherrypy.expose
	def query(self, bed_file=None,bed_data=None, background_file=None,background_data=None, 
				genomicfeature_file=None, niter=10, name="", strand="",run_annotation=False, default_background = "",db_version=None,jstree_gfs="",**kwargs):
		global results_dir, uploads_dir, sett
		# Assign a random id
		id = ''.join(random.choice(string.lowercase+string.digits) for _ in range(32))
		while (os.path.exists(os.path.join(uploads_dir,id))):
			id = ''.join(random.choice(string.lowercase+string.digits) for _ in range(32))
		res_dir = os.path.join(results_dir,str(id))
		upload_dir = os.path.join(uploads_dir,str(id))
		os.mkdir(upload_dir)
		os.mkdir(os.path.join(upload_dir,"fois"))
		os.mkdir(os.path.join(upload_dir,"gfs"))
		res_dir = os.path.join(results_dir,str(id))
		os.mkdir(res_dir)
		#if sett['group'] != "":
		#	gid = grp.getgrnam(sett['group']).gr_gid
		#	#os.chown(res_dir, -1,gid)
		#os.chmod(res_dir,0x777)
		
		fois = os.path.join(upload_dir,".fois") # contains a list of the paths to fois to run through the analysis
		gfs = os.path.join(upload_dir,".gfs") # contains a list of the paths to the gfs to run the fois against
		list_gfs = []
		data_dir = os.path.split(sett["data_dir"][db_version])[0]
		cherrypy.response.timeout = 3600
		try:
			jobname = kwargs["jobname"]
		except Exception, e:
			jobname = ""
			logger.error("id={}".format(id) + str(e))	

		# load the FOI data
		bed_filename,data = "",""
		demo_fois_dir =  kwargs["demo_fois"] if "demo_fois" in kwargs else ""  # If the user selects a demo set of FOIs to run, this will contain the directory
		if demo_fois_dir == "":
			try:
				with open(fois,"wb") as out_fois:
					# bed files uploaded
					if not isinstance(bed_file,(list)): bed_file = [bed_file] # makes a list if only one file uploaded
					if bed_file[0] and bed_file[0].filename != "":						
						for b in bed_file:
							bed_filename = b.filename
							f = os.path.join(upload_dir, "fois",bed_filename)
							extension = bed_filename.split(".")[-1]
							if not os.path.exists(f):
								with open(f, "wb") as out:
									if b != None and b.filename != "":
										logger.info("Received uploaded FOI file (name={}, id={})".format(bed_filename, id))
										while True:
											data = b.file.read(8192)
											# TODO find empty lines
											#data = os.linesep.join([s for s in data.splitlines() if s ])

											# strips out new lines not compatible with bed tools
											if extension not in ["gz","bb"]: data.replace("\r","")
											if not data:
												break
											out.write(data)			
										out_fois.write(f+"\n")
								# use dos2unix to remove \r from end of lines
								out = subprocess.Popen(['dos2unix {}'.format(f)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
								out.wait()								
								#script = "sort -k1,1 -k2,2n -k3,3n " + path +" | bgzip -c > " + outpath + ".gz.temp"
								#out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE)
								#out.wait()
							else:
								logger.error("id={} Upload file already exists at {}".format(id,f))
								print "id={} Upload file already exists at {}".format(id,f)
					# custom data entered	
					elif bed_data:
						f = os.path.join(upload_dir,"fois", "custom.bed")
						with open(f, "wb") as out:
							bed_filename = "custom.bed"
							logger.info('Received raw text  FOI data (id={})'.format(id))
							data = bed_data
							data = os.linesep.join([s for s in data.splitlines() if s])
							out.write(data)		
						out_fois.write(f+"\n")
						# use dos2unix to remove \r from end of lines
						out = subprocess.Popen(['dos2unix {}'.format(f)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
						out.wait()	
					else:
						return "ERROR: Please upload a feature of interest file."

			except Exception, e:
				logger.error("id={}".format(id) + str(e))
				return "ERROR: upload a file please"
		elif demo_fois_dir != "":
			# gather the FOI files in the demo directory
			ls_foi = [os.path.join(demo_fois_dir,f) for f in os.listdir(demo_fois_dir) if os.path.isfile(os.path.join(demo_fois_dir,f))]
			with open(fois,"wb") as writer:
				for f in ls_foi:
					writer.write(f+"\n")					
		else:
			return "Feature of Interest files not detected.  Please upload or choose Feature of Interest to run."

		

		# uploads custom genomic features
		try:
			with open(gfs,"a") as out_gfs:
				# bed files uploaded
				if genomicfeature_file:
					if not isinstance(genomicfeature_file,(list)): genomicfeature_file = [genomicfeature_file] # makes a list if only one file uploaded
					for b in genomicfeature_file:
						gfbed_filename = b.filename
						f = os.path.join(upload_dir, "gfs", gfbed_filename)
						extension = gfbed_filename.split(".")[-1]
						if not os.path.exists(f):
							with open(f, "wb") as out:
								if b != None and b.filename != "":
									logger.info("Received uploaded GF file (name={}, id={})".format(gfbed_filename, id))
									while True:
										data = b.file.read(8192)
										# TODO find empty lines
										#data = os.linesep.join([s for s in data.splitlines() if s ])

										# strips out new lines not compatible with bed tools
										if extension not in ["gz","bb"]: data.replace("\r","")
										if not data:
											break
										out.write(data)	
									out_gfs.write(f+"\n")
									# use dos2unix to remove \r from end of lines
							out = subprocess.Popen(['dos2unix {}'.format(f)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
							out.wait()
							list_gfs.append(base_name(f))		
						else:
							logger.error("id={} Uploaded GF file already exists at {}".format(id,f))
		except Exception, e:
			logger.error("id={}".format(id) + str(e))
			return "ERROR: Unable to process custom Genome annotation feature"

		
		organism,run,run_random = "",[],False

		# add genomic feature tracks loaded via the JStree control
		with open(gfs,"a") as out_gfs:
			for k in jstree_gfs.split(','):
				if k.startswith('file:'):
					g = k.split(":")[-1]
					if (base_name(g) not in list_gfs): 
						# check if score and/or strand filtered GF data exists
						g = verify_score_strand(g,kwargs['pct_score'],strand,data_dir)
						out_gfs.write(g+"\n")
					list_gfs.append(base_name(g))
		for k,v in kwargs.items():
			# organism to use
			if "organism:" in v:
				organism = v.split(":")[-1]
			# which tests to run
			if "run:" in k and v=="on":
				run.append(k.split(":")[-1])
			# append custom list to be run
			if "grouprun:" in k and v == "on":
				gp_gfs_dir = k.split(":")[-1]
				ls_foi = [os.path.join(gp_gfs_dir,f) for f in os.listdir(gp_gfs_dir) if os.path.isfile(os.path.join(gp_gfs_dir,f))]
				with open(gfs,"a") as writer:
					for f in ls_foi:
						if (base_name(f) not in list_gfs): writer.write(f.rstrip(".tbi")+"\n")
						list_gfs.append(base_name(f))	
			if k.startswith("run_random") and v == "on": run_random = True
			if k.startswith("run_annot") and v == "on": run_annotation = True


		# Create annotation folder, used by the server to check if annotation is going to be run
		if run_annotation:
			annot_outdir = os.path.join(res_dir,"annotations")
			if not os.path.exists(annot_outdir): os.mkdir(annot_outdir)
		

		# load the background data if uploaded
		background_name = ""
		try:
			if background_file != None and background_file.filename != "":
				print(background_file)
				b = os.path.join(upload_dir,background_file.filename)
				logger.info('Received uploaded background file (id={})'.format(id))
				background_name = background_file.filename
				extension = background_name.split(".")[-1]
				with open(b, "wb") as out:
					while True:
						data = background_file.file.read(8192)
						# TODO find empty lines
						#data = os.linesep.join([s for s in data.splitlines() if not s.isspace() ])

						# strips out new lines not compatible with bed tools
						if extension not in ["gz","bb"]: data = data.replace("\r","")
						if not data:
							break
						out.write(data)
				# use dos2unix to remove \r from end of lines
				out = subprocess.Popen(['dos2unix {}'.format(b)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
				out.wait()	
			elif background_data != None and background_data != "":
				b = os.path.join(upload_dir,"custom_background.bed")
				background_name = "custom.bed"
				with open(b, "wb") as out:
					logger.info('Received raw text background data (id={})'.format(id))
					data = os.linesep.join([s for s in background_data.split("\n") if s != ""])
					out.write(data)
				# use dos2unix to remove \r from end of lines
				out = subprocess.Popen(['dos2unix {}'.format(b)],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
				out.wait()
			else:				
				background_name = default_background.split("/")[-1].split(".")[0] 
				b = default_background
		except Exception, e:
			logger.error("id={}".format(id) + str(e))
			return "ERROR: unable to upload background"

		# make paths relative, needed for remote celery workers to function correctly
		for f in [fois,gfs]:
			list_foi = open(f).read().replace(uploads_dir,'/uploads').replace(results_dir,'/results').replace(data_dir,"")
			with open(f,'wb') as writer:
				writer.write(list_foi)
		b = b.replace(data_dir,'').replace(os.path.split(uploads_dir)[0],"").lstrip("/")

		# write the enrichment settings.
		path = os.path.join(res_dir, ".settings")
		set_info = {"Jobname:": str(id),
					"Time:": strftime("%Y-%m-%d %H:%M:%S", gmtime()),
					"Background:": background_name,
					"Organism:": organism,
					"Database version:":db_version,
					"% Score threshold:": str(kwargs['pct_score'])+"%",
					"Strand:": strand}

		with open(path, 'wb') as sett_files:
			for k,v in set_info.iteritems():
				sett_files.write(k+"\t"+v+"\n")

		gfs_count = 0
		if open(gfs).read() == "": 
			return "ERROR: No Genomic Features selected/uploaded."
		else:
			gfs_count = len([x for x in open(gfs).read().split("\n") if x != ""])

		# copy over .foi to results folder, used for the enrichment results
		shutil.copy(fois,os.path.join(res_dir,".fois"))

		# This starts the enrichment analysis in another OS process.
		# We know it is done when a file appears in the "results" directory
		#p = Process(target=grquery.run_hypergeom,
		#		args=(fois,gfs,b,res_dir,id,True,os.path.join(sett["data_dir"],organism,"bkg_overlaps.gr"),sett["data_dir"],run_annotation,run_random))				
		#p.start()

		# run using celery queues.  Uncomment to use.  Also uncomment in celeryconfiguration.
		# TODO find out why all jobs are getting sent to one worker.
		#if gfs_count > 10:
		#	print "LONG RUN STARTED"
		#	grsnp.worker_hypergeom4.run_hypergeom.apply_async(args=[fois,gfs,b,res_dir,id,True,os.path.join(sett["data_dir"],organism,"bkg_overlaps.gr"),sett["data_dir"],run_annotation,run_random],
		#														  queue='long_runs')
		#else:
		#	print "SHORT RUN STARTED"
		#	grsnp.worker_hypergeom4.run_hypergeom.apply_async(args=[fois,gfs,b,res_dir,id,True,os.path.join(sett["data_dir"],organism,"bkg_overlaps.gr"),sett["data_dir"],run_annotation,run_random],
		#														  queue='short_runs')
		
		try:
			grsnp.worker_hypergeom4.run_hypergeom.delay(fois,gfs,b,id,True,os.path.join(sett["data_dir"][db_version],organism,"bkg_overlaps.gr"),run_annotation,run_random,pct_score=kwargs['pct_score'],organism=organism,id=id,db_version=db_version)
		except Exception, e:
			print "WORKER ERROR"
		raise cherrypy.HTTPRedirect("result?id=%s" % id)

	@cherrypy.expose
	def result(self, id):
		global results_dir, uploads_dir, sett
		path = os.path.join(results_dir, id)
		params = {}
		params["run_id"] = id
		params["detailed"] = "Results not yet available"
		params["matrix"] = "Results not yet available"		
		tmpl = lookup.get_template("master.mako")

		# Loads the progress file if it exists
		p = {"status":"","curprog":0,"progmax":0}
		progress_path = os.path.join(path,".prog")
		if os.path.exists(progress_path):
			with open(progress_path) as f:
				p = json.loads(f.read() )
		params["log"] = "###Run Settings###\n"
		sett_path = os.path.join(path,".settings")
		organism = ""
		if os.path.exists(sett_path):
			with open(sett_path) as f:
				tmp = f.read()	
				params["log"] = params["log"]+ tmp
				organism = [x.split("\t")[1] for x in tmp.split("\n") if x.split("\t")[0] == "Organism:"][0]
		params["organism"] = organism
		params["log"] = params["log"] + "\n###Run Log###\n"
		debug_path = os.path.join(path,".log")
		if os.path.exists(debug_path):
			with open(debug_path) as f:
				params["log"] = params["log"] + f.read()

		# loads results from results file		
		detailed_path = os.path.join(path,"detailed.gr")
		if os.path.exists(detailed_path):
			with open(detailed_path) as f:
				params["detailed"] = f.read()
		
		foi_names_path = os.path.join(os.path.join(results_dir, id),".fois")
		if os.path.exists(foi_names_path):
			with open(foi_names_path) as f:
				params["fois"] = [basename(x).split(".")[0] for x in f.read().split("\n") if x != ""]
		else:
			params["fois"] = ""

		params["zipfile"] = os.path.join("results",id,"GR_{}.tar.gz").format(id)
		params["run_annotation"] = True if os.path.exists(os.path.join(results_dir,id,"annotations")) else  False
		params.update(p)
		try:
			rend_template = tmpl.render(body=lookup.get_template("results.mako").render(**params),script= lookup.get_template("results.js").render(**params))
			print "LOADED TEMPLATE"
		except Exception, e:
			traceback = MakoTraceback()
			str_error = ""
			for (filename, lineno, function, line) in traceback.traceback:
				str_error +=  "File %s, line %s, in %s" % (os.path.split(filename)[-1], lineno, function)
				str_error += "\n"
				str_error += line + "\n"
				str_error += "%s: %s" % (str(traceback.error.__class__.__name__), traceback.error)
			print str_error
			rend_template = str_error
		return rend_template

	@cherrypy.expose
	def results_shiny(self, id):
		global results_dir, uploads_dir, sett
		path = os.path.join(results_dir, id)	
		params = {}	
		params['run_id'] = id
		try:
			tmp = lookup.get_template("master.mako")
			script = lookup.get_template("results_shiny.js").render(run_id=id)
			rend_template = lookup.get_template("results_shiny.mako").render(run_id=id,script=script)

			#rend_template = lookup.get_template("master.mako").render(body = lookup.get_template("results_shiny.mako").render(), 
			#	script = lookup.get_template("results_shiny.js").render(run_id=id))
			#rend_template = lookup.get_template("master.mako").render(body=lookup.get_template("results_shiny.mako").render(run_id=id),script=lookup.get_template("results_shiny.js"))
			#rend_template = lookup.get_template("results_shiny.mako").render(script= lookup.get_template("results_shiny.js").render(**params),**params)
			

		except Exception, e:
			traceback = MakoTraceback()
			str_error = ""
			for (filename, lineno, function, line) in traceback.traceback:
				str_error +=  "File %s, line %s, in %s" % (os.path.split(filename)[-1], lineno, function)
				str_error += "\n"
				str_error += line + "\n"
				str_error += "%s: %s" % (str(traceback.error.__class__.__name__), traceback.error)
			print str_error
			rend_template = str_error
		return rend_template

	@cherrypy.expose
	def gf_descriptions(self,db_version,organism):
		# Use mako to render index.html
		tmpl = lookup.get_template("master.mako")
		body = lookup.get_template("gf_descriptions.mako").render()
		script = lookup.get_template("gf_descriptions.js").render(db_version=db_version,organism=organism)
		return tmpl.render(body=body,script=script)


	@cherrypy.expose
	def get_heatmaps(self, run_id, organism):
		"""	Returns clustered and PCC matrix if they exist.
		'organism': is used to load detailed labels for the GFs.
		"""
		global results_dir, uploads_dir, sett
		cherrypy.response.headers['Content-Type'] = 'application/json'
		trackdb = uscsreader.load_tabledata_dumpfiles(os.path.join(sett["data_dir"],organism,"trackDb"))
		results = {}
		path = os.path.join(results_dir, run_id)
		matrix_path = os.path.join(path,"matrix.txt")
		# Load clustered matrix
		matrix_clust_path  = os.path.join(path,"clustered.txt")
		if os.path.exists(matrix_path):
			with open(matrix_path) as f:
				results["matrix"] = f.read().replace("\"","")
		if os.path.exists(matrix_clust_path):
			with open(matrix_clust_path) as f:
				d = f.read()
				d = d.replace("\"","")
				if d[:6] != "ERROR:":
					results["matrix_data"] = d.replace("\"","")
					# d3 requires "gene_name" to be inserted into the first column
					tmp_data =  results["matrix_data"].split("\n")
					tmp = tmp_data[:] # this copy is passed onto the results page
					# insert the description column if it does not exist
					tmp_matrix = [x.split("\t") for x in tmp]
					if tmp_matrix[0][-1] != "Genomic Feature Description":
						tmp_matrix[0] += ["Genomic Feature Description"]
						for i in range(1,len(tmp)):
							description = [x["longLabel"] for x in trackdb if x["tableName"] == tmp_matrix[i][0]]
							if len(description) is not 0: description = description[0]
							else: description = ""
							tmp_matrix[i] += [description]					

					results["matrix_data"] = "\n".join(["\t".join(["gene_name",tmp[0]])]+tmp[1:])  
					results["matrix_data"] = results["matrix_data"]
					results["matrix_data_gf_description"] = "\t".join([x[-1] for x in tmp_matrix[1:]])
				else: results["matrix_data"] = d[6:]
		else: 
			results["matrix_data"] = "Heatmap will be available after the analysis is complete."
			results["matrix_data_gf_description"] = ""

		# Pearson's matrix results
		matrix_cor_path = os.path.join(path,"pcc_matrix.txt")
		if os.path.exists(matrix_cor_path):
			with open(matrix_cor_path) as f:
				d = f.read()
				d = d.replace("\"","")
				if d[:6] != "ERROR:":
					results["matrix_cor_data"] =  d.replace("\"","")
					results["matrix_cor"] = d.replace("\"","")
					# d3 requires "gene_name" to be inserted into the first column
					tmp =  results["matrix_cor"].split("\n")
					results["matrix_cor"] = "\n".join(["\t".join(["gene_name",tmp[0]])]+tmp[1:])  
					results["matrix_cor"] = results["matrix_cor"]
				else: results["matrix_cor"],results["matrix_cor_data"] = d[6:],""				

		else:
			results["matrix_cor_data"] = ""
			results["matrix_cor"] = "PCC heatmap will be available after the clustered matrix is created."
		pvalue_path = os.path.join(path,"pcc_matrix_pvalue.txt")
		if os.path.exists(pvalue_path):
			with open(pvalue_path) as f:
				d = f.read()
				results["matrix_cor_pvalues"] = d.replace("\"","")
				# d3 requires "gene_name" to be inserted into the first column
				tmp =  results["matrix_cor_pvalues"].split("\n")
				results["matrix_cor_pvalues"] = "\n".join(["\t".join(["gene_name",tmp[0]])]+tmp[1:])   
				results["matrix_cor_pvalues"] = results["matrix_cor_pvalues"]			
		else: 
			results["matrix_cor_pvalues"] = ""
		return simplejson.dumps(results)

	@cherrypy.expose
	def get_annotation(self,run_id,foi_name):
		annotation_path = os.path.join(results_dir,run_id,"annotations",foi_name + ".txt")
		results = []
		if os.path.exists(annotation_path):
			with open(annotation_path) as f:
				# skip the comment lines
				cols = f.readline().rstrip()
				while cols[0] == "#":
					cols = f.readline().rstrip()
				cols = cols.split("\t")	
				results.append(cols)			
				for foi in f:
					if foi.strip() != "":
						results.append(foi.rstrip().split("\t"))
		return simplejson.dumps(results)

	@cherrypy.expose
	def get_enrichment(self,run_id,foi_name):
		global results_dir, uploads_dir, sett
		enrichment_path = os.path.join(results_dir,run_id,"enrichment",foi_name + ".txt")
		results = []
		if os.path.exists(enrichment_path):
			with open(enrichment_path) as f:
				# skip the comment lines
				cols = f.readline().rstrip()
				while cols[0] == "#":
					cols = f.readline().rstrip()
				cols = cols.split("\t")	
				results.append(cols)			
				for foi in f:
					if foi.strip() != "":
						results.append(foi.rstrip().split("\t"))
		return simplejson.dumps(results)

	@cherrypy.expose
	def get_gf_descriptions(self,db_version,organism):
		global results_dir, uploads_dir, sett
		descriptions_path = os.path.join(sett["data_dir"][db_version],organism,"gf_descriptions.txt")		
		results = []
		if os.path.exists(descriptions_path):
			with open(descriptions_path) as f:
				# skip the comment lines
				cols = f.readline().rstrip()
				while cols[0] == "#":
					cols = f.readline().rstrip()
				cols = cols.split("\t")	
				results.append(cols)			
				for foi in f:
					if foi.strip() != "":
						results.append(foi.rstrip().split("\t"))
		return simplejson.dumps(results)


	@cherrypy.expose
	def get_cluster(self,run_id):
		global results_dir, uploads_dir, sett
		mat_path = os.path.join(results_dir,run_id,"clustered.json")
		if os.path.exists(mat_path):
			with open(mat_path) as f:
				return f.read()
		return simplejson.dumps([{"matrix": "\"INFO: No Results for the matrix found. Were enough features analyzed?"}])
	
	@cherrypy.expose
	def get_pcc(self,run_id):
		global results_dir, uploads_dir, sett
		mat_path = os.path.join(results_dir,run_id,"pcc_matrix.json")
		if os.path.exists(mat_path):
			with open(mat_path) as f:
				return f.read()
		return simplejson.dumps([{"matrix": "\"INFO: No Results for the matrix found. Were enough features analyzed?"}])			

	@cherrypy.expose
	def meta(self, tbl,organism,db_version):
		"""Returns the html description from the trackDb file for the specified organism.
		"""
		global results_dir, uploads_dir, sett
		try:
			trackdb = uscsreader.load_tabledata_dumpfiles(os.path.join(sett["data_dir"][db_version],organism,"trackDb"))
			html = trackdb[map(itemgetter('tableName'),trackdb).index(tbl)]['html']
		except Exception, e:
			return "<h3>(No data found for {}.)</h3>".format(tbl)
		if html=='':
			return "<h3>(No data found for {}.)</h3>".format(tbl)
		else:
			return html

	@cherrypy.expose
	def get_detailed(self,run_id):
		""" loads results from detailed results file
		"""
		global results_dir, uploads_dir, sett
		detailed_path,results = os.path.join(results_dir, run_id,"detailed.txt"),{"detailed": ""}		 
		if os.path.exists(detailed_path):
			with open(detailed_path) as f:
				results["detailed"] = f.read()
		return simplejson.dumps(results)

	@cherrypy.expose
	def get_progress(self, run_id):
		# Loads the progress file if it exists
		global results_dir, uploads_dir, sett
		p = {"status":"","curprog":0,"progmax":0}
		progress_path = os.path.join(os.path.join(results_dir, run_id),".prog")
		if os.path.exists(progress_path):
			with open(progress_path) as f:
				p = json.loads(f.read())
		return simplejson.dumps(p)

	@cherrypy.expose
	def get_log(self,run_id):
		global results_dir, uploads_dir, sett
		results = {"log": ""}
		log_path = os.path.join(os.path.join(results_dir, run_id),"gr_log.txt")
		if os.path.exists(log_path):
			with open(log_path) as f:
				results["log"] = f.read()
		return simplejson.dumps(results)

	@cherrypy.expose
	def get_gfdescriptions(self,organism,db_version):
		return open(os.path.join(sett["data_dir"][db_version],organism,"gf_descriptions.txt")).read()
	@cherrypy.expose
	def get_checkboxtree(self,organism,db_version):
		return open(os.path.join(sett["data_dir"][db_version],organism,"treeview.json")).read()


	@cherrypy.expose
	def enrichment_log(self, id):
		global results_dir
		with open(os.path.join(results_dir,id+".log")) as sr:
			x = sr.read()
			return "<p>{}</p>".format(x.replace("\n","<br/>"))

	@cherrypy.expose
	def cite(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("cite.mako").render(),script="")

	@cherrypy.expose
	def news(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("news.mako").render(),script="")

	@cherrypy.expose
	def overview(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("overview.mako").render(),script="")

	@cherrypy.expose
	def demo(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("demo.mako").render(),script="")

	@cherrypy.expose
	def roadmap(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("roadmap.mako").render(),script="")

	@cherrypy.expose
	def help(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("help.mako").render(),script="")

def base_name(k):
    return os.path.basename(k).split(".")[0]

def verify_score_strand(gf_path,pct_score,strand,data_dir):
    ''' Checks if a score and/or strand filtered version of gf_path exists in the database and
    returns the appropriate path if it does.
    '''
    global results_dir, uploads_dir, sett
    gf_path = os.path.join(data_dir,gf_path.lstrip("/"))
    gf_score_strand_path = gf_path.replace('/grsnp_db/','/grsnp_db_{}_{}/'.format(pct_score,strand))
    gf_score_path = gf_path.replace('/grsnp_db/','/grsnp_db_{}/'.format(pct_score))
    gf_strand_path = gf_path.replace('/grsnp_db/','/grsnp_db_{}/'.format(strand))
    if os.path.exists(gf_score_strand_path):
        return gf_score_strand_path
    elif os.path.exists(gf_score_path):
        return gf_score_path
    elif os.path.exists(gf_strand_path):
    	return gf_strand_path
    else:
    	return gf_path

def main():
	global sett, results_dir, uploads_dir
	root_dir = os.path.dirname(os.path.realpath(__file__))
	static_dir = os.path.abspath(os.path.join(root_dir, "frontend/static"))
	media = os.path.abspath(os.path.join(".","frontend/media"))
	parser = argparse.ArgumentParser(prog="python -m grsnp.server", description="Starts the GenomeRunner SNP server. Example: python -m grsnp.server -d /home/username/db_#.##_#.##.####/ -g hg19 -p 8000", epilog="Use GenomeRunner SNP: http://localhost:8000/gr")
	parser.add_argument("--data_dir" , "-d", nargs="?",type=str, help="Set the directory containing the database. Required. Use absolute path. Example: /home/username/db_#.##_#.##.####/.", required=True)
	parser.add_argument("--run_files_dir" , "-r", nargs="?", help="Set the directory where the server should save results. Required. Use absolute path. Example: /home/username/run_files/.", required=True)
	parser.add_argument("--organism" , "-g", nargs="?", help="The UCSC code for the organism to use. Default: hg19 (human). Data for the organism must exist in the database directory. Use dbcreator to make the database, if needed.", default="hg19")
	parser.add_argument("--port","-p", nargs="?", help="Socket port to start server on. Default: 8000", default=8080) 
	parser.add_argument("--num_workers", "-w", type=int, help="The number of celery workers to start. Default: 1", default=1)
	parser.add_argument("--group", "-z", type=str, help="The group to change results folder permission to", default="")	

	args = vars(parser.parse_args())
	port = args["port"]
	list_data_dir = args["data_dir"].split(",")

	if list_data_dir == "":
		print "ERROR: data_dir is a required argument"
		sys.exit()

	if args['run_files_dir'] == "":
		print "ERROR: run_files_dir is a required argument"
		sys.exit()
	data_dir = {}
	for db_dir in list_data_dir:
		if db_dir.strip()[-1] == "/":
			data_dir.update({os.path.split(db_dir[:-1])[1]:os.path.join(db_dir.strip(),"grsnp_db")})
		else:
			data_dir.update({os.path.split(db_dir)[1]:os.path.join(db_dir.strip(),"grsnp_db")})

	# setup the logging
	hdlr = logging.FileHandler(os.path.join(args['run_files_dir'], 'genomerunner_server.log'))
	hdlr_std = StreamHandler()
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.addHandler(hdlr_std)
	logger.setLevel(logging.INFO)

	# global settings used by GR
	sett = {"data_dir":data_dir,
				"run_files_dir": args["run_files_dir"],
				"default_organism":args["organism"],
				"group": args['group']	
				}

	
	#validate data directory
	for k,v in data_dir.items():
		if not os.path.exists(v):
			print "ERROR: {} does not exist. Please run grsnp.dbcreator or use --data_dir.".format(v)
			sys.exit()
		if not os.path.exists(os.path.join(v,sett["default_organism"])):
			print "ERROR: Database for default organism {} does not exist. Either change the default organism or add data for that organism to the database at {} using the bedfilecreator".format(sett["default_organism"],v)
			sys.exit()

	# validate run_files directory
	if not os.path.exists(sett["run_files_dir"]): os.mkdir(sett["run_files_dir"])
	results_dir = os.path.join(sett["run_files_dir"],"results")
	uploads_dir = os.path.join(sett["run_files_dir"],"uploads")	
	if not os.path.exists(results_dir):
		os.mkdir(results_dir)
	if not os.path.exists(uploads_dir):
		os.mkdir(uploads_dir)
	if port:
		cherrypy.server.max_request_body_size = 0
		cherrypy.config.update({
			"server.socket_port": int(port),
			"server.socket_host":"0.0.0.0"})
		conf = {"/static": 
					{"tools.staticdir.on": True,
					"tools.staticdir.dir": static_dir},
				"/results": 
					{"tools.staticdir.on": True,
					"tools.staticdir.dir": os.path.abspath(results_dir)}
				}
		# gather all of the data directories and make them available as static directories
		for k,v in data_dir.items():
			conf.update({"/" + k: 
							{"tools.staticdir.on": True,
							"tools.staticdir.dir": v}}) 
		
		# start redis server
		script = ["redis-server", "--port", str(celeryconfiguration.redis_port)]
		fh = open(os.path.join(sett["run_files_dir"],"redis.log"),"w")
		out = subprocess.Popen(script,stdout=fh,stderr=fh)
		#script = "ps auxww | grep  -E 'worker.*grsnp_LOCAL' | awk '{print $2}' | xargs kill -9"
		#out = subprocess.Popen(script,shell=True)
		#out.wait()
		app = Celery('grsnp')
		app.config_from_object('grsnp.celeryconfiguration')
		print "Checking for existing celery workers..."
		if app.control.inspect().ping() == None:
			for i in range(args["num_workers"]):
				print "Starting Celery worker[s]..."
				fh = open("worker{}.log".format(i),"w")
				script = ["celery","worker", "--app", "grsnp.worker_hypergeom4", "--loglevel", "INFO", "-n", "grsnp_LOCAL{}.%h".format(i),'-r',sett['run_files_dir'],'-d',args["data_dir"]]
				out = subprocess.Popen(script,stdout=fh,stderr=fh)
		else:
			workers = app.control.inspect().ping()
			pids = [str(app.control.inspect().stats()[j]['pid']) for j in workers.keys()]
			print "Celery workers already running. Pids:" + ",".join(pids) 
		print "Redis backend URL: ", celeryconfiguration.CELERY_RESULT_BACKEND
		cherrypy.config.update({'tools.sessions.timeout': 60})
		cherrypy.quickstart(WebUI(), "/", config=conf)

	else:
		print "WARNING: No port given. Server not started. Use --port flag to set port."


if __name__ == "__main__":	
	main()

