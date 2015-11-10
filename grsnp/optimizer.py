import sys
import logging
from logging import FileHandler,StreamHandler
import os
import gzip
import subprocess
import argparse
from grsnp import worker_optimizer as worker_opt
import pdb
from celery import group
import celeryconfiguration_optimizer
from time import sleep
import atexit

# connection information for the ucsc ftp server
logger = logging.getLogger()




def create_bkg_gf_overlap_db(gf_dir,background_dir,data_dir):
	""" Used to precalculate the overlapStatistics for the GFs against each of the
	default backgrounds.
	"""

	gf_bg_stats,list_completed_gfs = {},[]
	all_gfs = []
	backgrounds,gfs=[],[]
	db_path = os.path.join(data_dir,gf_dir,"bkg_overlaps.gr")
	full_gf_dir = os.path.join(data_dir,gf_dir)
	full_background_dir = os.path.join(data_dir,background_dir)
	print db_path

	# Read in all completed GF in the partially completed bkg_overlaps file
	if os.path.exists(db_path):
		list_completed_gfs = [x.split("\t")[0] for x in open(db_path).read().split("\n")]
		print [x for x in list_completed_gfs if "wgEncodeHaibTfbsEcc1Tcf12V0422111PkRep2" in x]
	# gather all directories (groups) in the database
	dirs = [name for name in os.listdir(full_gf_dir)
		if os.path.isdir(os.path.join(full_gf_dir, name))]

	# Gather backgrounds paths
	backgrounds = [os.path.join(background_dir, f) for f in os.listdir(full_background_dir) if f.endswith(('.gz', '.bb',".txt",".bed"))]

	cur_prog,prog_max = 1,_count_gfs(full_gf_dir)
	# Process each category of GFs
	for d in dirs:
		# Gather gfs paths
		gfs = []
		for base, d, files in os.walk(os.path.join(full_gf_dir,d)):
				partial_base = base.replace(data_dir,"") # get the relative path of the database
				gfs += [os.path.join(partial_base, f) for f 
					in files if f.endswith(('.gz', '.bb'))]	
		all_gfs += [x for x in gfs if x not in list_completed_gfs]
		prog_gf  = 1
	# Run overlap analysis using Celery
	logger.info("Running overlapStatistics for all GFs in {}".format(gf_dir))
	results = group(worker_opt.calculate_bkg_gf_overlap.s(gf_path=g,list_bkg_paths=backgrounds).set(queue='optimizer.group') for g in all_gfs)()
	while not results.ready():
		sys.stdout.write("{} of {} completed...\r".format(results.completed_count(),len(all_gfs)))
		sys.stdout.flush()
		sleep(5.0)

	results = results.join()
	write_results(results,db_path)

def write_results(results,outputpath):
	''' Results are written out to a temporary file which replaces the existing file if after successfully
	being written.

	NOTES:
	Results data structure is as follows:
		[{gf1_path: [{background1_overlapstatistics},{background2_overlapstatistics} ...]},
		 {gf2_path: [{background1_overlapstatistics},{background2_overlapstatistics} ...]} ...
		]
	'''

	with open(outputpath,'a') as writer:
		for res in results:
			if isinstance(res,str):
				logger.error(res)
				continue
			gf = res.keys()[0]
			stats = res[gf]
			stat_line = [x["queryfile"]+":"+str(x["intersectregions"])+":"+str(x["queryregions"]) for x in stats]
			stat_line = ",".join(stat_line) + "\n"
			writer.write(gf+"\t"+stat_line)

def _count_gfs(grsnp_db):
	x = 0
	for root, dirs, files in os.walk(grsnp_db):
		for f in files:
			if f.endswith(".bed.gz"):
				x = x+1
	return x

def shutdown_workers():
	# kill all existing optimizer workers
	print "Stopping local workers..."
	script = "ps auxww | grep  -E '*grsnp_optimizerLOCAL*' | awk '{print $2}' | xargs kill -9"
	out = subprocess.Popen(script,shell=True)
	out.wait()
	print "Removing leftover optimizer jobs from Celery queue 'optimizer.group'..."
	script = "celery amqp queue.purge optimizer.group --broker " +  celeryconfiguration_optimizer.CELERY_RESULT_BACKEND
	out = subprocess.Popen(script,shell=True)
	out.wait()


def main():
	parser = argparse.ArgumentParser(prog="python -m grsnp.optimizer", description="""Pre calculates the overlapStatistics for each of the backgrounds in <db_path>/custom_data/backgrounds/<organism> and genomic features in <db_path>/grsnp_db/<organism>. Example: python -m grsnp.optimizer -d /home/username/grs_db/ -g mm9""", epilog="""Creates a file  <db_path>/grsnp_db/<organism>/bkg_overlap.gr, automatically used by the server to speed up the analyses""")
	parser.add_argument('--data_dir','-d', nargs="?", help="Set the directory containing the database. Required. Use absolute path. Example: /home/username/db_2.00_6.26.2014/.", required=True)
	parser.add_argument('--organism','-g', nargs="?", help="The UCSC code for the organism to use. Default: hg19 (human). Data for the organism must exist in the database directory. Use dbcreator to make the database, if needed.", required=True, default="hg19")
	parser.add_argument("--num_workers", "-w", type=int, help="The number of local celery workers to start. Default: 1", default=1)	
	args = vars(parser.parse_args())

	hdlr = logging.FileHandler(os.path.join(args["data_dir"],'genomerunner_optimizer.log'))
	hdlr_std = StreamHandler()
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	hdlr_std.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.addHandler(hdlr_std)
	logger.setLevel(logging.INFO)

	atexit.register(shutdown_workers) # shutdown workers on termination
	logger.info('Running optimization')
	# start redis server
	script = ["redis-server", "--port", str(celeryconfiguration_optimizer.redis_port)]
	fh = open("redis.log","w")
	out = subprocess.Popen(script,stdout=fh,stderr=fh)
	for i in range(args["num_workers"]):
		fh = open("worker{}.log".format(i),"w")
		script = ["celery","worker", "--app", "grsnp.worker_optimizer","--loglevel", "INFO", '-Q','optimizer.group', "-n", "grsnp_optimizerLOCAL{}.%h".format(i),'--data_dir',args['data_dir']]
		out = subprocess.Popen(script,stdout=fh,stderr=fh)
	print "Redis backend URL: ", celeryconfiguration_optimizer.CELERY_RESULT_BACKEND

	# find all the folders with GF data including those filtered by score
	grdb_dirs = [os.path.join(args["data_dir"],name) for name in os.listdir(args["data_dir"])
			if os.path.isdir(os.path.join(args["data_dir"], name)) and "grsnp_db" in name]
	for gr_dir in grdb_dirs:

		background_dir = os.path.join("custom_data","backgrounds",args["organism"])
		gfs_dir = os.path.join(os.path.split(gr_dir)[1],args["organism"])
		if not os.path.exists(os.path.join(args['data_dir'], gfs_dir)):
			continue
		if not os.path.exists(os.path.join(args['data_dir'],background_dir)):
			print "ERROR: No backgrounds found in default background directory {}.  Please add backgrounds.".format(background_dir)
			sys.exit()
		logger.info("Pre calculating statistics for GR database in")
		create_bkg_gf_overlap_db(gf_dir=gfs_dir,background_dir=background_dir,data_dir=args['data_dir'])
	logger.info("Finished creating optimization files.")


if __name__ == "__main__":
	main()