from celery import Celery
from celery import signals
import grsnp.hypergeom4
from celery.bin.base import Option
from celery.exceptions import MaxRetriesExceededError #,Reject
import os
import traceback
import json

# celery
app = Celery('grsnp')
app.config_from_object('grsnp.celeryconfiguration')
app.user_options['preload'].add(
    Option('-d', '--data_dir', default='',
           help='Set the directory containing the database. Pass multiple directories for multiple database versions. Required. Use absolute path. Example: /home/username/db_#.##_#.##.####/.'),
)
app.user_options['preload'].add(
    Option('-r', '--run_files_dir', default='',
           help="Set the directory where the server should save results. Required. Use absolute path. Example: /home/username/run_files/."),
)

sett = {}

@app.task(ignore_result=False)
def run_hypergeom(fois, gfs, bg_path,job_name="",zip_run_files=False,bkg_overlaps_path="",run_annotation=False,run_randomization_test=False,pct_score="",organism="",id="",db_version=None,stat_test=None):
	global sett
	try:
		if db_version not in sett['root_data_dir'].keys():
			raise Exception("{} does not exist in the worker's database. Databases available: {}".format(db_version,",".join(sett['root_data_dir'].keys())))
		outdir=os.path.join(sett['run_files_dir'],'results',str(id))	
		# write out absolute gfs and fois file paths and pass these to the worker	
		for f_path in [fois, gfs]:
			list_f = [x for x in open(f_path).read().split("\n") if x!= ""]
			with open(f_path+'_full','wb') as writer:
				for f in list_f:
					if f.lstrip("/").startswith("grsnp_db") or f.lstrip("/").startswith("custom_data"):
						writer.write(os.path.join(sett['root_data_dir'][db_version],f.lstrip('/'))+"\n")
					else:
						writer.write(os.path.join(sett['run_files_dir'],f.lstrip("/"))+"\n")
		# make background path absolute
		if bg_path.lstrip("/").startswith("custom_data"):
			bg_path = os.path.join(sett['root_data_dir'][db_version],bg_path.lstrip("/"))
		else:
			bg_path = os.path.join(sett['run_files_dir'],bg_path.lstrip("/"))
		print "Worker starting job for {}".format(id)
		grsnp.hypergeom4.run_hypergeom(fois+"_full", gfs+"_full", bg_path,outdir,job_name,zip_run_files,bkg_overlaps_path,sett['root_data_dir'][db_version],run_annotation,run_randomization_test,pct_score,organism,stat_test=stat_test)
	except Exception, e:
		print traceback.print_exc()
		_write_progress("ERROR: Run crashed. Celery worker threw an error.",id,1,1)
		raise e


# process command line arguments if they exist
@signals.user_preload_options.connect
def cmd_options(options,**kwargs):
	# These settings set the location of the database that the celery worker should use.
	# The data_dir can be in a different location from that of the server's data_dir provided. Make sure to pass all version of the database that are used by the server.
	# it has the exact same GFs as the server.
	# run_files_dir MUST point to the same database that the server is using.
	global sett
	list_data_dir = options["data_dir"].split(",")
	if options['data_dir'] == '':
		raise Exception('data_dir is a required argument')	
	if options['run_files_dir'] == '':
		raise Exception('run_files_dir is a required argument')
	if not os.path.exists(options['run_files_dir']):
		raise Exception('{} does not exist'.format(options['run_files_dir']))

	root_data_dir = {}
	for db_dir in list_data_dir:
		# Extract db_version and use as key. Value is directory to db.
		root_data_dir.update({os.path.split(db_dir.rstrip("/").lstrip("/").strip())[1]:db_dir.strip().rstrip("/")})
	#validate data directory
	for k,v in root_data_dir.items():
		if not os.path.exists(v):
			raise Exception("{} does not exist. Please run grsnp.dbcreator or use --data_dir.".format(v))
	# save to global settings
	sett["run_files_dir"] = options["run_files_dir"].rstrip("/")
	sett["root_data_dir"] = root_data_dir

def _write_progress(line,id,curprog,progmax):
    """Saves the current progress to the progress file
    """
    global sett
    progress_outpath = os.path.join(sett['run_files_dir'],'results',id,".prog")
    if progress_outpath:
        progress = {"status": line, "curprog": curprog,"progmax": progmax}
        with open(progress_outpath,"wb") as progfile:
            progfile.write(json.dumps(progress))
