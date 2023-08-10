import argparse
import os
import sys
import yaml
import glob
import subprocess

def parse_arguments():
  parser = argparse.ArgumentParser(description = "Screen assemblies for junction sequences in order to characterize specific plasmid markers")
  parser.add_argument("-i", metavar = "--indir", dest = "indir", help = "Input path to assembly directory", required = True)
  parser.add_argument("-q", metavar = "--query_file", dest = "query_file", help = "Path to Query file", required = True)
  parser.add_argument("-T", metavar = "--theshold", dest = "threshold", help = "Cutoff threshold of match coverage and identity. Ignore threshold by setting to 0 or False. (Default 99.9)", default = 99.9)
  parser.add_argument("-b", dest = "top_only", help = "Include the best results, meaning only the top hit for each sample and their queries. (Default False)", action = "store_true")
  parser.add_argument("-o", metavar = "--outdir", dest = "outdir", help = "Output path to Results and Temporary files directory", required = True)
  parser.add_argument("-e", dest = "include_empty", help = "Include samples with not hits in the junctions results files. (Default False)", action = "store_true")
  parser.add_argument("-k", dest = "keep_config", help = "Preserve configuration file and exit pipeline for debugging purposes. (Default False)", action = "store_true")
  parser.add_argument("-t", metavar = "--threads", dest = "threads", help = "Number of threads to allocate for the pipeline. (Default 1)", default = 1)
  parser.add_argument("-F", dest = "force", help = "Force rerun of all tasks in pipeline. (Default False)", action = "store_true")
  parser.add_argument("-S", dest = "skipmake", help = "Skip Snakemake for requirering manual run of Snakemake. Config file will be generated.", action = "store_true")
  parser.add_argument("-n", dest = "dry_run", help = "Perform a dry run with Snakemake to see jobs but without executing them. (Default False)", action = "store_true")
  parser.add_argument("-d", dest = "debug", help = "Enable debug mode, stores snakemake object for inspection in R. (Default False)", action = "store_true")

  return(parser.parse_args())


def validate_snakemake(debug):
  here = os.listdir('.')
  workflow_here = 'workflow' in os.listdir('.')

  if (workflow_here):
    snakefile_here = "Snakefile" in os.listdir('workflow')

    if snakefile_here and debug:
      print("Snakefile detected")
    elif not snakefile_here:
      print('Error no Snakefile detected in workflow directory. Software is properbly corrupt, consider redownloading.')
      sys.exit(1)
  else:
    print('No workflow directory detected, are you sure you are running the script from the software folder?')
    sys.exit(1)


def generate_configfile(indir, query_file, threshold, top_only, outdir, include_empty, threads, debug, keep_config, skipmake):
  # Define config file
  config_file = "config/config.yaml"
  
  # Check for existing config file
  if not os.path.isfile(config_file):
    config_dir = os.path.dirname(config_file)
    
    # Check for existing config dir
    if not os.path.isdir(config_dir):
      print("No config dir detected, creating directory.")
      os.mkdir(config_dir)
  elif keep_config and not skipmake:
    print("--keep_config is set to true, exiting!")
    sys.exit(0)

  # Directing full paths
  in_path = os.path.abspath(indir).rstrip("/")
  out_path = os.path.abspath(outdir).rstrip("/")
  query_path = os.path.abspath(query_file).rstrip("/")
  
  config = {"indir" : in_path, "query_file" : query_path, "threshold" : threshold,  "outdir" : out_path, "top_only" : top_only, "include_empty" : include_empty, "debug" : debug}

  with open(config_file, "w") as config_yaml:
    yaml.dump(config, config_yaml)

# Derrive arguments
args = parse_arguments()

indir = args.indir
query_file = args.query_file
threshold = args.threshold
top_only = args.top_only
outdir = args.outdir
include_empty = args.include_empty
keep_config = args.keep_config
threads = args.threads
force = args.force
skipmake = args.skipmake
dry_run = args.dry_run
debug = args.debug

# Validate snakemake structure
validate_snakemake(debug)

# Prepare config file for snakemake
generate_configfile(indir, query_file, threshold, top_only, outdir, include_empty, threads, debug, keep_config, skipmake)

# Determine sample files
sample_files = glob.glob(pathname = "%s/*fasta" %indir)
samples = [os.path.basename(sample_file.rstrip("\.fastq\.gz")) for sample_file in sample_files]

if skipmake:
  print("Warning: Skipping Snakemake!")
else:
  snake_args = ""
  if force:
    snake_args += " -F "
  if dry_run:
    snake_args += " -n "

  snakemake_cmd = "snakemake --use-conda --cores %s%s" %(threads, snake_args) 
  if debug:
    print("Running command: %s" %snakemake_cmd)

  subprocess.Popen(snakemake_cmd, shell = True).wait()

