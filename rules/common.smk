# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd
import os, sys, pwd
import socket, platform
import subprocess
from datetime import datetime

# Ensure min Snakemake version
snakemake.utils.min_version("5.7")
basedir = workflow.basedir

cnvpipe_version = "0.0.1" #CNVPIPE_VERSION#

# Import config file
configfile: "config.yaml"

# Validate the format of config file TODO
#snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")

# Clean up the file name for the reference genome.
if config["data"]["genome"].endswith(".gz"):
    config["data"]["genome"] = os.path.splitext(config["data"]["genome"])[0]

# Store sample names into config['global']
if "global" in config:
    raise Exception("Config key 'global' cannot be defined in config file, please carefully check.")
else:
    config["global"] = {}

config["global"]["samples"] = pd.read_csv(
    config["data"]["samples"], sep='\t', dtype=str).set_index(["sample"], drop=False)

# Validate the format of sample tsv file TODO
#snakemake.utils.validate( config["global"]["samples"], schema="../schemas/samples.schema.yaml" )

# All control samples in sample tsv should start with 'control' keyword.
config["global"]["sample-names"] = list()
config["global"]["control-sample-names"] = list()
for index, row in config["global"]["samples"].iterrows():
    s = row["sample"]
    if s not in config["global"]["sample-names"] and not s.startswith('control'):
        config["global"]["sample-names"].append(s)
    else:
        config["global"]["control-sample-names"].append(s)

config["global"]["all-sample-names"] = config["global"]["sample-names"] + config["global"]["control-sample-names"]

def valid_filename(fn):
    # Only accept alnum, underscore, dash, and dot.
    return fn.replace('_', '').replace('-', '').replace('.', '').isalnum() and fn.isascii()

def valid_filepath(fn):
    # Only accept alnum, underscore, dash, dot, and slashes.
    clean = fn.replace('_', '').replace('-', '').replace('.', '').replace('/', '').replace('\\', '')
    return clean.isalnum() and clean.isascii()

# Check whether filename and filepath are valid
problematic_filenames = 0
for index, row in config["global"]["samples"].iterrows():
    if not valid_filename(row['sample']):
        raise Exception("Invalid sample name: "+ str(row["sample"]) + "we only allow alpha-numerical, dots, dashes, and underscores.")
    if not os.path.isfile(row['fq1']) or (not pd.isnull(row['fq2']) and not os.path.isfile(row['fq2'])):
        raise Exception(
            "Input fastq files listed in the input files table " + config["data"]["samples"] +
            " not found: " + str(row["fq1"]) + "; " + str(row["fq2"])
        )
    if not valid_filepath(row['fq1']) or (not pd.isnull(row['fq2']) and not valid_filepath(row['fq2'])):
        problematic_filenames += 1


# =================================================================================================
#     Pipeline User Output
# =================================================================================================

# Get a nicely formatted username and hostname
username = pwd.getpwuid( os.getuid() )[ 0 ]
hostname = socket.gethostname()
hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")

# Get the conda version, if available.
try:
    process = subprocess.Popen(['conda', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode('ascii')
    conda_ver = out[out.startswith("conda") and len("conda"):].strip()
    del process, out, err
    if not conda_ver:
        conda_ver = "n/a"
except:
    conda_ver = "n/a"

# Get the conda env name, if available.
# See https://stackoverflow.com/a/42660674/4184258
conda_env = os.environ['CONDA_DEFAULT_ENV'] + " (" + os.environ["CONDA_PREFIX"] + ")"
if conda_env == " ()":
    conda_env = "n/a"

# Get nicely wrapped command line
cmdline = sys.argv[0]
for i in range( 1, len(sys.argv)):
    if sys.argv[i].startswith("--"):
        cmdline += "\n                        " + sys.argv[i]
    else:
        cmdline += " " + sys.argv[i]

# Get abs paths of all config files
cfgfiles = []
for cfg in workflow.configfiles:
    cfgfiles.append( os.path.abspath(cfg) )
cfgfiles = "\n                        ".join(cfgfiles)

smpcnt = str(len(config['global']['sample-names']))
ctrsmpcnt = str(len(config['global']['control-sample-names']))

# Some helpful messages
logger.info("=====================================================================================")
logger.info("       ______  __   __ __        __ ______  ___   ______   _______ ")
logger.info("      /   ___\/  \ /  /\ \      / /|   _  \ \  \ |   _  \ /  ____/ ")
logger.info("     /   /    |   \|  | \ \    / / |  |_]  ||  | |  |_]  |  |___   ")
logger.info("    |   /     |       |  \ \  / /  |   ___/ |  | |   ___/|   ___|  ")
logger.info("     \  \_____|  |\   |   \ \/ /   |  |     |  | |  |    |  |____  ")
logger.info("      \______/|__| \__|    \__/    |__|     \__\ |__|    \_______\ ")
logger.info("")
logger.info("    Date:               " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
logger.info("    User:               " + username)
logger.info("    Host:               " + hostname)
logger.info("    Conda:              " + str(conda_ver))
logger.info("    Python:             " + str(sys.version.split(' ')[0]))
logger.info("    Snakemake:          " + str(snakemake.__version__))
logger.info("    CNVPipe:            " + str(cnvpipe_version))
logger.info("    Conda env:          " + str(conda_env))
logger.info("    Command:            " + cmdline)
logger.info("")
logger.info("    Base directory:     " + workflow.basedir)
logger.info("    Working directory:  " + os.getcwd())
logger.info("    Config file(s):     " + cfgfiles)
logger.info("    Samples:            " + smpcnt)
logger.info("    ControlSamples:     " + ctrsmpcnt)
logger.info("")
logger.info("=====================================================================================")
logger.info("")

# Warning about input names.
if problematic_filenames > 0:
    logger.warning(
        "In " + str(problematic_filenames) + " of the " + str(len(config["global"]["sample-names"])) +
        " input fastq files listed in the input files table " + config["data"]["samples"] +
        " contain problematic characters. We generally advise to only use alpha-numeric " +
        "characters, dots, dashes, and underscores. "
    )

# No need to have these output vars available in the rest of the snakefiles
del problematic_filenames
del username
del hostname
del conda_ver
del conda_env
del cmdline
del cfgfiles
del smpcnt
del ctrsmpcnt
