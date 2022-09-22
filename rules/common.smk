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

# import config file
configfile: "config.yaml"
#snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")    # Validate the format of config file

# store samples in config['global']
if "global" in config:
    raise Exception("Config key 'global' already defined. Someone messed with our setup.")
else:
    config["global"] = {}

config["global"]["samples"] = pd.read_csv(
    config["data"]["samples"], sep='\t', dtype=str).set_index(["sample"], drop=False
)
#snakemake.utils.validate( config["global"]["samples"], schema="../schemas/samples.schema.yaml" )

config["global"]["sample-names"] = list()
config["global"]["control-sample-names"] = list()
for index, row in config["global"]["samples"].iterrows():
    s = row["sample"]
    if s not in config["global"]["sample-names"] and not s.startswith('control'):
        config["global"]["sample-names"].append(s)
    else:
        config["global"]["control-sample-names"].append(s)

def valid_filename(fn):
    # Only accept alnum, underscore, dash, and dot.
    return fn.replace('_', '').replace('-', '').replace('.', '').isalnum() and fn.isascii()

def valid_filepath(fn):
    # Only accept alnum, underscore, dash, dot, and slashes.
    clean = fn.replace('_', '').replace('-', '').replace('.', '').replace('/', '').replace('\\', '')
    return clean.isalnum() and clean.isascii()

wildcard_constraints:
    sample="|".join(config["global"]["sample-names"]),
    control-sample="|".join(config["global"]["control-sample-names"])

print('samples:', config["global"]["sample-names"])
print('control samples:', config["global"]["control-sample-names"])

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = config["global"]["all-samples"].loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}


def is_single_end( sample, **kargs ):
    """Return True if sample-unit is single end."""
    return pd.isnull(config["global"]["all-samples"].loc[(sample), "fq2"])