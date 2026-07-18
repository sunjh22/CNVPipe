#!/usr/bin/env python3
"""PBS/Torque cluster status script for Snakemake."""

import sys
import subprocess
import xml.etree.cElementTree as ET

# Validate job ID (only allow alphanumeric, dots, and brackets for array jobs)
jobid = sys.argv[1]
if not all(c.isalnum() or c in '.-[]' for c in jobid):
    print("failed")
    sys.exit(0)

try:
    res = subprocess.run(["qstat", "-f", "-x", jobid], check=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    xmldoc = ET.ElementTree(ET.fromstring(res.stdout.decode())).getroot()
    job_state = xmldoc.findall('.//job_state')[0].text

    if job_state == "C":
        exit_status = xmldoc.findall('.//exit_status')[0].text
        if exit_status == '0':
            print("success")
        else:
            print("failed")
    else:
        print("running")

except (subprocess.CalledProcessError, IndexError) as e:
    print("failed")
