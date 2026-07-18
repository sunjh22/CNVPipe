#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess

from snakemake.utils import read_job_properties

parser=argparse.ArgumentParser(add_help=False)
parser.add_argument("--depend", help="Space separated list of ids for jobs this job should depend on.")
parser.add_argument("-a", help="Declare the time when the job becomes eligible for execution.")
parser.add_argument("-A", help="Define the account string.")
parser.add_argument("-b", help="PBS Server timeout.")
parser.add_argument("-c", help="Checkpoint options.")
parser.add_argument("-C", help="Directive prefix in script file.")
parser.add_argument("-d", help="Working directory to be used (default: ~). PBS_O_INITDIR")
parser.add_argument("-D", help="Root directory to be used. PBS_O_ROOTDIR")
parser.add_argument("-e", help="standard error path.")
parser.add_argument("-f", help="Fault tolerant.",action="store_true")
parser.add_argument("-h", help="Apply user hold at submission time",action="store_true")
parser.add_argument("-j", help="Merge standard error and standard out. (oe or eo)")
parser.add_argument("-l", help="Resource list.")
parser.add_argument("-m", help="Mail options.")
parser.add_argument("-M", help="Mail users.")
parser.add_argument("-N", help="Name for the job.")
parser.add_argument("-o", help="standard output path.")
parser.add_argument("-p", help="Set job priority.")
parser.add_argument("-P", help="Proxy user for job.")
parser.add_argument("-q", help="Set destination queue.")
parser.add_argument("-t", help="Array request.")
parser.add_argument("-u", help="Set user name for job.")
parser.add_argument("-v", help="Environment variables to export to the job.")
parser.add_argument("-V", help="Export all environment variables.",action="store_true")
parser.add_argument("-w", help="Set working directory. PBS_O_WORKDIR")
parser.add_argument("-W", help="Additional attributes.")
parser.add_argument("--help", help="Display help message.",action="store_true")

parser.add_argument("positional",action="append",nargs="?")
args = parser.parse_args()

if args.help :
    parser.print_help()
    sys.exit(0)

jobscript = sys.argv[-1]

job_properties = read_job_properties(jobscript)

atime=""
acc_string=""
pbs_time=""
chkpt=""
pref=""
dd=""
rd=""
se=""
ft=""
hold=""
j=""
resource=""
mail=""
mailuser=""
jname=""
so=""
priority=""
proxy=""
q=""
ar=""
user=""
ev=""
eall=""
wd=""
add=""
depend=""
resourceparams=""
extras=""


if args.depend:
	depend_ids = [m for m in args.depend.split() if m]
	if depend_ids:
		depend = "depend=afterok:" + ":".join(depend_ids)

if args.positional:
	for m in args.positional:
		extras = extras + " " + m

if args.a: atime = " -a " + args.a
if args.A: acc_string = " -A " + args.A
if args.b: pbs_time = " -b " + args.b
if args.c: chkpt = " -c " + args.c
if args.C: pref = " -C " + args.C
if args.d: dd = " -d " + args.d
if args.D: rd = " -D " + args.D
if args.e: se = " -e " + args.e
if args.f: ft = " -f"
if args.h: hold = " -h"
if args.j: j = " -j " + args.j
if args.l: resource = " -l " + args.l
if args.m: mail = " -m " + args.m
if args.M: mailuser = " -M " + args.M
if args.N: jname = " -N " + args.N
if args.o: so = " -o " + args.o
if args.p: priority = " -p " + args.p
if args.P: proxy = " -P " + args.P
if args.q: q = " -q " + args.q
if args.t: ar = " -t " + args.t
if args.u: user = " -u " + args.u
if args.v: ev = " -v " + args.v
if args.V: eall = " -V"
if args.w: wd = " -w " + args.w
if args.W: add = args.W

nodes=""
ppn=""
mem=""
walltime=""

if "threads" in job_properties:
    ppn = "ppn=" + str(job_properties["threads"])

if "resources" in job_properties:
    resources = job_properties["resources"]
    if "nodes" in resources: nodes="nodes=" + str(resources["nodes"])
    if ppn and not nodes : nodes="nodes=1"
    if "mem" in resources: mem="mem=" + str(resources["mem"])
    if "walltime" in resources: walltime="walltime=" + str(resources["walltime"])

if nodes or ppn or mem or walltime: resourceparams = " -l \""
if nodes: resourceparams = resourceparams + nodes
if nodes and ppn: resourceparams = resourceparams + ":" + ppn
if nodes and mem: resourceparams = resourceparams + ","
if mem: resourceparams = resourceparams + mem
if walltime and (nodes or mem): resourceparams = resourceparams + ","
if walltime: resourceparams = resourceparams + walltime
if nodes or mem or walltime: resourceparams = resourceparams + "\""

if "cluster" in job_properties:
    cluster = job_properties["cluster"]
    if "error" in cluster:
        os.makedirs(os.path.dirname(cluster["error"]), exist_ok=True)
        se = " -e " + cluster["error"]
    if "output" in cluster:
        os.makedirs(os.path.dirname(cluster["output"]), exist_ok=True)
        so = " -o " + cluster["output"]

# Build argument list safely (no shell=True to avoid injection)
cmd_args = ["qsub"]
if atime: cmd_args.extend(["-a", args.a])
if acc_string: cmd_args.extend(["-A", args.A])
if pbs_time: cmd_args.extend(["-b", args.b])
if chkpt: cmd_args.extend(["-c", args.c])
if pref: cmd_args.extend(["-C", args.C])
if dd: cmd_args.extend(["-d", args.d])
if rd: cmd_args.extend(["-D", args.D])
if se: cmd_args.extend(["-e", args.e])
if ft: cmd_args.append("-f")
if hold: cmd_args.append("-h")
if j: cmd_args.extend(["-j", args.j])
if resource: cmd_args.extend(["-l", args.l])
if mail: cmd_args.extend(["-m", args.m])
if mailuser: cmd_args.extend(["-M", args.M])
if jname: cmd_args.extend(["-N", args.N])
if so: cmd_args.extend(["-o", args.o])
if priority: cmd_args.extend(["-p", args.p])
if proxy: cmd_args.extend(["-P", args.P])
if q: cmd_args.extend(["-q", args.q])
if ar: cmd_args.extend(["-t", args.t])
if user: cmd_args.extend(["-u", args.u])
if ev: cmd_args.extend(["-v", args.v])
if eall: cmd_args.append("-V")
if wd: cmd_args.extend(["-w", args.w])
if add: cmd_args.extend(["-W", add])
# Resource list from job properties
if nodes or ppn or mem or walltime:
    res_list = []
    if nodes: res_list.append(nodes)
    if ppn: res_list.append(ppn)
    if mem: res_list.append(mem)
    if walltime: res_list.append(walltime)
    cmd_args.extend(["-l", ",".join(res_list)])
if depend: cmd_args.extend(["-W", depend])
if extras: cmd_args.append(extras.strip())
# Append the jobscript
cmd_args.append(jobscript)

try:
    res = subprocess.run(cmd_args, check=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

res = res.stdout.decode()
print(res.strip())
