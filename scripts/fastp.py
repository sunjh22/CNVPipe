__author__ = "SUN Jiahong"
__copyright__ = "Copyright 2022, SUN Jiahong"
__email__ = "jiahong.sun@connect.polyu.hk"
__license__ = "MIT"

from snakemake.shell import shell
shell.executable("bash")

samples = snakemake.input[0]
controls = snakemake.input[1]
reference = snakemake.output[0]
cns = snakemake.output[1]
cnr = snakemake.output[2]

outdir = snakemake.params.get('outdir', 'temp/cnvkit')
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

assert len(snakemake.input.control) >= 3, "At least 3 control samples are required."

shell(
    "(cnvkit.py batch {snakemake.input.sample} -n {snakemake.input.control} "
    "-m wgs -f {snakemake.params.ref} --access {snakemake.params.access} "
	"--target-avg-size {snakemake.params.bin_size} -p {snakemake.threads} "
	"--annotate {snakemake.params.refflat} --drop-low-coverage "
    "--output-reference {snakemake.output.reference} -d {outdir}) {log}"
)