# Develop a new CNV calling pipeline - CNVPipe

Objectives of CNVPipe

1. easy to use, high sensitivity and low FDR for CNV calling
2. integrate comprehensive plotting tools for visualization
3. integrate CNV annotation tools: pathogenicity prediction

## Run CNVPipe snakemake

    snakemake -n --directory analysis/ all_cnvfilter #--rerun-triggers mtime
    snakemake -n --directory ~/data3/project/CNVPipe/analysis/
    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data/biosoft/conda-env-cnvpipe --cores 100 --directory analysis/ all_cnvfilter
    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data/biosoft/conda-env-cnvpipe --cores 50 --directory ~/data3/project/CNVPipe/analysis/ all_cnvfilter
    snakemake --directory analysis/ all_cnvfilter --rulegraph | dot -Tpdf > dag.pdf

CNVPipe-token: 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

    git clone https://github.com/sunjh22/CNVPipe.git

The most important questions:
1. how to merge callset from different tools
2. how to simulate appropriate data for benchmarking
3. how to find real WGS dataset and their CNV set for benchmarking

## 01. Simulate data

For developing pipeline, we used 1X data with 6 samples and 6 controls.
For evaluating pipeline, we used 1X, 10X and 30X data with 6 samples and 6 controls.

### 01.1 SCNVSim + ART

`SCNVSim` is a tool to simulate somatic SVs and CNVs in tumors, here we use it to simulate CNVs.
It has two steps: 1. simulate germline SNV and INDELs to get a normal genome, which requires 
reference genome and its length file; `-s` assigns SNV rate, `-r` assigns INDEL rate; 2. simulate 
somatic SVs and CNVs based on the first step genome and an additional repeatmask file; `-s` assigns 
structure variation rate, `-l` assigns copy neutral LOH rate, `-c` assigns number of SVs simulated, 
`-g` assigns mean copy number segment size (default 1M), `-m` assigns minimal deletion size, `-x` 
assigns maximum size of tandem duplication, `-a` prefix.

In this simulation, we need germline SNPs for later SNP calling and BAF correction, we also set
SV rate to be relatively low to make sure only CNV happens in the genome, so we can fully test the
performance of CNV calling tools.

`normal` genome has SNP rate 0.001, while `normal1` genome has SNP rate 0.01.
`tumor1-6` in `10X/` corresponds to `normal`, while `tumor7-12` in `10X/` corresponds to `normal1`.
In the later reads simulation, we used `tumor7-12` as reference to simulate 1X, 10X and 30X data.

    curl -OL https://sourceforge.net/projects/scnvsim/files/latest/download
    unzip download
    cd data/
    time java -jar ~/data/biosoft/scnvsim_1.3.1/normgenomsim_1.3.1.jar -v ~/data/refs/hg38/analysisSet/hg38.analysisSet.size.rmMT.txt -n ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -o ./simulation/simuGenome -s 0.01 -a normal1 1>1.log 2>2.log &
    time java -jar ~/data/biosoft/scnvsim_1.3.1/tumorgenomsim_1.3.1.jar -v ~/data/refs/hg38/analysisSet/hg38.analysisSet.size.rmMT.txt -n ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -k ~/data/refs/hg38/bundle/SCNVSim/repeatmask.txt -i ./simulation/simuGenome/normal1_snvindelsim.vcf -o simulation/simuGenome -s 0.001 -l 0.001 -c 150 -m 500 -x 10000 -a tumor7 1>1.log 2>2.log &
    $ SCNVSim Version: 1.3.1

#### Use ART to simulate reads

`ART` was used to simulate paired-end sequencing reads. `-ss` model Illumina sequencing system, e.g
HSXn - HiSeqX PCR free (150bp), HS25 - HiSeq 2500 (125bp, 150bp), HSXt - HiSeqX TruSeq (150bp); `-f`
assigns coverage, `-l` assigns length, `-m` assigns mean size of DNA fragments in library; `-na` do 
not output ALN alignment file, `-s` standard deviation of DNA fragment size.

In three different depth data simulation, we all use `tumor7-12` genome as reference.
`sample1-6` are 1X depth

    parallel --dry-run "~/data/biosoft/art_bin_MountRainier/art_illumina -ss HSXn -f 1 -l 150 -m 500 -na -p -s 10 -i simulation/simuGenome/normal1_snvindelsim.fasta -o simulation/1X/control{}_ 1>1.log 2>&1" ::: 1 2 3 4 5 6 &
    time ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HSXn -f 1 -l 150 -m 500 -na -p -s 10 -i simulation/simuGenome/tumor7_svcnvsim_clone_0.fasta -o simulation/1X/sample1_ 1>2.log 2>&1 &
    parallel -j 6 -k gzip -q -1 simulation/1X/control{}_1.fq ::: 1 2 3 4 5 6 &
    parallel -j 6 -k gzip -q -1 simulation/1X/sample{}_1.fq ::: 1 2 3 4 5 6 &
    ART Version: v2.5.8

`sample7-12` are 10X depth.

    time ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HSXn -f 10 -l 150 -m 500 -na -p -s 10 -i simulation/simuGenome/normal1_snvindelsim.fasta -o simulation/10X/control7_ 1>1.log 2>2.log &
    time ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HSXn -f 10 -l 150 -m 500 -na -p -s 10 -i simulation/simuGenome/tumor7_svcnvsim_clone_0.fasta -o simulation/10X/sample7_ >1.log 2>2.log &
    parallel -j 6 -k gzip -q -1 simulation/10X/control{}_1.fq ::: 7 8 9 10 11 12 &
    parallel -j 6 -k gzip -q -1 simulation/10X/sample{}_1.fq ::: 7 8 9 10 11 12 &

`sample13-18` are 30X depth

    parallel --dry-run "~/data/biosoft/art_bin_MountRainier/art_illumina -ss HSXn -f 30 -l 150 -m 500 -na -p -s 10 -i simulation/simuGenome/normal1_snvindelsim.fasta -o simulation/30X/control{}_ 1>1.log 2>&1" ::: 13 14 15 16 17 18 &
    time ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HSXn -f 30 -l 150 -m 500 -na -p -s 10 -i simulation/simuGenome/tumor7_svcnvsim_clone_0.fasta -o simulation/30X/sample13_ >1.log 2>2.log &
    parallel -j 6 -k gzip -q -1 simulation/30X/control{}_1.fq ::: 13 14 15 16 17 18 &
    parallel -j 6 -k gzip -q -1 simulation/30X/sample{}_1.fq ::: 13 14 15 16 17 18 &

### 01.2 CNV-Sim + ART

Create a isolated Python2.7 environment for running CNV-Sim, installing all dependencies.

    cd ~/data3/biosoft
    git clone https://github.com/NabaviLab/CNV-Sim.git
    mamba create -n cnvsim python=2.7
    mamba activate cnvsim
    mamba install -c bioconda pysam
    mamba install -c bioconda biopython
    mamba install -c bioconda bedtools
    wget -O ART/art.tgz https://www.niehs.nih.gov/research/resources/assets/docs/artbingreatsmokymountains04.17.16linux64.tgz
    tar -xvzf ART/art.tgz -C ART/
    mkdir -p CNV-Sim/ART
    ln -s /home/jhsun/data3/biosoft/ART/art_bin_GreatSmokyMountains/art_illumina CNV-Sim/ART/
    mkdir ART

Simulate CNV and reads. This command must be ran under `CNV-Sim` directory.

    ./cnv-sim.py -o ~/data3/project/CNVPipe/simulation2/ -l 150 --coverage 1 -g 150 -r_min 5000 -r_max 5000000 -cn_min 1 -cn_max 5 genome ~/data3/refs/hg38/anal

I think the design of CNV-Sim still cannot fulfill my requirement about whole-genome haplotype-level
CNV simulation, so I decided to write a CNV simulator by myself. In this way, I can fully control
and understand every step of CNVPipe, thus give me more confidence about it.

### 01.3 CNV-simulator (Python, by Jiahong)

We wrote a CNV simulator by Python, which can simulate CNVs in whole-genome accessible regions and 
generate NGS reads by ART. The default CNV number is 200, average size of CNV is 200,000, 
amplification rate is 0.5, read length is 150. We use `art_illumina` `HSXn` (as well HiSeqX PCR 
free (150bp)) mode to simulate paired-end reads. Users could set coverage for read depth. 
The tool is in `~/data3/project/CNV-simulator`.

To simulate experimental samples with CNVs. This will simulate sample genomes and NGS reads,
concatenate and compress reads automatically.

Sample1-6 are 1X; sample7-12 are 10X; sample13-18 are 30X; sample19-24 are 0.1X; 
sample25-30 are 0.5X. sample31-36 are 5X.

    cd ~/data3/project/CNV-simulator
    ./cnv_simulator.py -o ~/data3/project/CNVPipe/simulation-CNVSimulator/1X -a sample1 -c 0.5 ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa ~/data3/refs/hg38/bundle/CNVKit/access-excludes.hg38.analysisSet.bed
    parallel --dry-run "cnv_simulator -o ~/data3/project/CNVPipe/simulation-CNVSimulator/1X -a sample{} -c 0.5 ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa ~/data3/refs/hg38/bundle/CNVKit/access-excludes.hg38.analysisSet.bed 1>cnv-simu.{}.log 2>&1" ::: 1 2 3 4 5 6 &

    parallel "cnv_simulator -o ~/data3/project/CNVPipe/simulation-CNVSimulator/0.1X -a sample{} -c 0.05 -e 1000000 -b 500000 -B 3000000 ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa ~/data3/refs/hg38/bundle/CNVKit/access-excludes.hg38.analysisSet.bed 1>cnv-simu.{}.log 2>&1" ::: 19 20 21 22 23 24

    parallel "cnv_simulator -o ~/data3/project/CNVPipe/simulation-CNVSimulator/5X -a sample{} -c 2.5 ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa ~/data3/refs/hg38/bundle/CNVKit/access-excludes.hg38.analysisSet.bed 1>cnv-simu.{}.log 2>&1" ::: 31 32 33 34 35 36 &

To simulate control samples, directly use reference genome.

    cd ~/data3/project/CNVPipe
    parallel --dry-run "art_illumina -ss HSXn -f 1 -l 150 -m 800 -na -p -s 10 -i ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa -o simulation-CNVSimulator/1X/control{}_ 1>1.log 2>&1" ::: 1 2 3 4 5 6 &
    parallel -j 6 -k gzip -q -1 simulation-CNVSimulator/1X/control{}_1.fq ::: 1 2 3 4 5 6 &
    art_illumina version 2.5.1
    GNU parallel 20221222

Run CNVPipe for 1X samples

    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data3/project/CNVPipe/envs/conda-env-cnvpipe --cores 150 --directory ~/data3/project/CNVPipe/analysis-CNVSimulator/ --rerun-triggers mtime

Run CNVPipe for 10X samples, let CNVPipe to determine a resolution.

    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data3/project/CNVPipe/envs/conda-env-cnvpipe --cores 150 --directory ~/data3/project/CNVPipe/analysis-CNVSimulator/ autobinBydepth

Run CNVPipe for 30X samples, let CNVPipe to determine a resolution, remember to remove the former
`log/autobin/binsize.txt` file.

Run CNVPipe for 0.1X samples, let CNVPipe to determine a resolution, remember to remove the former
`log/autobin/binsize.txt` file. Delly and Smoove both throw errors for too low coverage.

Run CNVPipe for 0.5X samples, let CNVPipe to determine a resolution, remember to remove the former
`log/autobin/binsize.txt` file.

## 02. Build CNVPipe step by step

We use Snakemake as beckend to build CNVPipe, separate the whole pipeline into several snakefiles,
each of them performs specific functions, like fastp processing, bwa alignment, CNV calling, SNP
calling, merging etc.

### 02.1 Main Snakefile and config file

Main file: include `.smk` file step by step, testing their availability during development.

Config file: 1. data: reference genome and resource bundle for different tools to run; 2. settings: 
tools available for specific tasks; 3. params: parameters for running the tool, basically are 
setting threads for different tools.

### 02.2 common.smk

Import config file, get fastq files (depending on se or pe), get patient sample and control sample 
names and set them to be global parameter, valid filename and filepath. CNVPipe interface, version.

### 02.3 pre-processing.smk

Build fasta index, bwa index (if not available) and dictionary for reference genome.

Includes another two snakefiles:
- `fastp.smk` for cleaning reads;
- `bwamem.smk` for mapping reads using bwa-mem, also includes marking duplicates and BQSR

### 02.4 CNV calling

We selected five tools for CNV calling, they are: CNVKit, CNVpytor, cn.MOPS, Smoove and Delly.
They will be run in parallel in CNVPipe. We will finally merge the results from them and assign
various scores to the final set to help refine important CNVs.

One of the most important parameter for CNV calling might be the resolution as well as the bin size.
Users could specifically define the resolution in config file, or alternatively, ask CNVPipe to
estimate one, but that will divide CNVPipe into two stages: one is for reads alignment and
resolution estimation, the other is CNV calling and downstream analysis.

If we want to estimate an optimal calling resolution (bin size) for WGS dataset by CNVPipe, 
we can run the following command. Otherwise, CNVPipe will use the resolution in config file.

    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data3/project/CNVPipe/envs/conda-env-cnvpipe --directory ~/data3/project/CNVPipe/analysis-bqsrTest/ autobinBydepth

We need to be clear about some features of these methods:
1. require normal samples? if yes, how many?
2. call germline CNV sample by sample or in batch?
3. exclude low-mappability region? require reference file?

#### 02.4.1 CNVKit

Require normal samples (no limit on the number), reference genome and exclude low-map region 
(pre-defined). Call germline CNV in batch.

Include four rules:
- `autobin` for estimating optimal bin resolution;
- `batch` for calling CNV for WGS data in batch mode;
- `segmetrics` for segmenting CNVs based on confidence interval;
- `call` for transforming log2 depth ratio into integer copy numbers.

And a simple shell command using `awk` to extract informative columns into another file. 
Columns include genomeic coordinates, integer CN, log2 ratio, depth, probe and weight.

The refFlat and access file should be provided for annotation and discarding reads in low-mappbility
 region.

*Principle of CNVKit*

The first step is to get accessible regions. Find low mappability track: 1. get from rqfu, in 
`~/data/refs/hg38/low-mappability-track/hg38.badRegions.bed`, which includes all centromeres, 
telomeres and heterochromatin regions; 2. download from 10x genomics - black list of SV 
<http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed> to
`~/data/refs/hg38/low-mappability-track/sv_blacklist.bed`, re-format it using
`grep -E -w "chr([0-9]+|[XY])" low-mappability-track/sv_blacklist.bed | sort -Vk 1 -k 2,3n > low-mappability-track/hg38.sv_blacklist.bed`.

Then use `cnvkit.py access hg38.analysisSet.fa -x hg38.badRegions.bed -x hg38.sv_blacklist.bed -o access-excludes.hg38.bed`
to get accessible regions, total length of these regions is 2,875,824,894. `cnvkit.py access` 
command will automatically find sequence of N's in genome and filter regions with long N's. 

Then we get target bin file by 
`cnvkit.py target refs/access-excludes.hg38.bed --annotate ~/data/refs/hg38/bundle/CNVKit/refFlat.txt --avg-size 20000 --split -o refs/target.bed`,
which will then be used to calculate coverage of samples in these bins.

#### 02.4.2 cnvpytor

Do not require normal sample and reference genome (a pre-defined GC file is in the package) and do 
not exclude low-map region. This method solely depends on comparing read depth between adjacent 
region (?) and regions with similar GC content. Call germline CNV sample by sample.

Include one rule with four commands:
- `rd` for extracting reads from sample bam file;
- `his` for calculating reads depth at certain resolution;
- `partition` for segmenting the bins;
- `call` for calling CNV.

The results include "CNVtype, CNVregion, CNVsize, CNVlevel, eval1, eval2, eval3, eval4, q0, pN, dG"
- `CNVlevel` - normalized read depth value, rough integer CN value should be `round(2*CNVlevel)`;
- `eval1` - e-value (p-value multiplied by genome size divided by bin size) calculated using t-test 
            statistics between RD statistics in the region and global - set to [0 - 0.00001];
- `eval2` - e-value from the probability of RD values within the region to be in the tails of a 
            gaussian distribution of binned RD - set to [0 - 0.00001];
- `q0` - fraction of reads mapped with q0 quality in call region - not useful for filtering since 
         most reads have good quality;
- `pN` - fraction of reference genome gaps (Ns) in call region - set to [0 - 0.5];
- `dG` - distance from closest large (>100bp) gap in reference genome - set to [>100kb];

Get read depth signal in 100bp level

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -chrom $(seq -f 'chr%g' 1 22) chrX chrY -rd sample1.bam &

Assign bin size, calculate read depth histogram

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -his 20000 &

Segment by mean-shift method

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -partition 20000 &

Call CNV

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -call 20000 > ../read-depth/cnvnator/sample1.2k.tsv &

Plot

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -plot rd 20000 -o ../read-depth/cnvnator/sample1.png &

Finally we use a python script to extract evalue, pN and dG into `info` column

`cnvpytor` accepts SNP data and call CNV based on both RD and BAF, but we want to separate them,
just use its RD function to call CNV, and filter by BAF after merging CNV by other tools.

*Principle of CNVpytor*

Firstly, extract RD signal in 100-bp intervals (pysam), get bin-level RD signal, perform GC biase
correction (100-bp GC content is pre-calculated), then merge and partition by **mean-shift**, CNV
call follows, all requires statistical test of RD between CNV segments and whole-genome regions.

Filter out HET SNPs that are inaccessible to short read, defined by strict mask from 1000G project,
BAF likelihood of bin is calculated by multiplying the likelihood of single SNPs in the bin (which
is modelled by a beta distribution ?).

Annotating CNV with genes (inside, covering or intersecting left/right breakpoints)

Merge calls over multiple regions: >50% reciprocal overlap

#### 02.4.3 cn.MOPS

Require at least 6 samples, do not require reference genome and do not exclude low-map 
region. Call CNV in batch. This method call CNV by comparing read depth between samples.

cn.mops is not applicable to very-low-coverage data, since Windows should contain 50-100 reads each.
Three steps in cn.mops:
1. calculate read counts in windows (bins) from bam file
2. call CNV
3. calculate integer copy numbers

Run cn.MOPS

    Rscript scripts/cnmops_wgs.R analysis/res/mops 20000 10 analysis/mapped/sample*.bam 2> analysis/logs/mops/call.log

Some unexpected error happened for test data, at least 6 samples are required for running mops.
I suspect it may due to the extremely low coverage of test data. The next step should be
increasing the sample number and coverage and test again. But now I have to stop at here because
Dr. Jin ask me to focus more on scCNV project. I will catch up again if I have time.
(问心无愧即可). Catch up here at Wed Oct 12 15:47:48 CST 2022.

#### 02.4.4 Lumpy and Smoove

`Lumpy` is written in C.

Five steps to call CNV by Lumpy. 1. get discordant reads (insert size over expected size or
mate pair mapped to different chromosomes); 2. get split reads by Lumpy script; 3. call CNV by
lumpyexpress; 4. genotype CNV by svtypes, a quality score was given (in QUAL filed of vcf file),
but how to select CNV based on score need further consideration; 5. extract CNV.

    samtools view -b -F 1294 analysis/mapped/sample1.bam | samtools sort -@ 8 -o analysis/temp/lumpy/sample1.discordant.bam -
    samtools view -h analysis/mapped/sample1.bam | extractSplitReads_BwaMem -i stdin | samtools sort -@ 8 -o analysis/temp/lumpy/sample1.split.bam -
    lumpyexpress -B analysis/mapped/sample1.bam -S analysis/temp/lumpy/sample1.split.bam -D analysis/temp/lumpy/sample1.discordant.bam -o analysis/temp/lumpy/sample1.vcf
    svtyper -i analysis/temp/lumpy/sample1.vcf -B analysis/mapped/sample1.bam -l analysis/temp/lumpy/sample1.bam.json > analysis/temp/lumpy/sample1.gt.vcf
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL\n' analysis/temp/lumpy/sample1.gt.vcf | egrep 'DUP|DEL' > analysis/res/lumpy/sample1.bed

A coverage calculation script was provided `scripts/get_coverages.py` for later reference.

A low mappbility track file provided by Heng Li is suggested, I downloaded and extracted
the file at Mon Oct 17 08:59:27 +08 2022, the file is in 
`~/refs/hg38/low-mappability-track/HengLi-lowMapTrach/hg38.exclude.hengli.bed`.
Overall length of low-map region is 228,061,378.


*Smoove* is a Lumpy wrapper, which is easier to use and faster.

A bad region file was downloaded according to the author of smoove. Overall length is 119,556,880,
which is almost half shorter than Heng Li's. But when we remove non-cannonical chromosomes, its 
length became 3,042,685, which is definitely not correct. 
The file is in `low-mappability-track/smoove/hg38.exclude.bed`.

    wget https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed

Our own bad region file is in `low-mappability-track/hg38.badRegions.bed`, which includes 
centromeres, telomeres and heterochromatins, the length is 199,765,358.

Another low-map track is downloaded from 10x genomics in 
`low-mappability-track/10x/hg38.sv_blacklist.bed`, its overall length is 212,765,070.
I think this one might be the *best* one to use.

Test `smoove` calling

    smoove call --outdir analysis/temp/smoove/ --exclude ~/data/refs/hg38/exclude.cnvnator_100bp.GRCh38.20170403.bed --name sample4 --fasta ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -p 1 --genotype analysis/mapped/sample4.bam

After comparing the results from Lumpy and Smoove, we found no big difference, so we will replace 
lumpy with smoove in our CNVPipe.

Note: `Lumpy` might produce multiple calls for the same region, which indicates a noisy region or 
complex SV, we will skip such regions when merging.

    python ~/data/project/CNVPipe/scripts/smooveFilter.py analysis/temp/smoove/sample17.bed analysis/res/smoove/sample17.bed
    for x in sample{1..18}; do python ~/data/project/CNVPipe/scripts/smooveFilter.py analysis/temp/smoove/${x}.bed analysis/res/smoove/${x}.bed; done

#### 02.4.5 Delly

Delly is written in C++. It is primarily designed for calling structural variant, both germline and 
somatic (cancers). It also supports CNV calling now, I suspect in this mode read-depth information 
is also used.

Mappbility files are mandantory for Delly to call CNV, we downloaded map files
from [here]<https://gear.embl.de/data/delly/> at Sat Oct 15 17:08:30 +08 2022. All three files are 
required.

Delly by default divide genome into 10kb-mappable bins, but we can set window size by `-i`
parameter. `-u` is for segmentation, `-l` is for using delly SV calling to refine breakpoints.

So we try to call SV first, then use its output as input for CNV calling. 
`delly call` only use paired-end and split-read information, read-depth is not used.

    delly call -g ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -x ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed -o analysis/temp/delly/sample1.sv.bcf analysis/mapped/sample1.bam

Use `delly cnv` to call CNV. It is better to call sv first and use it as input for refining the CNV
calling. However, Delly SV calling is very time-consuming, for a sample with 46G bam file, it takes
around 24 hours to finish. Therefore, we decide to skip this step and directly call CNV from bam
file.

    delly cnv -u -i 20000 -g ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -l analysis/temp/delly/sample1.sv.bcf -m ~/data/refs/hg38/bundle/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz -c analysis/temp/delly/sample1.cov.gz -o analysis/temp/delly/sample1.cnv.bcf analysis/mapped/sample1.bam

Same as `smoove`, we extract specific columns from `delly cnv` output, merge its results with that 
of other tools, and genotype with `duphold` in single run.

`delly classify` will select CNVs with 'PASS' (filter low-quality CNVs), bad thing is that it
transforms all quality score into 10000, so we decided to filter CNVs by grep and keep
the original quality score. One thing need to be noticed is that in some cases, value in `FILTER`
field is not matched with `FT` tag in `format` field. One could be PASS while the other one
could be LowQual. We just use the `FILTER`=="PASS" to select all good-quality CNVs. (strict mode)

Some R scripts were provided for plotting CNV, which could be a reference for later usage.

    Rscript ~/data/biosoft/delly/R/rd.R analysis/temp/delly/sample1.cov.gz analysis/temp/delly/sample1.segment.bed

    Version: Delly 1.1.5

#### 02.4.5 Control-FREEC (deprecated)

`Control-FREEC` is written in C++.

This tool can be run in two modes: 1. for WGS data without normal reference, a directory with fasta
of single chromosomes or a GC-content profile need to be provided; 2. for WES data, a normal 
reference must be provided. For mode 1, a `GEM` mappability file could be provided to remove bad
genomic regions. In my testing, control parameter is not working even if you provide it, do not know 
the reason.

Generate GC content profile for specific size of window. In this case, `chrFiles` should be set to
the directory with fasta sequences of single chromosomes, and `window` set to specific size. We can 
use very-low-depth data to do this.

    ~/data/biosoft/FREEC/src/freec -conf analysis/temp/freec/config_WGS.txt -sample analysis/mapped/sample1.bam

Run freec with specified GC-content profile and mappability file for high-coverage data

    ~/data/biosoft/FREEC/src/freec -conf analysis/temp/freec/config_WGS.txt -sample analysis/mapped/sample1.bam

After checking the results, we found that even a mappability file is provided, some telomere regions 
are called to be CNV, which either results from the problem of mappability file or from the tool
itself. At least, this result makes its all calling results not reliable.

*Principle of Control-FREEC*

Control-FREEC takes as an input aligned reads, then constructs and normalizes the copy number 
profile (by GC profile or normal read coverage), constructs the B-allele frequency (BAF) profile, 
segments both profiles (Lasso-based algorithm), ascribes the genotype status to each segment using 
both copy number and allelic frequency information, then annotates genomic alterations.

### 02.5 SNP calling

#### 02.5.1 freebayes

Call SNPs sample by sample (the problem is the accuracy). Then filter based on quality score and 
read depth, this is a critical step. We should set different threshold for different depth data. 
The accuracy of SNP calling will affect the correction of CNV in later step.

    freebayes -f ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa mapped/sample1.bam > snps/freebayes/sample1.raw.vcf
    bcftools filter -O v -o snps/freebayes/sample1.filtered.vcf -s LOWQUAL -e 'QUAL<10 || FMT/DP <5' --SnpGap 5 --set-GTs . snps/freebayes/sample1.raw.vcf
    bcftools view -v snps snps/freebayes/sample1.filtered.vcf > snps/freebayes/sample1.snp.vcf

#### 02.5.2 GATK

Download GATK bundle resource file from 
[here]<https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/>
at Thu Oct 27 11:29:27 +08 2022. I suspect bundle files in our server might have broken or I misuse
some of files in workflow at first. After downloading and checking with md5sum, I found they are 
totally identical.

When runing `variantRecalibrator`, an error happened: "Bad input: Found annotations with zero 
variance". They must be excluded before proceeding. This might be due to the very-low-coverage of
our testing data.
Another error is "Positive training model failed to converge.  One or more annotations (usually MQ) 
may have insufficient variance", this is caused by too less number of SNPs identified in GATK 
`HaplotypeCaller`. We need to consider re-simulate sequencing data with higher SNP rate.

In CNVPipe, we actually single-sample mode to call SNPs, which is not as accurate as multi-sample
mode, as the later one could use the SNP allele information from all samples, but single-sample mode
is easy to implement and runs faster. We might optimize this step if necessary.

    time gatk --java-options '-Xmx30G' HaplotypeCaller -R /data/jinwf/jhsun/refs/hg38/analysisSet/hg38.analysisSet.fa -I mapped/sample18.bam -O snps/gatk/sample18.raw_variants.vcf.gz --dbsnp /data/jinwf/jhsun/refs/hg38/bundle/gatk/dbsnp_146.hg38.vcf.gz 1>1.log 2>&1 &

    gatk --java-options '-Xmx30G' VariantRecalibrator -R /data/jinwf/jhsun/refs/hg38/analysisSet/hg38.analysisSet.fa -V snps/gatk/sample18.raw_variants.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/jinwf/jhsun/refs/hg38/bundle/gatk/hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=false,prior=12.0 /data/jinwf/jhsun/refs/hg38/bundle/gatk/1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 /data/jinwf/jhsun/refs/hg38/bundle/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz -an QD -an ReadPosRankSum -an FS -an SOR -mode SNP --max-gaussians 4 -O temp/gatk/sample18.recal --tranches-file temp/gatk/sample18.tranches --rscript-file temp/gatk/sample18.plots.R 1>1.log 2>&1

    gatk --java-options '-Xmx30G' ApplyVQSR -R /data/jinwf/jhsun/refs/hg38/analysisSet/hg38.analysisSet.fa -V snps/gatk/sample18.raw_variants.vcf.gz -O snps/gatk/sample18.vqsr.vcf.gz --truth-sensitivity-filter-level 99.9 --create-output-variant-index true --tranches-file temp/gatk/sample18.tranches --recal-file temp/gatk/sample18.recal -mode SNP

    gatk v4.2.2.0
    HTSJDK Version: 2.24.1
    Picard Version: 2.25.4

#### 02.5.3 samtools + bcftools (not implemented)

Almost the same as freebayes.

### 02.6 Merge the results from all CNV calling tools

#### Strategy 1

We tried to merge CNV results from all tools including cnvkit, delly, cnvpytor, smoove and mops.
During the process, we merge two CNVs iteratively, if there is overlap between two CNV region, 
their overlapped length will be accumulated, and the later CNV will dispear (only the first CNV 
will be kept). After the whole merging process, there should not exist any overlap CNVs. The
implementation is in `scripts/mergeCNV.py`. The priority to keep CNVs: cnvkit, delly, smoove, mops,
cnvpytor.

    cd analysis/
    python ~/data/project/CNVPipe/scripts/mergeCNV.py res/cnvkit/sample18.bed res/delly/sample18.bed res/smoove/sample18.bed res/mops/sample18.bed res/cnvpytor/sample18.bed ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed res/merge/sample18.bed >merge.log

    for x in sample{1..18}; do python ~/data/project/CNVPipe/scripts/mergeCNV.py res/cnvkit/${x}.bed res/delly/${x}.bed res/smoove/${x}.bed res/mops/${x}.bed res/cnvpytor/${x}.bed ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed res/merge/${x}.bed > logs/merge/mergeCNV.log; done
    (grep chromosome res/merge/sample18.bed; grep -v chromosome res/merge/sample18.bed | sort -k 7nr -k 6nr -Vk 1,1) | less

I created a test dataset in `~/data3/project/CNVPipe/analysis/res/merge/test/` to test this script. 
Testing is actually very important, I fixed some bugs after testing, which may not be identified if 
no testing because the script runs OK for real data.

    cd analysis/res/merge/test/
    python ~/data/project/CNVPipe/scripts/mergeCNV.py test1.bed test2.bed test3.bed test4.bed test5.bed test6.bed ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed test.merge.bed

In `analysis/res/merge2/`, we only merge cnvkit, delly mops and cnvpytor. The script is in 
`mergeCNV2.py`.

    for x in sample{1..18}; do python ~/data/project/CNVPipe/scripts/mergeCNV2.py res/cnvkit/${x}.bed res/delly/${x}.bed res/mops/${x}.bed res/cnvpytor/${x}.bed ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed res/merge2/${x}.bed > logs/merge/mergeCNV.2.log; done

#### Strategy 2

We next tried another merging strategy: if two CNVs have overlap, instead of taking only former one, 
this time we really merge two CNVs - extending the breakpoints. The script is in `mergeCNV3.py`. 
Notice that this is based on mergeCNV2.py, which means only four tools' results are merged.

Test first

    python ~/data/project/CNVPipe/scripts/mergeCNV3.py test1.bed test2.bed test3.bed test4.bed ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed test.merge.bed

Run for simulation data

    for x in sample{1..18}; do python ~/data/project/CNVPipe/scripts/mergeCNV3.py res/cnvkit/${x}.bed res/delly/${x}.bed res/mops/${x}.bed res/cnvpytor/${x}.bed ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed res/merge3/${x}.bed > logs/merge/mergeCNV.3.log; done

### 02.7 Score by duphold results

We sse `duphold` to do genotype for merged CNV set. Three additional fields will be added into 
`FORMAT` field per sample:
1. DHFC - fold-change for the variant depth relative to the rest of the chromosome where variant was
2. DHFFC - fold-change for the variant depth relative to Flanking regions
3. DHBFC - fold-change for the variant depth relative to bins in the genome with similar GC-content

*General usage of duphold*

`duphold` could be directly downloaded as an exectuable binary file, or used after installing 
`smoove`.

    duphold -t 4 -v analysis/temp/delly/sample4.filtered.bcf -b analysis/mapped/sample4.bam -f ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -o analysis/temp/delly/sample4.duphold.vcf
    smoove genotype -d -x -p 1 --name sample4 --outdir analysis/temp/smoove-genotype/ --fasta ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa --vcf analysis/temp/smoove/sample4-smoove.genotyped.vcf.gz analysis/mapped/sample4.bam
    Version: Smoove 0.2.8

In `INFO` field, `MP` means fraction of mappable positions; `GCF` means GC fraction (added by 
duphold); these will not be used for CNV filtering.

The authors of `duphold` suggest use DHFFC<0.7 to filter deletion and DHBFC>1.3 to filter duplication

    bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] > 1.3)' analysis/temp/smoove-genotype/sample4-smoove.genotyped.vcf.gz

Use `bcftools` to extract informative columns

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL[\t%DHFC\t%DHFFC\t%DHBFC]\n' analysis/temp/smoove-genotype/sample4-smoove.genotyped.vcf.gz

Follow author's suggestion to filter CNV based on DHFFC and DHBFC by `awk`

    awk '$4=="DEL" && $7<0.7 {print $0} $4=="DUP" && $8>1.3 {print $0}' analysis/res/smoove/sample1.bed | wcl

`duphold` also support SNP vcf file for filtering, just like `CNVfilteR`.

*Implementation of duphold in CNVPipe*

After merging, we got CNV set in bed file, but `duphold` only accept vcf format file, thus we 
transform the bed file into vcf file. Then we use it to calculate adjacent read depth for CNV 
region, and design a strategy to score CNV regions. Transforming script is in `scripts/bed2vcf.py`.

    cd analysis/res/merge
    python ~/data/project/CNVPipe/scripts/bed2vcf.py sample18.bed ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa.fai sample18.vcf
    time duphold -t 10 -v sample18.vcf -b ../../mapped/sample18.bam -f ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -o sample18.duphold.vcf &

    for x in sample{1..18}; do python ~/data/project/CNVPipe/scripts/bed2vcf.py ${x}.bed ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa.fai ${x}.vcf; done
    parallel --dry-run duphold -t 10 -v sample{}.vcf -b ../../mapped/sample{}.bam -f ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -o sample{}.duphold.vcf ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 &

Then use `bcftools` to extract columns out and score based on DHFFC and DHBFC.
The authors suggest to use DHFFC<0.7 to filter deletion and DHBFC>1.3 to filter duplication.
We score deletion with DHFFC<0.5 or duplication with DHBFC>1.5 100, and score deletion with 
DHFFC<0.7 or duplication with DHBFC>1.3 90. In this way, we can find highly confident CNVs.

`duphold` 0.2.1 keep throwing an error: "fatal.nim(49)            sysFatal
Error: unhandled exception: index -1 not in 0 .. 159345972 [IndexDefect]".
Version 0.2.3 works well.

    bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%CN\t%AS\t%GS\t%DHFFC\t%DHBFC]\t%INFO/TN\t%INFO/SAMPLE\n' sample18.duphold.vcf > sample18.duphold.bed
    python ~/data/project/CNVPipe/scripts/scoreDuphold.py sample18.duphold.bed sample18.duphold.score.bed
    (grep chromosome sample18.duphold.score.bed; grep -v chromosome sample18.duphold.score.bed | sort -k 7nr -k 8nr -k 6nr -Vk 1,1) | ct | less -SN

    cd analysis/res/merge
    for x in sample{1..18}; do bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%CN\t%AS\t%GS\t%DHFFC\t%DHBFC]\t%INFO/TN\t%INFO/SAMPLE\n' ${x}.duphold.vcf > ${x}.duphold.bed; done
    for x in sample{1..18}; do python ~/data/project/CNVPipe/scripts/scoreDuphold.py ${x}.duphold.bed ${x}.duphold.score.bed; done

### 02.8 Score by CNVfileR

Use `CNVfilteR` to do BAF correction

Run CNVfilteR:
1. load cnv file
2. load snv file
3. filterCNVs
4. plot single cnvs or all cnvs

CNVfilteR applies a scoring model, which is based on fuzzy logic, to filter false positive CNVs. 
Briefly, copy number deletion could be filtered if there is over 30% (default) heterozygous snvs in 
CNV region; For copy number duplication, scoring model is used: the variants with allele frequency 
close to 50% will be given positive value ([0,1]), and variants with allele frequency close to 33.3% 
or 66.6% will be given negative value ([-1,0]), a total score is calculated for CNV region by 
summing up all the score of snvs in that CNV region; finally, positive value indicates FP CNV and 
should be filtered.

This method requires as accurate SNP calling as possible. Is 10X data really suitable for this tool? 
We can rank the CNVs by score but not remove them.

The script is in `scripts/cnvFilteR.R`

For 10X data, too few SNPs were called by freebayes, so this tool only filter several CNV regions. 
However, for 30X data, even over 1037411 SNPs were called, this tool still only remove very few 
CNV regions. I am not sure it is the problem of `cnvfilter` or it is due to SNP calling, I will try 
using `GATK` to call SNPs from 30X data and check again. Use one sample to test first.

### 02.9 Score by overlap fraction with genomic bad regions and CNVs in normal population

When merging CNV results, we also calculate the overlap fraction of CNV region and low-map region, 
then assign a `goodScore` for each CNV, which is negativly correlated with overlap fraction. 
The higher overlap fraction, the lower score CNV region will get. Low-map region file is actually
the 10x sv blacklist file.

*High-frequency CNVs in normal population*

I downloaded dbVar Common SV set and conflict SV set (between common and pathogenic SV) into
`~/data/refs/hg38/genomic-variation-track/dbVar`. Then I filtered conflict SV entries in
Common SV set to get a common sv set with only non-pathogenic CNVs (or we can call it CNVs in 
normal population), the file is 
`~/data/refs/hg38/genomic-variation-track/dbVar/common_global_normal.bed`.

    cd ~/data3/refs/hg38/genomic-variation-track/
    wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/dbvarhub/hg38/common_global.bed
    wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/dbvarhub/hg38/conflict_pathogenic.bed
    python filterClinicalSV.py common_global.bed conflict_pathogenic.bed common_global_normal.bed

948 of 1043 conflict SVs could be found in curated commom SV list, theoratically all conflict SVs
should be in that list. This is weired.

The sum length of the common global sv list is 406,138,486.

    cut -f 1-3 common_global_normal.bed | CNVSum

### 02.10 Predict CNV pathogenicity

Use `ClassifyCNV` to predict pathogenicity.

`ClassifyCNV` cannot be downloaded from bioconda, and its result directory is weired, so we 
re-configured it to adjust to CNVPipe workflow. Tested by a small list of cnvkit result.

    cd ~/data3/project/CNVPipe/analysis-bqsrTest/res/cnvkit
    head sample1.bed > test2.bed
    python ~/data3/github-repo/CNVPipe/scripts/classifyCNV.py --absPath /home/jhsun/data3/github-repo/CNVPipe/ --infile test.classifycnv.bed --GenomeBuild hg38 --cores 1

### 02.11 Get final CNV list

The final CNV set file is in `res/merge/sample.final.bed`.

### 02.12 visualization and report

Tools pending for test:
1. wally - https://github.com/tobiasrausch/wally


## 03. Evaluate performance on simulation data

### 03.1 Data simulated by SCNVSim + ART

Construct ground truth set: we have two truth set when simulating data - one is CNV, the other is
SV, two types in SV - Deletion and TandemDup should also be considered as CNV, so we need to merge 
two set and get a non-overlapped ground truth set for benchmark. Here, we do not compare precise
copy number but their CNV type - deletion (CN<2) or duplication (CN>2). Besides, simulated CNVs that
are smaller than 5kb were filtered. The script is in `scripts/mergeTruthSet.py`, the test set is in
`~/data3/project/CNVPipe/simulation/simuGenome/test`.

    cd ~/data3/project/CNVPipe/simulation/simuGenome/
    for x in tumor{1..12}; do python ~/data/project/CNVPipe/scripts/mergeTruthSet.py ${x}_CNV_clone_0.bed ${x}_SV_Clone_0.txt ${x}.truth.tmp.bed; done
    for x in tumor{1..12}; do sort -Vk 1 -k 2,3n ${x}.truth.tmp.bed > ${x}.truth.bed; rm ${x}.truth.tmp.bed; done

I created another test dataset in `~/data3/project/CNVPipe/analysis/res/merge/test/` to make sure 
`evalPerform.py` script works appropriately.

    python ~/data/project/CNVPipe/scripts/evalPerform.py test/callSet.bed test/truthSet.bed single

For simulation data, sample13-18 corresponds to truth7-12, which is hard for parallel computing, 
so we copy truth7-12 to truth13-18 in `~/data3/project/CNVPipe/simulation/simuGenome/`.

    python ~/data/project/CNVPipe/scripts/evalPerform.py sample13.duphold.score.bed ~/data3/project/CNVPipe/simulation/simuGenome/tumor7.truth.bed merge
    python ~/data/project/CNVPipe/scripts/evalPerform.py ../cnvkit/sample13.bed ~/data3/project/CNVPipe/simulation/simuGenome/tumor7.truth.bed single

    for y in {1..18}; do echo -e 'sample' ${y}; for x in cnvkit smoove delly mops cnvpytor; do echo -e 'tool' ${x}; python ~/data/project/CNVPipe/scripts/evalPerform.py ../${x}/sample${y}.bed ~/data3/project/CNVPipe/simulation/simuGenome/tumor${y}.truth.bed single; echo -e '\n'; done; echo -e '\n'; done
    for x in {1..18}; do echo -e 'sample' ${x}; python ~/data/project/CNVPipe/scripts/evalPerform.py sample${x}.duphold.score.bed ~/data3/project/CNVPipe/simulation/simuGenome/tumor${x}.truth.bed merge; echo -e '\n'; done

To quickly evaluate performance under different parameters, I re-format `evalPerform.py`, only
output file need to be assigned now.

    python ~/data/project/CNVPipe/scripts/evalPerform.py ~/data3/project/CNVPipe/analysis/evaluation/sensitivity-FDR.v2.txt

Different merging strategies have corresponding evaluation scripts for them.

I re-simulated 30x control data and re-run CNVPipe, start at Sat Dec  3 14:29:53 +08 2022 with
120 cores, 

We used 50k resolution for 1x data, 5k resolution for 10x data and 1k resolution for 30x data.

Note:
1. For 30X-depth data, cnvkit and delly performs really good, with over 0.8 sensitivity and less 
than 0.05 FDR, mops and cnvpytor have moderate sensitivity, smoove follows, and freec perform the 
worse, in every sample, it has 0 sensitivity and 1 FDR, which means it did not identify any true 
CNV, this is abnormal, I need to check this.
2. We decided to remove freec in our workflow, not merging or evaluating it in this part, this 
improved the sensitivity of CNVPipe but its FDR is still higher than other tools.
3. For CNVPipe, we have some threshold to filter CNV regions, however, loose the cutoff will
increase sensitivity and FDR at the same time.
4. smoove shows high FDR

5. goodScore is useless in simulation study, should be set to negative infinity
v1: cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 0; goodScoreThe = 0; dupholdScoreThe = 30
v2: cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 0; goodScoreThe = -20; dupholdScoreThe = 0

6. increase accumScoreThe by 1 (0->1) could greatly reduce FDR
v3: cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 1; goodScoreThe = -20; dupholdScoreThe = 30

7. it looks like filtering based on dupholdScore (>=90) could not improve performance
v4: cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 1; goodScoreThe = -1000; dupholdScoreThe = 90
v5: cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 1; goodScoreThe = -1000; dupholdScoreThe = 0

8. Simply increase overlap proportion threshold will not improve performance
v6: cnvProp1 > 0.45 or cnvProp2 > 0.45; accumScoreThe = 1; goodScoreThe = -1000; dupholdScoreThe = 0
v7: cnvProp1 > 0.6 or cnvProp2 > 0.6; accumScoreThe = 1; goodScoreThe = -1000; dupholdScoreThe = 0

9. When set toolNum=2, both sensitivity and FDR of CNVPipe become comparable with cnvkit and delly.
v8: cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0, toolNum = 2
v9: cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe + goodScoreThe + dupholdScoreThe > 130

10. Even we take all CNVs from merged results, its sensitivity is not significantly better, which
should not be the case
v10: cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0

As beforehand results are not good enought, We next try to merge only cnvkit, delly, mops and 
cnvpytor.
v11: four tools, cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0

    python ~/data/project/CNVPipe/scripts/evalPerform2.py ~/data3/project/CNVPipe/analysis/evaluation/sensitivity-FDR.v11.txt

v12: four tools, cnvProp1 > 0.3 or cnvProp2 > 0.3; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0, toolNum>=2

We next try to merge cnvkit, delly, mops and cnvpytor by a new strategy: merge two CNVs if they
overlap rather than take the former one. And we require the called CNV at least cover 30% of truth
CNV.

v13: four tools, cnvProp1 > 0.3; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0, toolNum>=2

    python ~/data/project/CNVPipe/scripts/evalPerform3.py ~/data3/project/CNVPipe/analysis/evaluation/sensitivity-FDR.v13.txt

v14: four tools, cnvProp1 > 0.45; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0, toolNum>=2

v15: four tools, cnvProp1 > 0.6; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0, toolNum>=2

In v16, we got plausible result: higher sensitivity and much lower FDR for CNVPipe. We might keep
exploring better merging and filtering strategy.

v16: four tools, cnvProp1 > 0.8; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0, toolNum>=2

v17: re-simulate and analyze 1x data (50k), four tools, cnvProp1 > 0.8; accumScoreThe = 0; goodScoreThe = -1000; dupholdScoreThe = 0, toolNum>=2

### 03.1 Data simulated by CNV-simulator

For simulation data evaluation, we will only use AS and DS to refine CNV set, other scores will not
be used. All the results are in `~/data3/project/CNVPipe/analysis-CNVSimulator/evaluation`. 
The script is in `~/data3/github-repo/CNVPipe/scripts/evalPerform3.py`.

    python ~/data3/github-repo/CNVPipe/scripts/evalPerform3.py ~/data3/project/CNVPipe/analysis-CNVSimulator/evaluation/v1.txt

## 04. Benchmark by real WGS data

Download and install the latest SRAtoolkit

    $ wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
    $ tar -xzvf sratoolkit.tar.gz

Add it to PATH.

Then assign the user repository which will be the location to store downloaded data.

    $ vdb-config -i

Test whether working

    $ fastq-dump --stdout -X 2 SRR390728
    $ prefetch -V
    prefetch : 3.0.1

Check the size of SRR

    $ vdb-dump SRR3602759 --info
    size   : 61,416,019,469

### AK1 dataset

The original data is from a Nature paper (De novo assembly and phasing of a Korean human genome). 
Two WGS data are available: SRR3602738 and SRR3602759, but we will only download and analyze
SRR3602759 because it has higher read depth, basically these two datasets are from the same sample.
The raw data is in `~/data3/raw_data/ncbi/sra/`. 

AK1 benchmark CNV set is downloaded from another paper (Trost et al, 2018, AJHG), 
at `~/data3/project/CNVPipe/analysis/truthSet/AK1_hg38_CNV_benchmark.txt`, 
this benchmark set is specifically for SRR3602759.

Start downloading SRR3602759 at Thu Nov 24 09:58:58 +08 2022
Check downloaded state at Fri Nov 25 07:58:57 +08 2022.

    $ prefetch SRR3602759 --max-size 62000000000

Extract fastq files

    $ fasterq-dump SRR3602759

Compress fastq files

    $ time gzip SRR3602759_1.fastq &
    $ time gzip SRR3602759_2.fastq &
    $ seqkit stats *.gz > stats.txt
    file                   format  type     num_seqs         sum_len  min_len  avg_len  max_len
    SRR3602759_1.fastq.gz  FASTQ   DNA   432,004,180  65,232,631,180      151      151      151
    SRR3602759_2.fastq.gz  FASTQ   DNA   432,004,180  65,232,631,180      151      151      151

Data is ready, pending for CNVPipe running.

### NA12878

NA12878 belongs to a 17 member CEPH pedigree, and has been sequenced by a lot of big projects, like
1000 genome project, Illumina's Platinum Genomes (50X and 200X), Broad Institute (NA12878 clone 
reference sequences, cell line, in my understanding) and NIH GIAB.

NA12878 benchmark CNV set1 is downloaded from another paper (Trost et al, 2018, AJHG) into
`~/data3/project/CNVPipe/analysis/truthSet/NA12878_hg38_CNV_benchmark1.txt`, 
original file is from `svclassify`, Trost filtered it to keep only deletions>1kb.

NA12878 benchmark CNV set2 is downloaded from Sun et al (2021, BMC Medical Genomics) additional 
file2 table 10 Set1, but it is based on hg19, we need to use liftover to transform it into hg38 
coordinate. We extracted set2 into 
`~/data3/project/CNVPipe/analysis/truthSet/NA12878_hg38_CNV_benchmark2.txt`.

NA12878 benchmark CNV set3 is also extracted from Sun et al (2021, BMC Medical Genomics) additional 
file2 table 10 Set3. `~/data3/project/CNVPipe/analysis/truthSet/NA12878_hg38_CNV_benchmark3.txt`.

Download NA12878 SVset from svclassify paper.
It is hg19 version, we then use liftover to transform it into hg38 version. (608 CNVs)

    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
    grep --color=never -v Chr Personalis_1000_Genomes_deduplicated_deletions.bed | awk '{OFS="\t"; print "chr"$1,$2,$3}' > NA12878-SVset.svclassify.hg19.bed
    ~/data3/biosoft/liftOver NA12878-SVset.svclassify.hg19.bed ~/data3/refs/hg38/bundle/liftover/hg19ToHg38.over.chain.gz NA12878-SVset.svclassify.hg38.bed umap.bed

1. 1000G NA12878 sequencing data 
(https://www.internationalgenome.org/data-portal/sample/NA12878)

    wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram
    wget -c -o 1000G.download.r2.log ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622457/SRR622457_2.fastq.gz &
    wget -c -o 1000G.download.r1.log ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622457/SRR622457_1.fastq.gz &

2. Illumina Platinum NA12878 data (50X)
(https://www.ebi.ac.uk/ena/browser/view/PRJEB3381?show=reads)
Note: 200X data is also available, but not useful in this project.
Start downloading at Thu Nov 24 19:54:22 +08 2022.
Check downloading finished at Sat Nov 26 19:40:14 +08 2022.

    wget -c -o illumina.download.log ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz &
    wget -c -o illumina.download.log ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz &
    $ seqkit stats *.gz > stats.txt
    file                  format  type     num_seqs         sum_len  min_len  avg_len  max_len
    ERR194147_1.fastq.gz  FASTQ   DNA   787,265,109  79,513,776,009      101      101      101
    ERR194147_2.fastq.gz  FASTQ   DNA   787,265,109  79,513,776,009      101      101      101

3. Broad Institute clone reference sequences (not sure if it fits our need)
(ftp://ftp.broadinstitute.org/pub/crd/NA12878_clones/)
    
    wget ftp://ftp.broadinstitute.org/distribution/crd/NA12878_clones/raw_data/A2925.1.Solexa-127359.aligned.duplicates_marked.bam
    wget ftp://ftp.broadinstitute.org/distribution/crd/NA12878_clones/raw_data/A2925.1.Solexa-127365.aligned.duplicates_marked.bam
    
4. NIH Genome in a Bottle (GIAB). Total reads depth reach 300X for NA12878 (14 libraries together)
(ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/).
`MGISEQ` under this directory might contain original fastq files of NA12878, but no README in this 
directory.
`BGISEQ500` contains PE50 and PE100 reads for NA12878 sequenced on BGISEQ-500 platform.
We start downloading Illumina Hiseq reads

    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam.bai
    
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/MGISEQ/NA12878_1/


5. BGI (paper PMID: 33849535), NA12878-1, sequenced in MGISEQ-2000, 197X.
(https://db.cngb.org/search/project/CNP0000813/)

TODO:
1. 找其他的有benchmark set的WGS数据集
2. 下载1000G中的normal sample
3. 其他的NA12878的benchmark set？( Comprehensive performance comparison of high-resolution array platforms for
genome-wide Copy Number Variation (CNV) analysis in humans)

### 1000G-normal

We decided to download 10 1000G WGS data as normal: HG00513, HG00525, HG00537, HG00556, HG00436,
HG00443, HG00448, HG00479, HG00619, HG00657, they are all "90 Han Chinese high coverage genomes".
In my understanding, this batch of data is generated by BGI-Shenzhen. Another batch of "Han Chinese
South" was produced by New York Genome Center for 1000 genome project phase3. The samples for these
two batch are the same, which are from cell lines produced by Coriell Institute.

    ascp -i /home/jhsun/.aspera/connect/etc/asperaweb_id_dsa.openssh -Tr -Q -l 100M -P33001 -L- era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR104/000/ERR1044780/ERR1044780_1.fastq.gz ./

Unfortunately, this file only have 80M reads, cannot reach 30X coverage, do not understand why it
is called high-coverage WGS. Over 12 paired-end fastq files for one single sample, their add-up 
coverage should be over 30X, but it is inconvenient to download and merge them all. Then we try to
download cram file from "Han Chinese South" of 1000G, because they only provide cram files,
an exact identical reference is needed if we want to extract reads from cram files.

    ascp -i /home/jhsun/.aspera/connect/etc/asperaweb_id_dsa.openssh -Tr -Q -l 100M -P33001 -L- era-fasp@fasp.sra.ebi.ac.uk:vol1/run/ERR323/ERR3239484/NA12778.final.cram ./ 2>download.log

Download reference genome, its index, dict and bwa index used in 1000G phase3.

    wget -c -o download.log ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa &

or use ascp

    ascp -i /home/jhsun/.aspera/connect/etc/asperaweb_id_dsa.openssh -Tr -Q -l 100M -P33001 -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa ./ 2>download.log

Transform cram to bam, then to fastq, compress with gzip. There are over 700M reads in NA12778.

    samtools fastq --reference ~/data/refs/hg38/1000G-ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -1 NA12778_1.fq -2 NA12778_2.fq NA12778.final.cram
    parallel samtools fastq --reference ~/data/refs/hg38/1000G-ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -1 {}_1.fq -2 {}_2.fq {}.final.cram ::: HG00557 HG00479 HG00524 HG00449 HG00421 HG00472 HG00595 HG00534 HG00580 HG00593

Download all 10 control samples.

    cat 1000G-CHS-10sampleList.txt | parallel -j 3 ascp -i /home/jhsun/.aspera/connect/etc/asperaweb_id_dsa.openssh -Tr -Q -l 100M -P33001 -L- era-fasp@fasp.sra.ebi.ac.uk:vol1/run/{} ./ 2>>download.{}.log

### CHM13

    prefetch -a "/home/jhsun/.aspera/connect/bin/ascp|/home/jhsun/.aspera/connect/etc/asperaweb_id_dsa.putty" SRR3986881 --max-size 50000000000

Download SV set from dbVar-nstd137. (1008 CNVs)

    wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd137.GRCh38.variant_call.tsv.gz
    cut -f 3,8,11,14 nstd137.GRCh38.variant_call.tsv | grep --color=never '^deletion' | awk '$4-$3>1000 {OFS="\t"; print "chr"$2,$3,$4,$1}' | sort -Vk 1 -k 2,3n > CHM13-SVset.1kb.bed

### GIAB - HG002

We downloaded HG00514 SV set from dbVar database with accession number nstd175. Only deletions can
be found in this file, no duplications. (538 CNVs)

    wget -c -r -o download.log ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392

    wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd175.GRCh38.variant_call.tsv.gz
    $cut -f 3 nstd175.GRCh38.variant_call.tsv | sortrn
        7359 insertion
        5500 deletion
        729 delins
        1 variant_call_type
        1 #VARIANT CALL
    cut -f 3,8,11,14 nstd175.GRCh38.variant_call.tsv | grep --color=never -e '^deletion' | awk '$4-$3>1000 {OFS="\t"; print "chr"$2,$3,$4,$1}' | sort -Vk 1 -k 2,3n


### HGSVC - HG00514,HG00733,NA19240

We downloaded HG00514 SV set from dbVar database with accession number nstd152. In this file, there
are deletions, duplications and copy number variations, in my opinion, the former two belong to the
later one. Therefore, here we will only include deletions and duplications

    wget -c -o download.log -i IGSR.fastq.ftp.txt &

    wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd152.GRCh38.variant_call.tsv.gz
    $cut -f 3 nstd152.GRCh38.variant_call.tsv | sortrn
        80281 deletion
        62457 insertion
        19549 complex substitution
        15763 alu insertion
        14733 alu deletion
        9243 duplication
        3016 sequence alteration
        2184 line1 insertion
        1968 inversion
        1896 line1 deletion
        1467 sva insertion
        675 copy number variation
        584 mobile element deletion
        472 mobile element insertion
        420 sva deletion
        110 herv deletion
        99 herv insertion
    grep --color=never HG00514 nstd152.GRCh38.variant_call.tsv | cut -f 3,8,11,14 | egrep --color=never '^(deletion|duplication)' | awk '$4-$3>1000 {OFS="\t"; print "chr"$2,$3,$4,$1}' | sort -Vk 1 -k 2,3n | grep -v inversion > HG00514-SVset.1kb.bed
    grep --color=never HG00733 nstd152.GRCh38.variant_call.tsv | cut -f 3,8,11,14 | egrep --color=never '^(deletion|duplication)' | awk '$4-$3>1000 {OFS="\t"; print "chr"$2,$3,$4,$1}' | sort -Vk 1 -k 2,3n | grep -v inversion > HG00733-SVset.1kb.bed
    grep --color=never NA19240 nstd152.GRCh38.variant_call.tsv | cut -f 3,8,11,14 | egrep --color=never '^(deletion|duplication)' | awk '$4-$3>1000 {OFS="\t"; print "chr"$2,$3,$4,$1}' | sort -Vk 1 -k 2,3n | grep -v inversion > NA19240-SVset.1kb.bed


### Run CNVPipe

Create soft symbolic link to `~/data3/project/CNVPipe/readAnalysis/samples`, copy config file to 
the directory, and create `samples.tsv` by `scripts/generate-table.py`, which was obtained from 
grenepipe.

    python ~/data/project/CNVPipe/scripts/generate-table.py /home/jhsun/data3/project/CNVPipe/realAnalysis/samples/ realAnalysis/samples.tsv

samples : NA12878 NA12778 HG00421 HG00449 HG00472 AK1
controls: HG00479 HG00524 HG00534 HG00557 HG00580 HG00593

Unfortunately, the results are very bad, none of the benchmark set CNV was identified by any of 
calling tools, which indicating either the benchmark set is wrong or the downloaded fastq data is 
not correct.

In addition, using 1000G normal as control to identify CNVs in AK1 and NA12878 may not be a good
idea, we should use simualted data directly from the reference genome.




## 05. Real patient WGS data analysis

### 05.1 Preparing for CNVPipe running

Put all sequencing data (PE150) under `/ubda/home/19044464r/data/`, then create soft link to `/ubda/home/19044464r/project/CNV-calling/analysis/fastq`. Then generate a table of samples by `~/project/CNVPipe/scripts/generate-table.py ~/project/CNV-calling/analysis/fastq/ ~/project/CNV-calling/analysis/samples.tsv`.

Here we also need to add control samples, we need add 'R' to label paired-end reads otherwise the script will throw an error, because E251_1 can be matched with E251_2 and E252_1. Finally we generated a table with both patient samples and control samples, which are ready for CNVPipe running.

#### Prepare CNVPipe config file.

Adjust the path of all reference files, including reference genome (analysis set from UCSC), CNVKit access (generated by `cnvkit.py access`) and refFlat file, 10x genomics sv blacklist for smoove to exclude (downloaded by `wget http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed`, and processed by `grep -E -w "chr([0-9]+|[XY])" low-mappability-track/sv_blacklist.bed | sort -Vk 1 -k 2,3n > low-mappability-track/hg38.sv_blacklist.bed`), delly mappibility file (downloaded from https://gear.embl.de/data/delly/, all three files are required), and finally gatk bundle files (downloaded from Google bitbuckit).

#### Installing snakemake by mamba

	wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
	bash Mambaforge-Linux-x86_64.sh
	mamba create -c conda-forge -c bioconda -n snakemake snakemake
	conda activate snakemake

### 05.2 CNVPipe running

Run CNVPipe, make sure `conda` and `mamba` is in PATH.

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --cluster qsub --jobs 6 --directory ~/project/CNV-calling/analysis/

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --cluster qsub --jobs 1 --directory ~/project/CNV-calling/analysis/ mapped/12719.bam

Only one core was used to do mapping, which is too slow, how to increase the threads in running?

I tested following commands, none of them worked, in all situation only one core was provided by UBDA.

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --cluster "qsub -t {threads}" --default-resources --jobs 1 --directory ~/project/CNV-calling/analysis/ mapped/12719.bam

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --cluster qsub --default-resources --jobs 1 --directory ~/project/CNV-calling/analysis/ mapped/12719.bam

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --cluster qsub --jobs 1 --cores 10 --directory ~/project/CNV-calling/analysis/ mapped/12719.bam

	snakemake --cluster qsub -j 1 --cores 30 --directory ~/project/CNV-calling/analysis/ mapped/12719.bam

Then we tried to use "profile", UBDA cluster management system is pbs-torque+MAUI. We download the snakemake pbs-torque profile template first and try run CNVPipe again. The torque profile is in `~/.config/snakemake/torque/`.

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --cluster qsub --jobs 1 profile ~/.config/snakemake/torque --directory ~/project/CNV-calling/analysis/ mapped/12719.bam

It still failed. We then tried the following one, finally it worked, we need to be careful that maximum cores in one node in UBDA queue q2s01 is 26, so we can not assign number of cores larger than this to snakemake rules.

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --profile ~/.config/snakemake/torque --directory ~/project/CNV-calling/analysis/ mapped/12719.bam

Run for all samples

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --profile ~/.config/snakemake/torque --directory ~/project/CNV-calling/analysis/ res/smoove/12719.bed

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --profile ~/.config/snakemake/torque --jobs 10 --directory ~/project/CNV-calling/analysis/

Final version for running all samples

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --profile torque --directory ~/project/CNV-calling/analysis/

Encountered an error: 

	$ cat logs/fastp/30196.log
	ERROR: igzip: encountered while decompressing file: /ubda/home/19044464r/project/CNV-calling/analysis/fastq/30196_2.fq.gz

Then checked integrity of file `30196_2.fq.gz`, it seems that this file is corrupted, but `30196_1.fq.gz` is good, which is weired. I tried to uncompress this file, same error was throw out.

	$ gzip -t fastq/30196_2.fq.gz
	gzip: fastq/30196_2.fq.gz: invalid compressed data--crc error
	gzip: fastq/30196_2.fq.gz: invalid compressed data--length error

Check the integrity of all fastq files

	for x in `cat sampleName`; do gzip -t fastq/${x} 2>>gzip.log; done

`29861_1.fq.gz` is also corrupted, so I skipped mapping them and used mapped BAM from previous results.

When doing mark duplicates by GATK, some samples encounter errors:
mapped/control-E264.raw.bam
mapped/23539.raw.bam
mapped/26473.raw.bam

The log files are normal, it is probably the problem of RAM. Then I decreased maximum RAM of GATK to 10G, everything goes well.

The final command line is as follow, as we are updating CNVPipe, we need to add `--rerun-triggers mtime` to trigger CNVPipe only by time stamp of specific files rather than by modification of rules.

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --profile torque --directory ~/project/CNV-calling/analysis/ all_bwamem --rerun-triggers mtime

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --profile torque --directory ~/project/CNV-calling/analysis/ all_mops --rerun-triggers mtime

	snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --profile torque --directory ~/project/CNV-calling/analysis/ --rerun-triggers mtime

Issues for rule `duphold_score`: the duphold version in smoove0.2.8 is too slow, thus we need to install the latest duphold and svtyper in `smoove.yaml`; second issue is that bam files of some samples are copied from previous mapping, the read group has changed, for example, sample 30237, previously is called 9-30237, duphold cannot work well on this name inconsistency. So we need to replace the read group information for these samples.

	gatk AddOrReplaceReadGroups -I input.bam -O output.bam -LB ${sample} -PL ILLUMINA -PU sample -SM sample -ID sample
	
Another is issue is that rule `delly_call_sv` runs too slow - it took almost 24 hours to run for one sample, which is unacceptable now, so I skip this step and directly call CNVs.

When using GATK VariantRecalibrator, sample 22562 has an error: A USER ERROR has occurred: Positive training model failed to converge. Basic idea is that the variance is not big enough, other samples work well, I think if we use gVCF calling, there will be no problem. I am trying to use less gaussians cores and see whether it can work.

It worked when I use `--max-gaussians 3`.


### 05.3 Test for SNP calling

One sample need around 4 hours

    snakemake --use-conda --conda-frontend mamba --conda-prefix /ubda/home/19044464r/biosoft/conda-env-cnvpipe --profile torque --directory ~/project/CNV-calling/analysis/ snps/gatk/12719.vqsr.vcf.gz --rerun-triggers mtime

Call SNPs for all samples.

Error occured in samples 29355,30021,28561, they all belong to previously mapped samples, while
samples mapped in this time runs well, why? Do I need to re-align the previous samples? That's very
resource-consuming, but it seems we have to do so.

### 05.4 Quality control

Check the quality of sequencing data by multiqc. Sample 30196 and 29861 is abnormal, we made up
their fastp json file. We will use previous QC results.

    cd ~/project/CNV-calling/analysis/cleaned
    multiqc -i fastp_result -n all_sample_quality ./
    multiqc, version 1.9

From multiqc results, we can see that sample 30361, 30382 and 30390 only have thousand reads, other
samples all have over 100M reads. It turns out I mis-linked the original fastq file, these samples
should be in Run10-2_V350034138 Lane4 with Adaptor 1-4, but I used Adaptor 13-16 to do soft link.
So I relink and rerun the calling and downstream analysis.
















## Principles of Varbin

We firstly downloaded the `k100.umap.bed` file from <https://bismap.hoffmanlab.org/>,
(Umap, Human, hg38, Single-read). This is a file showing goodzones of the genome, which
means the position can just be uniquely mapped. We then count the mappable positions in
each chromosome using `bin/varbin/mappable.py` and get `refs/k100.mappable.bed`. Then we construct
a target bin file using `bin/varbin/bin.boundary.py`, then sort this file based on chromosome using
`sort -Vk 1 -k 2,3n refs/bin.5k.test.txt > refs/bin.5k.boundary.txt`. Next, calculate the
gc content of target bins by using `bin/varbin/gc.content.py` and get `refs/bin.5k.boundary.gc.txt`.
Then we calculate the coverage of samples in target bins and correct the log2 ratio using
`bin/varbin/coverage.py`, in which `pysam` `bedcov` utils was used to calclulate the basecount in bin
regions. For example,

    time python bin/varbin/coverage.py refs/bin.5k.boundary.gc.txt data/simulation/1X/sample1.mkdup.bam analysis/read-depth/coverage/sample1.coverage.bed 1>1.log 2>2.log &
