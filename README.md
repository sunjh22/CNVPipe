# Develop a new CNV calling method - CNVPipe

Objectives of CNVPipe

1. easy to use, high sensitivity and low FDR for CNV calling
2. integrate comprehensive plotting tools for visualization
3. integrate CNV annotation tools: pathogenicity prediction

## Run CNVPipe snakemake

    snakemake -np --directory analysis/ all_lumpy --rerun-triggers mtime
    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data/biosoft/conda-env-cnvpipe --cores 10 --directory analysis/ all_lumpy --rerun-triggers mtime

CNVPipe-token: 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

    git clone https://github.com/sunjh22/CNVPipe.git

Use bio1 server to develop with small data

Use bio2 server to run big data

## 01. Prepare test data - simulation

We want to simulate three patient samples with CNV in them and three control samples,
all downsampled to 0.01X for fast testing of our pipeline.

### Use SCNVSim to simulate SV in genome

    curl -OL https://sourceforge.net/projects/scnvsim/files/latest/download
    unzip download
    cd data/
    time java -jar ~/data/biosoft/scnvsim_1.3.1/normgenomsim_1.3.1.jar -v ~/data/refs/hg38/analysisSet/hg38.analysisSet.size.rmMT.txt -n ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -o ./simulation/simuGenome -s 0.00000001 -r 0.000000005 2>2.log &
    time java -jar ~/data/biosoft/scnvsim_1.3.1/tumorgenomsim_1.3.1.jar -v ~/data/refs/hg38/analysisSet/hg38.analysisSet.size.rmMT.txt -n ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -k ../refs/repeatmask.txt -i ./simulation/simuGenome/normal_snvindelsim.vcf -o simulation/simuGenome -s 0.001 -l 0.001 -c 100 -a tumor1 1>1.log 2>2.log &
    $ SCNVSim Version: 1.3.1

### Use ART to simulate reads

    ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HS25 -f 1 -i simulation/simuGenome/normal_snvindelsim.fasta -l 150 -m 500 -na -o simulation/1X/control_ -p -s 10 > 1.log 2> 2.log &
    ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HS25 -f 1 -i simulation/simuGenome/tumor1_svcnvsim_clone_0.fasta -l 150 -m 500 -na -o simulation/1X/sample2_ -p -s 10 > 1.log 2> 2.log &
    gzip --fast sample_1.fq
    ART Version: v2.5.8

### Use seqkit to downsample

    cd data/simulation/1X
    parallel seqkit sample -p 0.01 -s 100 control{}_1.fq.gz -o ../../../analysis/samples-control/control{}_1.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 control{}_2.fq.gz -o ../../../analysis/samples-control/control{}_2.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 sample{}_1.fq.gz -o ../../../analysis/samples/sample{}_1.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 sample{}_2.fq.gz -o ../../../analysis/samples/sample{}_2.fq.gz ::: 1 2 3 &
    Seqkit Version: 0.10.1

`seqkit stat analysis/samples-control/*`, 116619 reads in control samples, 201149 in test samples.

## 02. Prepare snakemake config file

data: samples, control-samples, genome, refFlat

## 03. Prepare main Snakefile

Include `.smk` file step by step, testing their availability during development.

## 04. Write common.smk

Import config file, get fastq files (depending on se or pe), get patient sample and 
control sample names and set them to be global parameter. wildcard_constraints (sample 
names and control-sample names), valid filename and filepath.

## 05. Write pre-processing.smk

This script includes another two scripts:
- `fastp.smk` for cleaning reads, use snakemake wrapper;
- `bwamem.smk` for mapping reads using bwa-mem, also use snakemake wrapper.

## 06. Write smk for Read-depth based CNV-calling methods

We need to be clear about some features of these methods:
1. require normal samples? if yes, how many?
2. call germline CNV sample by sample or in batch
3. exclude low-mappability region? require reference file?

### 06.1. CNVKit

Require normal samples (no limit on the number), reference genome and exclude low-map region.
Call germline CNV in batch.

Include four rules:
- `autobin` for estimating optimal bin resolution;
- `batch` for calling CNV for WGS data in batch mode;
- `segmetrics` for segmenting CNVs based on confidence interval;
- `call` for transforming log2 depth ratio to integer copy numbers.

And a simple shell command using `awk` to extract informative columns
into another file. Columns include genomeic coordinates, integer CN,
log2 ratio, depth, probe and weight.

The refFlat and access file should be provided for annotation and
discarding reads in low-mappbility region.

#### Principle of CNVKit

The first step is to get accessible regions. Find low mappability track: 1. get from
rqfu, in `~/data/refs/hg38/low-mappability-track/hg38.badRegions.bed`,
which includes all centromere, telomere and heterochromatin regions; 2. download from
10x genomics - black list of SV <http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed> to
`~/data/refs/hg38/low-mappability-track/sv_blacklist.bed`, re-format it using
`grep -E -w "chr([0-9]+|[XY])" low-mappability-track/sv_blacklist.bed | sort -Vk 1 -k 2,3n > low-mappability-track/hg38.sv_blacklist.bed`.
Then use `cnvkit.py access hg38.analysisSet.fa -x hg38.badRegions.bed -x hg38.sv_blacklist.bed -o access-excludes.hg38.bed`
to get accessible regions. `cnvkit.py access` command will automatically find
sequence of N's in genome and filter regions with long N's. Then we get target bin file
by `cnvkit.py target refs/access-excludes.hg38.bed --annotate ~/data/refs/hg38/bundle/CNVKit/refFlat.txt --avg-size 20000 --split -o refs/target.bed`,
which will then be used to calculate coverage of samples in these bins.

### 06.2. cnvpytor

Do not require normal sample and reference genome and do not exclude low-map region.
This method solely depends on comparing read depth between adjacent region (?) and regions
with similar GC content. Call germline CNV sample by sample.

Include one rule with four commands:
- `rd` for extracting reads from sample bam file;
- `his` for calculating reads depth at certain resolution;
- `partition` for segmenting the bins;
- `call` for calling CNV.

The results include "CNVtype, CNVregion, CNVsize, CNVlevel, eval1, eval2, eval3, eval4, q0, pN, dG"
- `CNVlevel` - normalized read depth value, estimated integer CN value should be `round(2*CNVlevel)`;
- `eval1` - e-value (p-value multiplied by genome size divided by bin size) calculated using t-test statistics between RD statistics in the region and global - set to [0 - 0.00001];
- `eval2` - e-value from the probability of RD values within the region to be in the tails of a gaussian distribution of binned RD - set to [0 - 0.00001];
- `q0` - fraction of reads mapped with q0 quality in call region - not useful for filtering since most reads have good quality;
- `pN` - fraction of reference genome gaps (Ns) in call region - set to [0 - 0.5];
- `dG` - distance from closest large (>100bp) gap in reference genome - set to [>100kb];

Run cnvpytor

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -chrom $(seq -f 'chr%g' 1 22) chrX chrY -rd sample1.bam &
    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -his 20000 &
    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -partition 20000 &
    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -call 20000 > ../read-depth/cnvnator/sample1.2k.tsv &
    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -plot rd 20000 -o ../read-depth/cnvnator/sample1.png &

### 06.3. Control-FREEC

It claims requiring control samples, but I cannot include them into analysis. Require reference genome,
do not exclude low-map region. Call germline CNV sample by sample.

After reading developer's manual and testing for several times, I found it hard to include normal sample in this software,
which is really weired. So now I just need to focus on step 2 and 3. Build a config file first and run the software for
each sample.

Run freec

    ~/data/biosoft/FREEC/src/freec -conf analysis/temp/freec/config_WGS.txt -sample analysis/mapped/sample1.bam
    ~/data/biosoft/FREEC/src/freec -conf analysis/temp/freec/config_WGS.txt -sample analysis/mapped/sample1.bam -control analysis/mapped/control1.bam

`assess_significance.R` is not working - may need some debug.

`makeGraph.R` scirpt plot the dot plot of copy number in each window in each chromosome. gain is red, loss is blue.

### 06.4. cn.MOPS

Require at least 6 normal samples, do not require reference genome and do not exclude low-map region.
Call CNV in batch. This method call CNV by comparing read depth between patient sample and control sample.

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

## 07. Merge the results from four CNV calling tools

Except merging the results from four calling tools, we try to remove genomic bad regions 
(centromere, telomere and heterochromatin). The bad region file is in 
`~/data/refs/hg38/low-mappability-track/hg38.badRegions.bed`. Besides, we need to make sure the 
merged file adapt to the input format of `CNVfilteR`.

We tried to merge CNV results from all tools including cnvkit, cnvpytor, freec, mops, smoove 
and delly.

??? What is heterochromatin

## 08. Call SNPs from data

The reason we do not use GATK for this task is that it
requires many known-variant files, which increase the
preparation load if for clinical usages. I will consider
add that option in future development.

### 08.1 freebayes

Call SNP sample by sample (the problem is the accuracy). Then filter based on quality score
and read depth, this is a critical step. We should set different threshold for different depth
data. The accuracy of SNP calling will affect the correction of CNV in later step.

    freebayes -f ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa mapped/sample1.bam > snps/freebayes/sample1.raw.vcf
    bcftools filter -O v -o snps/freebayes/sample1.filtered.vcf -s LOWQUAL -e 'QUAL<10 || FMT/DP <5' --SnpGap 5 --set-GTs . snps/freebayes/sample1.raw.vcf
    bcftools view -v snps snps/freebayes/sample1.filtered.vcf > snps/freebayes/sample1.snp.vcf

### 08.2 samtools + bcftools

Almost the same as freebayes.

### 08.3 GATK

## 9. Call CNV based on read-pair and split-read

### 9.1 Lumpy

Lumpy is written by C.
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

### 9.2 Smoove

There is a Lumpy wrapper called `smoove`, which is claimed to be faster. Just have try and see
if it can replace lumpy.

An bad region file was downloaded according to the author of smoove. Overall length is 119,556,880,
which is almost half shorter than Heng Li's. But when we remove non-cannonical chromosomes, its 
length became 3,042,685, which is definitely not correct. 
The file is in `low-mappability-track/smoove/hg38.exclude.bed`.

    wget https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed

Our own bad region file is in `low-mappability-track/hg38.badRegions.bed`, which includes centromere,
telomere and heterochromatin, the length is 199,765,358.

Another low-map track is downloaded from 10x genomics in 
`low-mappability-track/10x/hg38.sv_blacklist.bed`, its overall length is 212,765,070.
I think this one might be the best one to use.

Test `smoove` calling

    smoove call --outdir analysis/temp/smoove/ --exclude ~/data/refs/hg38/exclude.cnvnator_100bp.GRCh38.20170403.bed --name sample4 --fasta ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -p 1 --genotype analysis/mapped/sample4.bam

After comparing the results from Lumpy and Smoove, we found no big difference, so we will use smoove
to replace lumpy in our CNVPipe.

Use `duphold` to genotype. Three additional fields will be added into `FORMAT` per sample:
DHFC - "fold-change for the variant depth relative to the rest of the chromosome the variant was found on"
DHFFC - "fold-change for the variant depth relative to Flanking regions"
DHBFC - "fold-change for the variant depth relative to bins in the genome with similar GC-content"

    smoove genotype -d -x -p 1 --name sample4 --outdir analysis/temp/smoove-genotype/ --fasta ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa --vcf analysis/temp/smoove/sample4-smoove.genotyped.vcf.gz analysis/mapped/sample4.bam

The authors suggest use DHFFC<0.7 to filter deletion and DHBFC>1.3 to filter duplication

    bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] > 1.3)' analysis/temp/smoove-genotype/sample4-smoove.genotyped.vcf.gz

Use `bcftools` to extract informative columns

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL[\t%DHFC\t%DHFFC\t%DHBFC]\n' analysis/temp/smoove-genotype/sample4-smoove.genotyped.vcf.gz

Follow author's suggestion to filter CNV based on DHFFC and DHBFC by `awk`

    awk '$4=="DEL" && $7<0.7 {print $0} $4=="DUP" && $8>1.3 {print $0}' analysis/res/smoove/sample1.bed | wcl

`duphold` also support SNP vcf file for filtering, just like `CNVfilteR`.

### 9.3 Delly

Delly is written in C++. It is primarily designed for calling structural variant, both
germline and somatic (cancers). It also supports CNV calling now.

Mappbility files are mandantory for Delly to call CNV, we downloaded map files
from [here]<https://gear.embl.de/data/delly/> at Sat Oct 15 17:08:30 +08 2022. All three
files are required.

Delly by default divide genome into 10kb-mappable bins, but we can set window size by `-i`
parameter. `-u` is for segmentation, `-l` is for using delly SV calling to refine breakpoints.

So we try to call SV first, and use its output as input for CNV calling. 
`delly call` only use paired-end and split-read information, and read-depth is not used.

    delly call -g ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -x ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed -o analysis/temp/delly/sample1.sv.bcf analysis/mapped/sample1.bam

Use `delly cnv` to call CNV

    delly cnv -u -i 20000 -g ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -m ~/data/refs/hg38/bundle/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz -c analysis/temp/delly/sample1.cov.gz -o analysis/temp/delly/sample1.cnv.bcf analysis/mapped/sample1.bam

Use `duphold` to do genotyping

    duphold -v analysis/temp/delly/sample4.filtered.bcf -b analysis/mapped/sample4.bam -f ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa -o analysis/temp/delly/sample4.duphold.vcf

`delly classify` will select CNVs with 'PASS' (filter low-quality CNVs), bad thing is that it
transforms all quality score into 10000, so we decided to filter CNVs by grep and keep
the original quality score. One thing need to be noticed is that in some cases, value in `FILTER`
field is not matched with `FT` tag in `format` field. One could be PASS while the other one
could be LowQual. We just use the `FILTER`=="PASS" to select all good-quality CNVs. (strict mode)

Use `bcftools` to extract informative columns

    bcftools query -f '%FILTER\t%CHROM\t%POS\t%INFO/END[\t%RDCN\t%DHFC\t%DHFFC\t%DHBFC]\n' analysis/temp/delly/sample4.duphold.vcf | grep 'PASS' | cut -f 2- > analysis/res/delly/sample4.bed

In `INFO` field, `MP` means fraction of mappable positions; `GCF` means GC fraction (added by 
duphold); these should not be used for CNV filtering.

We also need to filter CNV based on DHFFC and DHBFC as stated in `smoove` section.

    awk '$4<2 && $7<0.7 {print $0} $4>2 && $8>1.3 {print $0}' analysis/res/delly/sample1.bed | wcl

Some R scripts were provided for plotting CNV, which could be a reference for later usage.

    Rscript ~/data/biosoft/delly/R/rd.R analysis/temp/delly/sample1.cov.gz analysis/temp/delly/sample1.segment.bed

    Version: Delly 1.1.5

## 10. Merge results from Lumpy and Delly



## 11. BAF-based CNV correction

`CNVfilteR`

How to install and use it in snakemake? Or I need to write a simplified version in Python by myself. 
Or I need to find another tool that was implemented in Python. Directly download in R script if not 
installed.

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

## 12. Combine results from read-depth based and read-pair based methods

## 13. CNV pathogenicity prediction

### 13.1 ClassifyCNV

ClassifyCNV cannot be downloaded from bioconda, and its result directory is weired,
we need to figure out how to make it work in our CNVPipe.

### 13.2 X-CNV

XCNV only supports hg19, which is outdated, so we give up using this tool to predict
CNV pathogenicity.

## 14. visualization


TODO:
1. merge some conda envs yaml files, put common tools together, only those require specific environment get a single env























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
