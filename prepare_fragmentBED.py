# -*- coding: utf-8 -*-

import os
import re
import time
from multiprocessing import Pool

def file_path(tmp_path):
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)


genebody_dir = "data/bed/genebody_5hmC"
promoter_dir = "data/bed/promoter_5hmC"

fastqc = "fastqc"
bowtie = "bowtie2"
picard="picard.jar"
samtools = "samtools"
bedtools = "bedtools"
# macs2 = "macs2"
bdg2bw = "bedGraphToBigWig"

index = "data/ref/GRCh37.p13.genome"
blacklist = "/home/asus/Documents/reference/hg19_consensusBlacklist.bed"
fai_file = "data/ref/GRCh37.p13.genome.fa.fai"


# os.system("source ~/.bash_profile")

def Fastqc(fastq1, fastq2, outdir):
    os.system("{fastqc} {fastq1} -o {outdir} -t 2 ".format(fastqc=fastqc, fastq1=fastq1, outdir=outdir))
    os.system("{fastqc} {fastq2} -o {outdir} -t 2".format(fastqc=fastqc, fastq2=fastq2, outdir=outdir))

def makePEbed(sample, bed_file, workdir):
    print sample, ": now start to make PEbed", time.asctime()
    Bed = open(bed_file, "r")
    PEbed = os.path.join(workdir, "paired_dedup_fragment_%s.bed" % sample)
    outBed = open(PEbed, "w")
    line = Bed.readline()
    info = line.strip().split("\t")
    chrID = info[0]
    start = info[1]
    end = info[2]
    readname = info[3].split("/")[0]
    line = Bed.readline()
    while line:
        info = line.strip().split("\t")
        newchrID = info[0]
        newstart = info[1]
        newend = info[2]
        newreadname = info[3].split("/")[0]
        if chrID == newchrID and newreadname == readname:
            if float(newstart) < float(start):
                start = newstart
            if float(newend) > float(end):
                end = newend
            newline = chrID + '\t' + start + "\t" + end + "\t" + readname + "\t" + info[-2] + "\t" + info[-1] + "\n"
            outBed.write(newline)
            line = Bed.readline()
            if line:
                info = line.strip().split("\t")
                chrID = info[0]
                start = info[1]
                end = info[2]
                readname = info[3].split("/")[0]
            else:
                break
        else:
            chrID = newchrID
            start = newstart
            end = newend
            readname = newreadname
        line = Bed.readline()
    Bed.close()
    outBed.close()
    os.system("sort -k1,1 -k2,2n %s -o %s" % (PEbed, PEbed))

def mapping(fastq1, fastq2, samplename):
    print samplename, ": now start to mapping", time.asctime()
    workdir = os.path.join(bam_dir, samplename)
    file_path(workdir)
    os.chdir(workdir)
    os.system("{bowtie} --threads 10 --rg \"PL:ILLUMINA\" --rg \"SM:{sampleID}\" -x {index} -1 {fq1} -2 {fq2} -S {sampleID}.sam".format(
            bowtie=bowtie, sampleID=samplename, index=index, fq1=fastq1, fq2=fastq2))
    os.system("{samtools} view -f 2 -F 1548 -q 30 -h {sampleID}.sam > filtered_{sampleID}.sam".format(samtools=samtools, sampleID=samplename))
    os.system("{samtools} view -b -S filtered_{sampleID}.sam > filtered_{sampleID}.bam".format(samtools=samtools, sampleID=samplename))
    os.system("{samtools} view -b -S {sampleID}.sam > {sampleID}.bam".format(samtools=samtools, sampleID=samplename))
    os.system("java -jar {picard} SortSam \
                    INPUT=filtered_{sampleID}.bam \
                    OUTPUT=sorted_filtered_{sampleID}.bam \
                    SORT_ORDER=coordinate".format(picard=picard, sampleID=samplename))
    os.system("java -jar {picard} MarkDuplicates \
                    REMOVE_DUPLICATES=true \
                    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 \
                    INPUT=sorted_filtered_{sampleID}.bam \
                    OUTPUT=dedup_sorted_filtered_{sampleID}.bam \
                    METRICS_FILE=dedup_sorted_filtered_{sampleID}.metrics".format(picard=picard, sampleID=samplename))
    bedfile = workdir + "/" + "dedup_sorted_filtered_%s.bed" % samplename
    os.system("{bedtools} bamtobed -i filtered_{sampleID}.bam > filtered_{sampleID}.bam.bed".format(bedtools=bedtools, sampleID=samplename))
    os.system("{bedtools} bamtobed -i dedup_sorted_filtered_{sampleID}.bam > {bedfile}".format(bedtools=bedtools, sampleID=samplename, bedfile=bedfile))
    os.system("sort -k4,4 -k1,1 -k2,2n %s -o %s" % (bedfile, bedfile))
    os.system("{samtools} flagstat {sampleID}.bam > {sampleID}.bam.flagstat".format(samtools=samtools, sampleID=samplename))
    os.system("{samtools} flagstat filtered_{sampleID}.bam > filtered_{sampleID}.bam.flagstat".format(samtools=samtools, sampleID=samplename))
    os.system("{samtools} flagstat dedup_sorted_filtered_{sampleID}.bam > dedup_sorted_filtered_{sampleID}.bam.flagstat".format(
                samtools=samtools, sampleID=samplename))
    os.system("rm {sampleID}.sam filtered_{sampleID}.sam".format(sampleID=samplename))
    print samplename, "mapping have been done!", time.asctime()
    makePEbed(samplename, bedfile, workdir)
    PEbed = os.path.join(workdir, "paired_dedup_fragment_%s.bed" % samplename)

def trimming(samplename, fastq1, fastq2, cutlength):
    Fastqc(fastq1, fastq2, fastqc_dir)
    filter_fq1 = os.path.join(cutadapter_dir, "%s_1_paired.fq.gz" % samplename)
    filter_fq2 = os.path.join(cutadapter_dir, "%s_2_paired.fq.gz" % samplename)
    print samplename, ": now start to cutapter!"
    os.system("java -jar /mnt/raid/tools/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -threads 2 {read1} {read2} \
                    {filter_fq1} {cutadapter}/{name}_1_unpaired.fq.gz \
                    {filter_fq2} {cutadapter}/{name}_2_unpaired.fq.gz \
                    ILLUMINACLIP://mnt/raid/cfDNA20190116/all_TruSeq3-PE.fa:2:30:10 LEADING:20 \
                    TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 HEADCROP:{cutlength}".format(read1=fastq1, read2=fastq2,
                                                                                 filter_fq1=filter_fq1,
                                                                                 filter_fq2=filter_fq2,
                                                                                 cutadapter=cutadapter_dir,
                                                                                 name=samplename, cutlength=cutlength))
    print samplename, ": now start to fastqc!"
    Fastqc(filter_fq1, filter_fq2, cutadapter_QC)
    mapping(filter_fq1, filter_fq2, samplename)

def make_sample_info(sample_info_file):
    sample_list = {}
    f = open(sample_info_file, "r")
    head = f.readline()
    line = f.readline()
    while line:
        rawname = line.split()[0]
        newname = line.split()[1]
        type = line.split()[4]
        sample_list[rawname] = [newname, type]
        line = f.readline()
    f.close()
    return(sample_list)


projectdir = "data"
file_path(projectdir)
fastqc_dir = os.path.join(projectdir, "fastqc")
file_path(fastqc_dir)
cutadapter_dir = os.path.join(projectdir, "cutadapter_data")
file_path(cutadapter_dir)
cutadapter_QC = os.path.join(projectdir, "cutadapter_QC")
file_path(cutadapter_QC)
bam_dir = os.path.join(projectdir, "mapping_result")
file_path(bam_dir)


sample_list = make_sample_info("data/sampleInfo.csv")


pools = Pool(8)
for sample in os.listdir("data/raw"):
    rawname = sample.split("_BKDL")[0]
    if sample_list.has_key(rawname) :
        samplename = sample_list[rawname][0]
        if sample_list[rawname][1] == "5hmC":
            print '5hmC',rawname,samplename
            fastq1 = "data/raw" + "/" + sample + "/" + sample + "_1.clean.fq.gz"
            fastq2 = "data/raw" + "/" + sample + "/" + sample + "_2.clean.fq.gz"
            # trimming(samplename, fastq1, fastq2, "15")
            pools.apply_async(trimming, args=(samplename, fastq1, fastq2, "15"))
            # break
    else:
        print rawname,'is not included.'

pools.close()
pools.join()
del pools