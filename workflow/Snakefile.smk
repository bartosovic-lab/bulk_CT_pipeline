include: 'prep.smk'

blacklist_path = download_blacklist(config)[0]

rule all:
  input:
    expand("{out}/{sample}/fastqc/{sample}_{lane}_{read}_fastqc.log", out = config['general']['output_dir'], sample = sample_list, lane = config['general']['lanes'], read=config['general']['reads']),
    expand("{out}/{sample}/{sample}_merged_nodup.bam", out = config['general']['output_dir'], sample = sample_list, lane = config['general']['lanes']),
    expand("{out}/{sample}/{sample}_merged_nodup.bw",out = config['general']['output_dir'], sample = sample_list),
    expand("{out}/{sample}/{sample}_peaks.broadPeak",out = config['general']['output_dir'], sample = sample_list),
    
rule fastqc:
  input:
    lambda wildcards: config['samples'][wildcards.sample][wildcards.lane][wildcards.read]
  output:
    "{out}/{sample}/fastqc/{sample}_{lane}_{read}_fastqc.log"
  params:
    outdir = "{out}/{sample}/fastqc/"
  conda: "../envs/bulkCT_fastqc.yaml"
  shell:
    "fastqc {input} --outdir {params.outdir} --nogroup > {output}"

rule trim_trim_galore:
  input:
    lambda wildcards: list(config['samples'][wildcards.sample][wildcards.lane].values())
  output:
    temp("{out}/{sample}/trimming/{sample}_{lane}_val_1.fq.gz"),
    temp("{out}/{sample}/trimming/{sample}_{lane}_val_2.fq.gz"),
  params:
    outdir = "{out}/{sample}/trimming/"
  conda: "../envs/bulkCT_trim.yaml"
  threads: 8
  resources:
    mem_mb = 8000
  shell:
    "trim_galore --cores {threads} --fastqc --paired -o {params.outdir} --basename {wildcards.sample}_{wildcards.lane} {input}"

rule map_bowtie2:
  input:
    read1 = "{out}/{sample}/trimming/{sample}_{lane}_val_1.fq.gz",
    read2 = "{out}/{sample}/trimming/{sample}_{lane}_val_2.fq.gz",
  output:
    bam = temp("{out}/{sample}/{sample}_{lane}_mapped.bam"),
    log = "{out}/{sample}/logs/{sample}_{lane}_bowtie2_map.log"
  params:
    index = config['general']['bowtie2_index']
  conda: "../envs/bulkCT_map.yaml"
  threads: 16
  resources:
    mem_mb= 32000
  shell:
    """
    bowtie2 --threads {threads} \
            --dovetail \
            -x {params.index} \
            -1 {input.read1} \
            -2 {input.read2} 2> {output.log} | samtools view -bS > {output.bam} 
    """
rule download_blacklist:
  output:
    blacklist_path
  conda: "../envs/bulkCT_download.yaml"
  shell:
    download_blacklist(config)[1]

rule bam_remove_blacklist_reads:
  input:
    bam = "{out}/{sample}/{sample}_{lane}_mapped.bam",
    blacklist = blacklist_path
  output:
    bam = temp("{out}/{sample}/{sample}_{lane}_mapped_blacklist_removed.bam")
  conda: "../envs/bulkCT_bedtools.yaml"
  threads: 1
  resources:
    mem_mb=32000
  shell:
    "bedtools intersect -v -abam {input.bam} -b {input.blacklist} > {output.bam}"

rule bam_sort_and_index:
  input:
    "{out}/{sample}/{sample}_{lane}_mapped_blacklist_removed.bam",
  output:
    bam_sorted = temp("{out}/{sample}/{sample}_{lane}_sorted.bam"),
    bam_index  = temp("{out}/{sample}/{sample}_{lane}_sorted.bam.bai"),
  threads: 16
  conda: "../envs/bulkCT_map.yaml"
  resources:
    mem_mb=32000
  shell:
    "samtools sort -o {output.bam_sorted} -@ {threads} {input} && samtools index {output.bam_sorted} "
    

rule merge_mapped:
  input:
    lambda wildcards: expand("{out}/{sample}/{sample}_{lane}_sorted.bam",out=config['general']['output_dir'],sample = wildcards.sample, lane = config['general']['lanes'])
  output:
    bam   = temp("{out}/{sample}/{sample}_merged.bam"),
    index = temp("{out}/{sample}/{sample}_merged.bam.bai"),
  threads: 16
  conda: "../envs/bulkCT_map.yaml"
  resources:
    mem_mb=32000
  shell:
    "samtools merge -@ {threads} {output.bam} {input} && samtools index {output.bam}"

rule picard_rmduplicates:
  input:
    bam = "{out}/{sample}/{sample}_merged.bam",
  output:
    bam     = "{out}/{sample}/{sample}_merged_nodup.bam",
    index   = "{out}/{sample}/{sample}_merged_nodup.bam.bai",
    metrics = "{out}/{sample}/{sample}_marked_dup_metrics.txt"
  conda: "../envs/bulkCT_picard.yaml"
  resources:
    mem_mb=32000
  shell:
    "PICARD=`which picard`;"
    "$PICARD MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES ;"
    "samtools index {output.bam}"

rule bam_to_bigwig:
  input:
    "{out}/{sample}/{sample}_{suffix}.bam",
  output:
    "{out}/{sample}/{sample}_{suffix}.bw",
  conda: "../envs/bulkCT_deeptools.yaml"
  threads: 16
  resources:
    mem_mb=32000
  shell:
    "bamCoverage -b {input} -o {output} -p {threads} --normalizeUsing RPKM --binSize 10 --ignoreDuplicates --extendReads 200 --centerReads --smoothLength 200 --effectiveGenomeSize 2652783500"


rule macs:
  input:
    bam = "{out}/{sample}/{sample}_merged_nodup.bam",
  output:
    peaks = "{out}/{sample}/{sample}_peaks.broadPeak",
  conda:
    "../envs/bulkCT_macs.yaml"
  threads: 1
  resources:
    mem_mb=32000
  shell:
    "macs2 callpeak -t {input} -f BAMPE -g hs -n {wildcards.sample} --outdir {wildcards.out}/{wildcards.sample} --broad --broad-cutoff 0.1 --nomodel --extsize 200 --keep-dup all --verbose 3"



