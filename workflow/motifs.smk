import os 

include: "complexity.smk"
include: "prep.smk"
include: "Snakefile.smk" 

rule all_motifs:
    input:
        expand("{out}/{sample}/motifs/meme_out/meme.html", out = config['general']['output_dir'], sample = sample_list),

rule macs_narrow:
  input:
    bam = "{out}/{sample}/{sample}_merged_nodup.bam",
  output:
    peaks = "{out}/{sample}/{sample}_peaks.narrowPeak",
  conda:
    "../envs/bulkCT_macs.yaml"
  threads: 1
  resources:
    mem_mb=32000
  shell:
    "macs2 callpeak -t {input} -f BAMPE -g hs -n {wildcards.sample} --outdir {wildcards.out}/{wildcards.sample} --nomodel --extsize 200 --keep-dup 1 --verbose 3"

rule resize_peaks:
    input:
        "{out}/{sample}/{sample}_summits.bed",
    output:
        "{out}/{sample}/motifs/{sample}_peaks_resized.bed",
    conda:
        "../envs/bulkCT.yaml"
    params:
        peak_size = 200,
        npeaks = 10000,
    shell:
        """
        set +o pipefail
        awk 'BEGIN{{OFS="\\t"}} {{$2 = $2 - {params.peak_size}/2; $3 = $3 + {params.peak_size}/2; print $0}}' {input}  | 
        keep_standard_chromosomes.sh | 
        filter_invalid_bed.sh | 
        sort -k5,5nr | 
        head -{params.npeaks} | 
        sort -k1,1 -k2,2n > {output}
        """

rule get_fasta_resized_peaks:
    input:
        "{out}/{sample}/motifs/{sample}_peaks_resized.bed",
    output:
        "{out}/{sample}/motifs/{sample}_peaks_resized.fasta",
    conda:
        "../envs/bulkCT.yaml"
    params:
        genome_fa = config['general']['genome_fa'],
    shell:
        """
        bedtools getfasta -fi {params.genome_fa} -bed {input} -fo {output} -name -s
        """

rule run_meme:
    input:
        "{out}/{sample}/motifs/{sample}_peaks_resized.fasta",
    output:
        "{out}/{sample}/motifs/meme_out/meme.html",
    conda:
        "../envs/bulkCT_meme.yaml"
    params:
        meme_out = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),
        jaspar_db = config['general']['jaspar_db'],
    shell:
        """
        meme-chip -oc {params.meme_out} -db {params.jaspar_db} -maxw 15 -meme-nmotifs 10 -dna {input};
        touch {output}
        """