# bulk CUT&Tag data analysis pipeline 
Standard pipeline to map CUT&amp;Tag data and generate genome browser tracks
Marek Bartosovic 

The pipeline assumes the standard Illumina naming convention of the fastq files:
```P27054_1001_S1_L001_R1_001.fastq.gz```

For more information in Illumina naming convention visit:
https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm

## Step 1 - clone the pipeline 

cd into directory with the project and clone the repo

```
cd /PATH/TO/MY_PROJECT
git clone https://github.com/bartosovic-lab/bulk_CT_pipeline
```

## Step 2 - generate config file 

Create config file in yaml format 
Minimal config file should look like this:
```
general:
  lanes:
    - L001
    - L002
  reads:
    - R1
    - R2
  bowtie2_index: /proj/snic2022-23-547/private/marek/reference/GRCh38/bowtie2_illumin_iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome # Your path to bowtie2 index
samples:
  P27054_1001:
    L001:
      R1:   /crex/proj/snic2022-23-547/private/marek/Thermo_Abs/P27054/P27054_1001/02-FASTQ/220817_A00621_0730_AH7F2GDRX2/P27054_1001_S1_L001_R1_001.fastq.gz   # Your path to L001 R1 file
      R2:   /crex/proj/snic2022-23-547/private/marek/Thermo_Abs/P27054/P27054_1001/02-FASTQ/220817_A00621_0730_AH7F2GDRX2/P27054_1001_S1_L001_R2_001.fastq.gz   # Your path to L001 R2 file
    L002:
      R1:   /crex/proj/snic2022-23-547/private/marek/Thermo_Abs/P27054/P27054_1001/02-FASTQ/220817_A00621_0730_AH7F2GDRX2/P27054_1001_S1_L002_R1_001.fastq.gz   # Your path to L002 R1 file
      R2:   /crex/proj/snic2022-23-547/private/marek/Thermo_Abs/P27054/P27054_1001/02-FASTQ/220817_A00621_0730_AH7F2GDRX2/P27054_1001_S1_L002_R2_001.fastq.gz   # Your path to L002 R2 file


```

You can use pipeline script to generate config:

```python3 bulk_CT_pipeline/scripts/generate_config.py /PATH/TO/FOLDER_WITH_FASTQ_FILES > config.yaml```

The script will search for all ```*.fastq.gz ``` files in the nested folder and from _LXXX_ and _RXXX_ determine R1 and R2 

#### Important:  After generating the config check that all needed files are present in thefile
#### Important: Change the path to the bowtie2 index within the pipeline config file
You can always edit the paths or any mistakes in the config file manualy or write your own script to generate config file

## Step 3 - setup slurm profile for batch job submission (Optional)

If running on rackham, create slurm profile for snakemake so the jobs can be submitted into the slurm queue by snakemake

For instructions how to do this follow:
https://github.com/Snakemake-Profiles/slurm

Here is how the setup looks for Rackham
```angular2html
profile_dir="${HOME}/.config/snakemake"
template="gh:Snakemake-Profiles/slurm"

mkdir -p "$profile_dir"

# Here you get prompted several times
cookiecutter --output-dir "$profile_dir" "$template"
You've downloaded /home/marek/.cookiecutters/slurm before. Is it okay to delete and re-download it? [yes]: yes
profile_name [slurm]: slurm_rackham
Select use_singularity:
1 - False
2 - True
Choose from 1, 2 [1]: 2
Select use_conda:
1 - False
2 - True
Choose from 1, 2 [1]: ^CAborted!
(base) [marek@rackham4 Thermo_Abs]$ cookiecutter --output-dir "$profile_dir" "$template"
You've downloaded /home/marek/.cookiecutters/slurm before. Is it okay to delete and re-download it? [yes]: yes
profile_name [slurm]: slurm_rackham
Select use_singularity:
1 - False
2 - True
Choose from 1, 2 [1]: 1
Select use_conda:
1 - False
2 - True
Choose from 1, 2 [1]: 2
jobs [500]:    
restart_times [0]: 0
max_status_checks_per_second [10]: 
max_jobs_per_second [10]: 
latency_wait [5]: 
Select print_shell_commands:
1 - False
2 - True
Choose from 1, 2 [1]: 2
sbatch_defaults []: account=snic2022-22-1062 time=0-24:00 partition=core ntasks=1 output=logs/slurm_%j.out error=logs/slurm_%j.err
cluster_sidecar_help [Use cluster sidecar. NB! Requires snakemake >= 7.0! Enter to continue...]: 
Select cluster_sidecar:
1 - yes
2 - no
Choose from 1, 2 [1]: 2
cluster_name []:        
cluster_jobname [%r_%w]: 
cluster_logpath [logs/slurm/%r/%j]: 
cluster_config_help [The use of cluster-config is discouraged. Rather, set snakemake CLI options in the profile configuration file (see snakemake documentation on best practices). Enter to continue...]: 
cluster_config []:
```
A folder with the config should appear at 
```$HOME/.config/snakemake/```

## Step 4 - Run the pipeline 
The pipeline is optimized to run on Rackham cluster

It is recommended to run the pipeline inside tmux session, for detailed instructions on tmux follow:
https://github.com/tmux/tmux/wiki/Getting-Started


Run the pipeline :

```snakemake --snakefile bulk_CT_pipeline/workflow/Snakefile.smk --cores 20 --jobs 100 -p --configfile config.yaml --rerun-incomplete --use-conda --profile slurm```