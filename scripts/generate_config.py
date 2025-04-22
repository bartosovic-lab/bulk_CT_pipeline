#!/usr/bin/env/python3

import os, sys
import glob
import argparse

parser = argparse.ArgumentParser(description="Generate a config file for the pipeline.")
parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input directory containing the fastq files.')
parser.add_argument('-r', '--read_R3_as_R2', action= 'store_true', help='Use this flag to report R3 as R2 if sequenced together with single-cell libraries on the same flowcell')
parser.add_argument('-e', '--extension', type=str, default='*.fastq.gz', help='File extension of the input files. Default is *.fastq.gz')
args = parser.parse_args()


class seq_file:
  def __init__(self, path):
    self.path        = os.path.abspath(path)
    self.basename    = os.path.basename(self.path)
    self.dirname     = os.path.dirname(self.path)
    self.seq_id      = self.basename.split("_")[0]
    self.seq_no      = self.basename.split("_")[1]
    self.sample_name = self.seq_id + "_" + self.seq_no
    self.sample_no   = self.basename.split("_")[2]
    self.lane        = self.basename.split("_")[3]
    self.read        = self.basename.split("_")[4]
    self.gziped      = self.is_gzip()
  def is_gzip(self):
    if self.basename.endswith(".gz"):
      return True
    else:
      return False

def generate_config_dict(files_list):
  config_dict = {}
  for f in files_list:
    try:
      config_dict[f.sample_name].append(f)
    except KeyError:
      config_dict[f.sample_name] = [f]
    
  return(config_dict)


def format_config(config_dict):
  print("samples:")
  for key in sorted(config_dict.keys()):
    print("  {}:".format(key))
    curr_files = config_dict[key]
    lanes      = sorted([x.lane for x in curr_files])
    for lane in sorted(set(lanes)):
      print("    {}:".format(lane))
      reads_ls = sorted([curr_files[x].read for x in range(0,len(curr_files)) if curr_files[x].lane == lane])
      for read in reads_ls:
        if args.read_R3_as_R2:
          if read == "R1":
            read_report = "R1"
          elif read == "R2":
            continue
          elif read == "R3":
            read_report = "R2"
        else:
          read_report = read
        print("      {}:".format(read_report),end = " ")
        print("  {}".format("\n".join([curr_files[x].path for x in range(0,len(curr_files)) if curr_files[x].lane == lane and curr_files[x].read == read])))


def print_general_config():
  print("general:")
  print("  lanes:\n    - L001\n    - L002")
  print("  reads:\n    - R1\n    - R2")
  print("  bowtie2_index: FILL_IN_PATH_TO_BOWTIE2_INDEX_HERE")

def main(path, extension):
  print(path+extension)
  files_list  = [seq_file(x) for x in glob.glob(path + "/**/" + extension,recursive=True)]
  config_dict = generate_config_dict(files_list)
  print_general_config()
  format_config(config_dict)

main(args.input, args.extension)


