#!/usr/bin/env/python3

import os, sys
import glob
# import yaml

path = sys.argv[1]
extension = "/**/*.fastq.gz"

class seq_file:
  def __init__(self, path):
    self.basename    = os.path.basename(path)
    self.dirname     = os.path.dirname(path)
    self.path        = path
    self.seq_id      = self.basename.split("_")[0]
    self.sample_name = self.basename.split("_")[1]
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
      config_dict[f.seq_id + "_" + f.sample_name].append(f)
    except KeyError:
      config_dict[f.seq_id + "_" + f.sample_name] = [f]
    
  return(config_dict)


def format_config(config_dict):
  for key in sorted(config_dict.keys()):
    print("{}:".format(key))
    curr_files = config_dict[key]
    lanes      = [x.lane for x in curr_files]
    for lane in sorted(set(lanes)):
      print("  {}:".format(lane))
      reads_ls = [curr_files[x].read for x in range(0,len(curr_files)) if curr_files[x].lane == lane]
      for read in reads_ls:
        print("    {}:".format(read))
        print("      {}".format("\n".join([curr_files[x].path for x in range(0,len(curr_files)) if curr_files[x].lane == lane and curr_files[x].read == read])))
      



def main(path, extension):
  files_list  = [seq_file(x) for x in glob.glob(path + extension,recursive=True)]
  config_dict = generate_config_dict(files_list)
  format_config(config_dict)

main(path, extension)


