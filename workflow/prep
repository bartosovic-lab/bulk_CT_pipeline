import os
import itertools

configfile: "/data/proj/GCB_MB/CT_version2/analysis/config/config.yaml"


sample_list         = [x for x in config['samples'].keys()]
files_dict          = {x: [config['samples'][x][l][r] for l in config['samples'][x].keys() for r in config['samples'][x][l].keys()] for x in sample_list}
files_dict_basename = {key: [os.path.basename(files_dict[key][i]) for i in range(0,len(files_dict[key]))] for key in files_dict }
files_dict_dirname  = {key: [os.path.dirname(files_dict[key][i]) for i in range(0,len(files_dict[key]))] for key in files_dict }

#print(sample_list)
#print(files_dict)
#print(files_dict_basename)
