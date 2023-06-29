import os

sample_list         = [x for x in config['samples'].keys()]
lanes               = config['general']['lanes']
reads               = config['general']['reads']
files_dict          = {x: [config['samples'][x][l][r] for l in lanes for r in reads] for x in sample_list}
files_dict_basename = {key: [os.path.basename(files_dict[key][i]) for i in range(0,len(files_dict[key]))] for key in files_dict }
files_dict_dirname  = {key: [os.path.dirname(files_dict[key][i]) for i in range(0,len(files_dict[key]))] for key in files_dict }

# print(sample_list)
# print(files_dict)
# print(files_dict_basename)
