import os

sample_list         = [x for x in config['samples'].keys()]
lanes               = config['general']['lanes']
reads               = config['general']['reads']
files_dict          = {x: [config['samples'][x][l][r] for l in lanes for r in reads] for x in sample_list}
files_dict_basename = {key: [os.path.basename(files_dict[key][i]) for i in range(0,len(files_dict[key]))] for key in files_dict }
files_dict_dirname  = {key: [os.path.dirname(files_dict[key][i]) for i in range(0,len(files_dict[key]))] for key in files_dict }

def download_blacklist(config):
    blacklist_url = {   # This can be changed if needed
        'hg19': 'https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz',
        'hg38': 'https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz',
        'mm10': 'https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz'
    }

    if config['general']['genome'] in blacklist_url.keys():
        blacklist_url = blacklist_url[config['general']['genome']]
        blacklist_file = os.path.join(config['general']['output_dir'], 'blacklist_{}.bed.gz'.format(config['general']['genome']))
        cmd = 'wget -O {} {}'.format(blacklist_file, blacklist_url)
    else:
        print('*** Error: No blacklist available for genome {}'.format(config['general']['genome']))
        blacklist_file = None
        sys.exit(1)
    return blacklist_file, cmd