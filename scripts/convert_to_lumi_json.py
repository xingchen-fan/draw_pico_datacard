#!/bin/env python
import pickle
import collections
import json
import os
import argparse

def get_dataset_name(nanoaod_filename):
  return os.path.basename(nanoaod_filename).split("__")[0]

def combine_runs_lumiblocks(new_runs_lumiblocks, combined_runs_lumiblocks):
  for run in new_runs_lumiblocks:
    new_lumiblocks = new_runs_lumiblocks[run]
    if run not in combined_runs_lumiblocks:
      combined_runs_lumiblocks[run] = new_lumiblocks
    else:
      combined_runs_lumiblocks[run] = combined_runs_lumiblocks[run].union(new_lumiblocks)

def reduce_lumiblock_range(lumiblocks):
  if len(lumiblocks) == 0: print("Crashing because lumiblocks is empty")
  sorted_lumiblocks = sorted(lumiblocks)

  reduced_lumiblock_ranges = []

  start_lumiblock = sorted_lumiblocks[0]
  track_lumiblock = sorted_lumiblocks[0]
  for lumiblock in sorted_lumiblocks[1:]:
    # continous
    if track_lumiblock+1 == lumiblock:
      track_lumiblock = lumiblock
    # Non-continous
    else:
      reduced_lumiblock_ranges.append([start_lumiblock,track_lumiblock])
      start_lumiblock = lumiblock
      track_lumiblock = lumiblock
  # Add last
  # Only one lumiblock
  # Multiple lumiblocks
  reduced_lumiblock_ranges.append([start_lumiblock,track_lumiblock])
  #print(sorted_lumiblocks)
  #print(reduced_lumiblock_ranges)
  return reduced_lumiblock_ranges

def expand_lumiblock_range(lumiblock_ranges):
  lumiblocks = set()
  for lumiblock_range in lumiblock_ranges:
    start_lumiblock = lumiblock_range[0]
    end_lumiblock = lumiblock_range[1]
    for lumiblock in range(start_lumiblock, end_lumiblock+1):
      lumiblocks.add(lumiblock)
  return lumiblocks
      
def get_runs_lumiblocks(json_filename):
  with open(json_filename) as json_file:
    raw_json = json.load(json_file)
  runs_lumiblocks = {}
  for run in raw_json:
    lumiblock_range = raw_json[run]
    runs_lumiblocks[int(run)] = expand_lumiblock_range(lumiblock_range)
  return runs_lumiblocks

def filter_runs_lumiblocks(runs_lumiblocks, golden_runs_lumiblocks):
  filtered_runs_lumiblocks = {}
  for run in runs_lumiblocks:
    if run not in golden_runs_lumiblocks: continue
    #print(run)
    #print(runs_lumiblocks[run])
    #print(golden_runs_lumiblocks[run])
    #print(runs_lumiblocks[run] - golden_runs_lumiblocks[run])
    #print(runs_lumiblocks[run].intersection(golden_runs_lumiblocks[run]))
    filtered_lumiblocks = runs_lumiblocks[run].intersection(golden_runs_lumiblocks[run])
    if len(filtered_lumiblocks) != 0:
      filtered_runs_lumiblocks[run] = filtered_lumiblocks
    #if len(golden_runs_lumiblocks[run] - runs_lumiblocks[run]) != 0:
    #  print("Missing run: "+str(run)+ " lumiblock: "+str(sorted(golden_runs_lumiblocks[run] - runs_lumiblocks[run])))
  return filtered_runs_lumiblocks

# Returns: a - b
def subtract_runs_lumiblocks(runs_lumiblocks_a, runs_lumiblocks_b):
  subtracted_runs_lumiblocks = {}
  for run in runs_lumiblocks_a:
    if run not in runs_lumiblocks_b: subtracted_runs_lumiblocks[run] = runs_lumiblocks_a[run]
    else: 
      lumiblocks = runs_lumiblocks_a[run] - runs_lumiblocks_b[run]
      if len(lumiblocks) != 0: subtracted_runs_lumiblocks[run] = lumiblocks
  return subtracted_runs_lumiblocks

def make_json_string(runs_lumiblocks):
  lumiblock_json_string = '{'
  for run in sorted(runs_lumiblocks.keys()):
    lumiblock_json_string += '"'+str(run)+'": '
    lumiblock_json_string += str(reduce_lumiblock_range(runs_lumiblocks[run]))+',\n'
    # print(expand_lumiblock_range(reduce_lumiblock_range(runs_lumiblocks[run])).symmetric_difference(runs_lumiblocks[run]))
  lumiblock_json_string = lumiblock_json_string[:-2]+'}'
  return lumiblock_json_string

def print_files_runs_lumiblocks(runs_lumiblocks, files_runs_lumiblocks):
  for run in runs_lumiblocks:
    for lumiblock in runs_lumiblocks[run]:
      filenames = []
      for filename in files_runs_lumiblocks:
        if run in files_runs_lumiblocks[filename]:
          if lumiblock in files_runs_lumiblocks[filename][run]:
            filenames.append(filename)
      if len(filenames) != 0:
        print("Run: "+str(run)+" lumiblock: "+str(lumiblock))
        for filename in filenames:
          print("  "+filename)

def get_golden_filename(year):
  if year== 2016: return 'runs_lumiblocks/golden/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'
  elif year == 2017: return 'runs_lumiblocks/golden/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
  elif year == 2018: return 'runs_lumiblocks/golden/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'

def get_old_golden_filename(year):
  if year== 2016: return 'runs_lumiblocks/old_golden/golden_Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16.json'
  elif year == 2017: return 'runs_lumiblocks/old_golden/golden_Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17.json'
  elif year == 2018: return 'runs_lumiblocks/old_golden/golden_Cert_314472-325175_13TeV_PromptReco_Collisions18.json'

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='''\
Converts pickle files that contain run and lumiblock numbers to a json file.
The pickle file is generated by get_all_lumiblocks.py or get_all_lumiblocks_from_pico.py.
The output json file can be used with LUMI POG's luminosity calculation tool.
--------
Example commands:
convert_to_lumi_json.py -y 2016 -i runs_lumiblocks -o runs_lumiblocks -g txt/golden/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
convert_to_lumi_json.py -y 2017 -i runs_lumiblocks -o runs_lumiblocks -g txt/golden/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
convert_to_lumi_json.py -y 2018 -i runs_lumiblocks -o runs_lumiblocks -g txt/golden/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
  
Command to use old json file:
convert_to_lumi_json.py -y 2018 -i runs_lumiblocks -o runs_lumiblocks -g txt/golden/golden_Cert_314472-325175_13TeV_PromptReco_Collisions18.json -t _old_golden
  
Command to use pico file:
convert_to_lumi_json.py -y 2016 -i runs_lumiblocks_from_pico -o runs_lumiblocks_from_pico -g txt/golden/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-y', '--year_string', required=True, help='Select year to put in the json file')
  parser.add_argument('-i', '--input_folder', required=True, help='Folder that contains the pickle file')
  parser.add_argument('-o', '--output_folder', required=True, help='Folder that contains the output json file')
  parser.add_argument('-g', '--golden_filename', required=True, help='Golden JSON that holds valid run and lumiblock numbers')
  parser.add_argument('-t', '--tag_output', default='', help='Add a postfix to the output json file')
  args = parser.parse_args()


  ## Options
  ## year can be 2016, 2017, 2018
  #year = 2016
  ## From pico or from NanoAOD
  #from_pico = False
  ## Use old golden Json file or new golden Json file
  #use_old_golden = True

  ## Set input files
  ## Requires pickle files from get_all_lumiblocks.py
  #if from_pico: input_folder = 'runs_lumiblocks_from_pico/'+str(year)
  #else: input_folder = 'runs_lumiblocks/'+str(year)
  ## Requires golden jsons files in runs_lumiblocks/golden/
  #if use_old_golden: golden_filename = get_old_golden_filename(year)
  #else: golden_filename = get_golden_filename(year)

  #old_string = ''
  #if use_old_golden: old_string = '_old_golden'

  year = int(args.year_string)
  input_folder = args.input_folder+'/'+str(year)
  golden_filename = args.golden_filename
  output_folder = args.output_folder
  tag_output = args.tag_output

  ## Set output files
  #if from_pico:
  #  all_output_filename = 'runs_lumiblocks_from_pico/nanoaod_'+str(year)+old_string+'.json'
  #  missing_all_output_filename = 'runs_lumiblocks_from_pico/nanoaod_'+str(year)+old_string+'_missing.json'
  #  met_output_filename = 'runs_lumiblocks_from_pico/nanoaod_'+str(year)+old_string+'_met.json'
  #  no_met_output_filename = 'runs_lumiblocks_from_pico/nanoaod_'+str(year)+old_string+'_no_met.json'
  #  dataset_output_folder = 'runs_lumiblocks_from_pico/datasets'+old_string+'/'
  #else:
  #  all_output_filename = 'runs_lumiblocks/nanoaod_'+str(year)+old_string+'.json'
  #  missing_all_output_filename = 'runs_lumiblocks/nanoaod_'+str(year)+old_string+'_missing.json'
  #  met_output_filename = 'runs_lumiblocks/nanoaod_'+str(year)+old_string+'_met.json'
  #  no_met_output_filename = 'runs_lumiblocks/nanoaod_'+str(year)+old_string+'_no_met.json'
  #  dataset_output_folder = 'runs_lumiblocks/datasets'+old_string+'/'

  all_output_filename = output_folder+'/nanoaod_'+str(year)+tag_output+'.json'
  missing_all_output_filename = output_folder+'/nanoaod_'+str(year)+tag_output+'_missing.json'
  met_output_filename = output_folder+'/nanoaod_'+str(year)+tag_output+'_met.json'
  no_met_output_filename = output_folder+'/nanoaod_'+str(year)+tag_output+'_no_met.json'
  dataset_output_folder = output_folder+'/datasets'+tag_output+'/'


  # Open files
  with open(input_folder + '/files_runs_lumiblocks.pickle', 'rb') as pickle_file:
    files_runs_lumiblocks = pickle.load(pickle_file)
  with open(input_folder + '/dataset_runs_lumiblocks.pickle', 'rb') as pickle_file:
    dataset_runs_lumiblocks = pickle.load(pickle_file)
  with open(input_folder + '/all_runs_lumiblocks.pickle', 'rb') as pickle_file:
    all_runs_lumiblocks = pickle.load(pickle_file)
  # files_runs_lumiblocks[filename][run] = set(lumiblock)
  #print(files_runs_lumiblocks)
  # dataset_runs_lumiblocks[dataset][run] = set(lumiblock)
  #print(dataset_runs_lumiblocks)
  # all_runs_lumiblocks[run] = set(lumiblock)
  #print(all_runs_lumiblocks)

  ## Regather all_runs_lumiblocks with only JetHT, MET, SingleElectron, SingleMuon, and EGamma
  #all_runs_lumiblocks = {}
  #primary_datasets = ['JetHT', 'MET', 'SingleElectron', 'SingleMuon', 'EGamma']
  #for filename in files_runs_lumiblocks:
  #  dataset_name = get_dataset_name(filename)
  #  if dataset_name not in primary_datasets: continue
  #  runs_lumiblocks = files_runs_lumiblocks[filename]
  #  combine_runs_lumiblocks(runs_lumiblocks, all_runs_lumiblocks)

  golden_runs_lumiblocks = get_runs_lumiblocks(golden_filename)

  # Filter runs_lumiblocks
  valid_all_runs_lumiblocks = filter_runs_lumiblocks(all_runs_lumiblocks, golden_runs_lumiblocks)

  # Filter datasets
  valid_datasets_runs_lumiblocks = {}
  for dataset in dataset_runs_lumiblocks:
    valid_datasets_runs_lumiblocks[dataset] = filter_runs_lumiblocks(dataset_runs_lumiblocks[dataset], golden_runs_lumiblocks)

  # Extra filters
  valid_met_runs_lumiblocks = filter_runs_lumiblocks(dataset_runs_lumiblocks['MET'], golden_runs_lumiblocks)
  # Find missing met lumiblocks respect to golden json
  no_met_runs_lumiblocks = subtract_runs_lumiblocks(golden_runs_lumiblocks, valid_met_runs_lumiblocks)
  # Find missing met lumiblocks respect to all_runs_lumiblocks
  no_met_runs_lumiblocks = subtract_runs_lumiblocks(valid_all_runs_lumiblocks, valid_met_runs_lumiblocks)
  # Find files for runs_lumiblocks
  #print_files_runs_lumiblocks(no_met_runs_lumiblocks, files_runs_lumiblocks)
  # Find missing lumiblocks
  missing_all_runs_lumiblocks = subtract_runs_lumiblocks(golden_runs_lumiblocks, valid_all_runs_lumiblocks)


  # Save files
  with open(all_output_filename,'w') as nanoaod_json: nanoaod_json.write(make_json_string(valid_all_runs_lumiblocks))
  with open(missing_all_output_filename,'w') as nanoaod_json: nanoaod_json.write(make_json_string(missing_all_runs_lumiblocks))
  with open(no_met_output_filename,'w') as nanoaod_json: nanoaod_json.write(make_json_string(no_met_runs_lumiblocks))

  if not os.path.exists(dataset_output_folder): os.makedirs(dataset_output_folder)
  for dataset in valid_datasets_runs_lumiblocks:
    with open(dataset_output_folder+"/nanoaod_"+str(year)+"_"+dataset+".json",'w') as nanoaod_json: 
      nanoaod_json.write(make_json_string(valid_datasets_runs_lumiblocks[dataset]))

