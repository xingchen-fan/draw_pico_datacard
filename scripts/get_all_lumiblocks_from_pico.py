#!/usr/bin/env python
import os
import glob
import ROOT
import pickle
import multiprocessing
import argparse

# Returns: runs_lumiblocks[run] = set(lumiblock)
def get_runs_and_lumiblocks(nanoaod_filename):
  #print(nanoaod_filename)
  runs_lumiblocks = {}

  # Slow method
  #nanoaod_file = ROOT.TFile.Open(nanoaod_filename, "READ")
  #for event in nanoaod_file.Events:
  #  if event.run not in runs_lumiblocks:
  #    runs_lumiblocks[event.run] = set()
  #  runs_lumiblocks[event.run].add(event.luminosityBlock)
  #print(runs_lumiblocks)

  # Fast method
  chain = ROOT.TChain("tree")
  chain.Add(nanoaod_filename)
  chain.SetEstimate(chain.GetEntries()+1)
  nEntries = chain.Draw("run:lumiblock","","goff")
  for iEntry in range(nEntries):
    run = int(chain.GetV1()[iEntry])
    lumiblock = int(chain.GetV2()[iEntry])
    #print(iEntry, run,lumiblock)
    if run not in runs_lumiblocks:
      runs_lumiblocks[run] = set()
    runs_lumiblocks[run].add(lumiblock)
  #print(runs_lumiblocks)
  
  return runs_lumiblocks

def pool_get_runs_and_lumiblocks(nanoaod_filename):
  print(nanoaod_filename)
  return [nanoaod_filename, get_runs_and_lumiblocks(nanoaod_filename)]

def get_dataset_name(nanoaod_filename):
  return os.path.basename(nanoaod_filename).split("__")[0].replace('raw_pico_','')

def combine_runs_lumiblocks(new_runs_lumiblocks, combined_runs_lumiblocks):
  for run in new_runs_lumiblocks:
    new_lumiblocks = new_runs_lumiblocks[run]
    if run not in combined_runs_lumiblocks:
      combined_runs_lumiblocks[run] = new_lumiblocks
    else:
      combined_runs_lumiblocks[run] = combined_runs_lumiblocks[run].union(new_lumiblocks)
# Example
#  runs_lumiblocks = {}
#  runs_lumiblocks[10] = set([0,1,2,3])
#  combine_runs_lumiblocks(runs_lumiblocks, all_runs_lumiblocks)
#  print(all_runs_lumiblocks)
#  runs_lumiblocks2 = {}
#  runs_lumiblocks2[10] = set([3,4,5,6])
#  combine_runs_lumiblocks(runs_lumiblocks2, all_runs_lumiblocks)
#  print(all_runs_lumiblocks)
#  runs_lumiblocks3 = {}
#  runs_lumiblocks3[11] = set([0,1,2])
#  combine_runs_lumiblocks(runs_lumiblocks3, all_runs_lumiblocks)
#  print(all_runs_lumiblocks)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='''\
Finds the run and lumiblock numbers for pico files. Saves the data to a pickle file in the OUTPUT_FOLDER. Uses multiple cores.
--------
Example commands:
get_all_lumiblocks_from_pico.py -p /net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath -y "2016,2017,2018" -o runs_lumiblocks_from_pico
get_all_lumiblocks_from_pico.py -p /net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath -y 2016 -o runs_lumiblocks_from_pico/2016
get_all_lumiblocks_from_pico.py -p /net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath -y 2017 -o runs_lumiblocks_from_pico/2017
get_all_lumiblocks_from_pico.py -p /net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath -y 2018 -o runs_lumiblocks_from_pico/2018
''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-p', '--pico_base_folder', required=True, default='/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath', help='Folder that contains the pico folders per year')
  parser.add_argument('-y', '--years_string', required=True, help='Years to find the run and lumiblock numbers. Split the years by a comma.')
  parser.add_argument('-o', '--output_folder', required=True, default='runs_lumiblocks', help='Output folder that holds pickle files of run and lumiblock numbers.')
  args = parser.parse_args()


  # Options
  #pico_base_path = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath"

  #years = [2016]
  #output_folder = 'runs_lumiblocks_from_pico/2016'
  #input_folder = 'runs_lumiblocks_from_pico/2016'
  #years = [2017]
  #output_folder = 'runs_lumiblocks_from_pico/2017'
  #input_folder = 'runs_lumiblocks_from_pico/2017'
  #years = [2018]
  #output_folder = 'runs_lumiblocks_from_pico/2018'
  #input_folder = 'runs_lumiblocks_from_pico/2018'

  pico_base_path = args.pico_base_folder
  years = args.years_string.split(',')
  output_folder = args.output_folder

  # Get pico files
  pico_filenames = []
  for year in years:
    pico_path = pico_base_path + "/" + str(year) + "/data/raw_pico"
    pico_filenames.extend(glob.glob(pico_path+"/*.root"))

  # files_runs_lumiblocks[filename][run] = set(lumiblock)
  files_runs_lumiblocks = {}
  # dataset_runs_lumiblocks[dataset][run] = set(lumiblock)
  dataset_runs_lumiblocks = {}
  # all_runs_lumiblocks[run] = set(lumiblock)
  all_runs_lumiblocks = {}

  #pico_filenames = pico_filenames[:10] # For testing

  # Multiprocess
  pool = multiprocessing.Pool()
  pool_result = pool.map(pool_get_runs_and_lumiblocks, pico_filenames)
  #print(pool_result)
  # Organize results
  for pico_filename, runs_lumiblocks in pool_result:
    files_runs_lumiblocks[pico_filename] = runs_lumiblocks
    dataset_name = get_dataset_name(pico_filename)
    if dataset_name not in dataset_runs_lumiblocks:
      dataset_runs_lumiblocks[dataset_name] = {}
    combine_runs_lumiblocks(runs_lumiblocks, dataset_runs_lumiblocks[dataset_name])
    combine_runs_lumiblocks(runs_lumiblocks, all_runs_lumiblocks)
    #print(pico_filename)
    #print(files_runs_lumiblocks)
    #print(dataset_runs_lumiblocks)
    #print(all_runs_lumiblocks)


  ## Single process
  #for pico_filename in pico_filenames:
  #  print(pico_filename)
  #  # Get run number and lumiblock
  #  runs_lumiblocks = get_runs_and_lumiblocks(pico_filenames[0])
  #  files_runs_lumiblocks[pico_filename] = runs_lumiblocks
  #  dataset_name = get_dataset_name(pico_filename)
  #  if dataset_name not in dataset_runs_lumiblocks:
  #    dataset_runs_lumiblocks[dataset_name] = {}
  #  combine_runs_lumiblocks(runs_lumiblocks, dataset_runs_lumiblocks[dataset_name])
  #  combine_runs_lumiblocks(runs_lumiblocks, all_runs_lumiblocks)
  #  #print(pico_filename)
  #  #print(files_runs_lumiblocks)
  #  #print(dataset_runs_lumiblocks)
  #  #print(all_runs_lumiblocks)

  # Save
  if not os.path.exists(output_folder):
    os.makedirs(output_folder)
  with open(output_folder+'/files_runs_lumiblocks.pickle', 'wb') as pickle_file:
    pickle.dump(files_runs_lumiblocks, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
  with open(output_folder+'/dataset_runs_lumiblocks.pickle', 'wb') as pickle_file:
    pickle.dump(dataset_runs_lumiblocks, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
  with open(output_folder+'/all_runs_lumiblocks.pickle', 'wb') as pickle_file:
    pickle.dump(all_runs_lumiblocks, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)

  ## Example
  ## Load
  #with open(input_folder + '/files_runs_lumiblocks.pickle', 'rb') as pickle_file:
  #  files_runs_lumiblocks = pickle.load(pickle_file)
  #with open(input_folder + '/dataset_runs_lumiblocks.pickle', 'rb') as pickle_file:
  #  dataset_runs_lumiblocks = pickle.load(pickle_file)
  #with open(input_folder + '/all_runs_lumiblocks.pickle', 'rb') as pickle_file:
  #  all_runs_lumiblocks = pickle.load(pickle_file)
  ##print(files_runs_lumiblocks)
  ##print(dataset_runs_lumiblocks)
  ##print(all_runs_lumiblocks)
