#!/bin/env python
import os
import sys
import argparse
import subprocess
import ROOT
import json
import yaml

def runCommand(command, log_filename = None, printOutput=True, noRun=False):
  print("Will run below command: ")
  print(command)
  out_string = ''
  err_string = ''
  if not noRun:
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    out_string = out.decode("utf-8").rstrip()
    err_string = err.decode("utf-8").rstrip()
    if printOutput:
      if out_string != '':
        print("######## Output log ########")
        print (out_string)
        print("######## Output log end #########")
      if err_string != '':
        print("######## Error log ########")
        print(err_string)
        print("######## Error log end ########")
    if log_filename:
      with open(log_filename,'w') as log_file:
        print("Saving out and err to "+log_filename)
        if out_string != '':
          log_file.write("######## Output log ########\n")
          log_file.write(out_string+"\n")
          log_file.write("######## Output log end #########\n")
        if err_string != '':
            log_file.write("######## Error log ########\n")
            log_file.write(err_string+"\n")
            log_file.write("######## Error log end #########\n")
  return out_string, err_string

def set_POI_parameter_names():
  POI_parameter_names = []
  xsig_variant = ['xsig0', 'xsig1']
  ysig_varient = ['ysig']
  met_varient = ['met0', 'met1', 'met2', 'met3']
  drmax_varient = ['drmax0', 'drmax1']
  for xsig in xsig_variant:
    for ysig in ysig_varient:
      for met in met_varient:
        for drmax in drmax_varient:
          POI_parameter_names.append('rp_'+'_'.join((xsig, ysig, met, drmax)))
  return POI_parameter_names

def set_prefit_mask_parameter_names():
  mask_parameter_names = []
  xsig_variant = ['xsig0', 'xsig1']
  ysig_varient = ['ysig']
  met_varient = ['met0', 'met1', 'met2', 'met3']
  drmax_varient = ['drmax0', 'drmax1']
  for xsig in xsig_variant:
    for ysig in ysig_varient:
      for met in met_varient:
        for drmax in drmax_varient:
          mask_parameter_names.append('mask_'+'_'.join((xsig, ysig, met, drmax))+"=1")
  return mask_parameter_names

def set_mask_parameter_names():
  mask_parameter_names = []
  xsig_variant = ['xsig0', 'xsig1', 'xbkg']
  ysig_varient = ['ysig', 'ybkg']
  met_varient = ['met0', 'met1', 'met2', 'met3']
  drmax_varient = ['drmax0', 'drmax1']
  for xsig in xsig_variant:
    for ysig in ysig_varient:
      for met in met_varient:
        for drmax in drmax_varient:
          mask_parameter_names.append('mask_'+'_'.join((xsig, ysig, met, drmax))+"=1")
  return mask_parameter_names

# Returns POI_results
# If fit result is not good or a POI parameter is missing, returns {}
def get_RooFitResult(root_filename, rooFitResult_name, POI_parameter_names):
  root_file = ROOT.TFile(root_filename)
  rooFitResult = root_file.Get(rooFitResult_name)

  # Check if fit is valid
  isValid = True
  nStatus = rooFitResult.numStatusHistory()
  for iFit in range(nStatus):
    if rooFitResult.statusCodeHistory(iFit) != 0:
      print("[Error] Fit status "+str(rooFitResult.statusLabelHistory(iFit))+" is "+rooFitResult.statusCodeHistory(iFit)+" for RooFitResult "+rooFitResult_name+" in "+root_filename)
      isValid = False
  if rooFitResult.covQual() != 3:
    #  Document: https://root.cern.ch/root/htmldoc/guides/minuit2/Minuit2.pdf
    #  covariance matrix meaning: https://root.cern/doc/master/classROOT_1_1Minuit2_1_1Minuit2Minimizer.html#a44b31a43d1eee371dee1792f5dfaee2b
    #    status = -1 : not available (inversion failed or Hesse failed)
    #    status = 0 : available but not positive defined
    #    status = 1 : covariance only approximate
    #    status = 2 : full matrix but forced pos def
    #    status = 3 : full accurate matrix
    print("[Error] Covariance quality is "+str(rooFitResult.covQual())+" for RooFitResult "+rooFitResult_name+" in "+root_filename)
    isValid = False
  if not isValid: return {}

  # Get POI results
  # POI_results[POI_parameter] = [fit value, down uncertainty (difference), diff up uncertainty (difference)]
  POI_results = {}
  fitParameters = rooFitResult.floatParsFinal()
  for POI_parameter_name in POI_parameter_names:
    if not fitParameters.find(POI_parameter_name):
      print("[Error] POI parameter: "+POI_parameter_name+" is not in RooFitResult "+rooFitResult_name+" in "+root_filename)
      return {}
    POI_parameter = fitParameters.find(POI_parameter_name)
    POI_results[POI_parameter_name] = [POI_parameter.getValV(), POI_parameter.getAsymErrorLo(), POI_parameter.getAsymErrorHi()]

  return POI_results

def rp(abcd_name, plane_name):
  return 'rp_'+abcd_to_xy(abcd_name)+'_'+plane_name
def kappa(abcd_name, plane_name):
  return 'kappa_'+abcd_to_xy(abcd_name)+'_'+plane_name
def yield_range(data_yield):
  if data_yield < 2: return '[0,20.0]'
  return '[0,'+str(data_yield*10)+']'

def abcd_to_xy(abcd_name):
  abcd_to_xy_dict = {'A0':'xsig0_ysig', 'A1':'xsig1_ysig', 'B0':'xsig0_ybkg', 'B1':'xsig1_ybkg', 'C': 'xbkg_ysig', 'D': 'xbkg_ybkg'}
  return abcd_to_xy_dict[abcd_name]

def xy_to_abcd(xy_name):
  abcd_to_xy_dict = {'A0':'xsig0_ysig', 'A1':'xsig1_ysig', 'B0':'xsig0_ybkg', 'B1':'xsig1_ybkg', 'C': 'xbkg_ysig', 'D': 'xbkg_ybkg'}
  xy_to_abcd_dict = dict([(value, key) for key, value in abcd_to_xy_dict.items()])
  return xy_to_abcd_dict[xy_name]

def convert_to_inverse_abcd_datacard(abcd_datacard_name, inverse_abcd_datacard_name):
  # Collect data yield for each bin
  # dataYields[xy_name] = yield
  dataYields = {}
  with open(abcd_datacard_name) as abcd_datacard:
    for line in abcd_datacard:
      if 'bin' in line:
        name_split = line.split()
      if 'Observation' in line:
        observation_split = line.split()
        break
    for iName, name in enumerate(name_split):
      dataYields[name] = observation_split[iName]

  # Replace equation to become inverse datacard
  inverse_abcd_datacard_string = ''
  with open(abcd_datacard_name) as abcd_datacard:
    for line in abcd_datacard:
      if 'rp_' in line[0:3]:
        rp_name = line.split()[0]
        xy_name = '_'.join(rp_name.split('_')[1:3])
        plane_name = '_'.join(rp_name.split('_')[3:])
        abcd_name = xy_to_abcd(xy_name)
        if abcd_name == 'D' or abcd_name == 'C':
          inverse_abcd_datacard_string += line
        elif abcd_name == 'B0':
          split_line = line.split()
          split_line[4] = '@0*@1/@2/@3'
          split_line[5] = rp('A0',plane_name)+','+rp('D',plane_name)+','+rp('C',plane_name)+','+kappa('A0',plane_name)
          inverse_abcd_datacard_string += ' '.join(split_line)+'\n'
        elif abcd_name == 'B1':
          split_line = line.split()
          split_line[4] = '@0*@1/@2/@3'
          split_line[5] = rp('A1',plane_name)+','+rp('D',plane_name)+','+rp('C',plane_name)+','+kappa('A1',plane_name)
          inverse_abcd_datacard_string += ' '.join(split_line)+'\n'
        elif abcd_name == 'A0':
          split_line = line.split()
          data_yield = float(dataYields[abcd_to_xy('A0')+'_'+plane_name])
          if data_yield == 0: data_yield = 1. # For fit stability
          split_line[4] = str(data_yield)
          split_line[5] = yield_range(data_yield)
          inverse_abcd_datacard_string += ' '.join(split_line)+'\n'
        elif abcd_name == 'A1':
          split_line = line.split()
          data_yield = float(dataYields[abcd_to_xy('A1')+'_'+plane_name])
          if data_yield == 0: data_yield = 1. # For fit stability
          split_line[4] = str(data_yield)
          split_line[5] = yield_range(data_yield)
          inverse_abcd_datacard_string += ' '.join(split_line)+'\n'
      else:
        inverse_abcd_datacard_string += line

  with open(inverse_abcd_datacard_name,'w') as inverse_abcd_datacard:
    inverse_abcd_datacard.write(inverse_abcd_datacard_string)
  print('Wrote '+inverse_abcd_datacard_name)

def get_limit(root_filename):
  root_file = ROOT.TFile(root_filename)
  tree = root_file.Get('limit')
  tree.GetEntry(0)
  return tree.limit
  

# process_inverse_datacard.py -d cards_1d_kappa_priority1/datacard-TChiHH_mChi-400_mLSP-0_Tune_2016,2017,2018_priority1_resolved.txt -i abcd_inverse_datacard.txt -o datacard_results.json
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Processes a inverse ABCD datacard. \n'\
  +'Outputs a JSON with following format: \n'\
  +'  POI_results[POI_parameter_name][prefit/postfit] = [fit value, down uncertainty (difference), diff up uncertainty (difference)]\n'\
  +'Can be read by following code in python:\n'\
  +'  with open(OUTPUT_JSON) as POI_result_file\n'\
  +'    POI_result = json.load(POI_result_file)'
  ,formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-d', '--datacard', required=True, help='ABCD datacard')
  parser.add_argument('-i', '--inverse_datacard', required=True, help='Set the inverse ABCD datacard. If filename does not exist, will make datacard')
  parser.add_argument('-o', '--output_json', required=True, help='Processed output')
  parser.add_argument('-f', '--fakeRun', action="store_true", help='Do not run commands. Commands are only printed.')
  args = parser.parse_args()

  #noRun = True
  noRun = args.fakeRun
  printOutput = False

  # POI_parameter_names: [POI parameter name]
  POI_parameter_names = set_POI_parameter_names()

  # Generate inverse datacard if not exist
  if not os.path.exists(args.inverse_datacard):
    print(args.inverse_datacard+' does not exists. Will generate inverse datacard from '+args.datacard)
    convert_to_inverse_abcd_datacard(args.datacard, args.inverse_datacard)

  # Do prefit
  print("\n\n****************Doing prefit****************")
  # Convert txt datacard to root datacard for prefit that requires masking
  # text2workspace.py doesn't seem to be able to detect errors
  runCommand('text2workspace.py '+args.inverse_datacard+' --channel-masks', None, printOutput, noRun)
  inverse_datacard_root = args.inverse_datacard.replace('.txt','.root')
  # prefit_mask_parameter_names: [mask parameter name]
  prefit_mask_parameter_names = set_prefit_mask_parameter_names()
  prefit_combine_command = 'combine -M MultiDimFit --saveFitResult '\
  +'--redefineSignalPOIs='+','.join(POI_parameter_names)\
  +' --setParameters r=0,'+','.join(prefit_mask_parameter_names)\
  +' --freezeParameters r '+inverse_datacard_root+' --name _prefit'
  prefit_log_filename = 'prefit_log.txt'
  runCommand(prefit_combine_command, prefit_log_filename, printOutput, noRun)
  # Get results from 'multidimfit_prefit.root'
  # prefit_POI_results[POI_parameter_name] = [fit value, down uncertainty (difference), diff up uncertainty (difference)]
  prefit_POI_results = get_RooFitResult('multidimfit_prefit.root', 'fit_mdf', POI_parameter_names)

  # Do postfit
  print("\n\n****************Doing postfit****************")
  postfit_combine_command = 'combine -M MultiDimFit --saveFitResult '\
  +'--redefineSignalPOIs='+','.join(POI_parameter_names)\
  +' --setParameters r=0'\
  +' --freezeParameters r '+args.inverse_datacard+' --name _postfit'
  postfit_log_filename = 'postfit_log.txt'
  runCommand(postfit_combine_command, postfit_log_filename, printOutput, noRun)
  # Get results from 'multidimfit_postfit.root'
  # postfit_POI_results[POI_parameter] = [fit value, down uncertainty (difference), diff up uncertainty (difference)]
  postfit_POI_results = get_RooFitResult('multidimfit_postfit.root', 'fit_mdf', POI_parameter_names)

  # Do significance per bin for datacard
  print("\n\n****************Doing significance per bin****************")
  # Convert txt datacard to root datacard for prefit that requires masking
  # text2workspace.py doesn't seem to be able to detect errors
  runCommand('text2workspace.py '+args.datacard+' --channel-masks', None, printOutput, noRun)
  datacard_root = args.datacard.replace('.txt','.root')
  mask_parameter_names = set_mask_parameter_names()
  # postfit_significance_results[POI_parameter_name] = [limit]
  postfit_significance_results = {}
  for POI_parameter_name in POI_parameter_names:
    xy_name = '_'.join(POI_parameter_name.split('_')[1:3])
    plane_name = '_'.join(POI_parameter_name.split('_')[3:])
    abcd_name = xy_to_abcd(xy_name)
    bin_unmasks = []
    bin_unmasks.append('mask_'+abcd_to_xy("C")+"_"+plane_name+'=1')
    bin_unmasks.append('mask_'+abcd_to_xy("D")+"_"+plane_name+'=1')
    if abcd_name == "A0":
      bin_unmasks.append('mask_'+abcd_to_xy("A0")+"_"+plane_name+'=1')
      bin_unmasks.append('mask_'+abcd_to_xy("B0")+"_"+plane_name+'=1')
    elif abcd_name == "A1":
      bin_unmasks.append('mask_'+abcd_to_xy("A1")+"_"+plane_name+'=1')
      bin_unmasks.append('mask_'+abcd_to_xy("B1")+"_"+plane_name+'=1')
    plane_mask = list(set(mask_parameter_names).difference(bin_unmasks))

    postfit_plane_combine_command = 'combine -M Significance '\
    +' --setParameters r=0,'+','.join(plane_mask)\
    +' --freezeParameters r '+datacard_root+' --name _postfit_'+plane_name\
    +' --rMin -10 --uncapped=1' # To allow negative significance
    postfit_plane_log_filename = 'postfit_'+abcd_name+'_'+plane_name+'_log.txt'
    runCommand(postfit_plane_combine_command, postfit_plane_log_filename, printOutput, noRun)

    postfit_significance_results[POI_parameter_name] = [get_limit('higgsCombine_postfit_'+plane_name+'.Significance.mH120.root'),0,0]

  if not noRun:
    # Collect results
    # POI_results[POI_parameter_name][prefit/postfit/significance] = [fit value, down uncertainty (difference), diff up uncertainty (difference)]
    POI_results = {}
    for POI_parameter_name in POI_parameter_names:
      if not POI_parameter_name in POI_results:
        POI_results[POI_parameter_name] = {}
      POI_results[POI_parameter_name]['prefit'] = prefit_POI_results[POI_parameter_name]
      POI_results[POI_parameter_name]['postfit'] = postfit_POI_results[POI_parameter_name]
      POI_results[POI_parameter_name]['significance'] = postfit_significance_results[POI_parameter_name]

    # Print
    for POI_parameter_name in POI_parameter_names:
      print(POI_parameter_name)
      print('  prefit: '+str(prefit_POI_results[POI_parameter_name]))
      print('  postfit: '+str(postfit_POI_results[POI_parameter_name]))
      print('  significance: '+str(postfit_significance_results[POI_parameter_name]))

    with open(args.output_json, 'w') as output_json:
      json.dump(POI_results, output_json)

    # Example of reading JSON
    with open(args.output_json) as POI_results_file:
      POI_results = json.load(POI_results_file)
    print(POI_results)
