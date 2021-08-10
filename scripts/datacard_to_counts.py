#!/usr/bin/env python
import ROOT
import math

# iPlane, abcd
def parse_name(full_bin_name):
  abcd_name = 'unkown'
  iPlane = -1
  if ('xbkg' in full_bin_name and 'ybkg' in full_bin_name): abcd_name = 'D'
  elif ('xbkg' in full_bin_name and 'ysig' in full_bin_name): abcd_name = 'C'
  elif ('xsig0' in full_bin_name and 'ybkg' in full_bin_name): abcd_name = 'B0'
  elif ('xsig0' in full_bin_name and 'ysig' in full_bin_name): abcd_name = 'A0'
  elif ('xsig1' in full_bin_name and 'ybkg' in full_bin_name): abcd_name = 'B1'
  elif ('xsig1' in full_bin_name and 'ysig' in full_bin_name): abcd_name = 'A1'
  else: print("Unknown bin: "+full_bin_name)
  if ('met0' in full_bin_name and 'drmax0' in full_bin_name): iPlane = 0
  elif ('met0' in full_bin_name and 'drmax1' in full_bin_name): iPlane = 1
  elif ('met1' in full_bin_name and 'drmax0' in full_bin_name): iPlane = 2
  elif ('met1' in full_bin_name and 'drmax1' in full_bin_name): iPlane = 3
  elif ('met2' in full_bin_name and 'drmax0' in full_bin_name): iPlane = 4
  elif ('met2' in full_bin_name and 'drmax1' in full_bin_name): iPlane = 5
  elif ('met3' in full_bin_name and 'drmax0' in full_bin_name): iPlane = 6
  elif ('met3' in full_bin_name and 'drmax1' in full_bin_name): iPlane = 7
  else: print("Unknown plane: "+full_bin_name)
  return iPlane, abcd_name

def addBin(binA, binB):
  if len(binA) == 0: binA=[0,0]
  if len(binB) == 0: binB=[0,0]
  binA_value = binA[0]
  binA_uncertainty = binA[1]
  binB_value = binB[0]
  binB_uncertainty = binB[1]
  sum_value = binA_value+binB_value
  sum_uncertainty = math.sqrt(binA_uncertainty**2+binB_uncertainty**2)
  return [sum_value, sum_uncertainty]

if __name__ == "__main__":

  # background_only, signal
  model = "signal"

  # counts[iPlane][iABCD][abcd_index] = [value, error]
  counts = {}

  # From combine -M FitDiagnostics --saveNormalizations --saveWithUncertainties --numToysForShapes 1000 datacard.txt
  #filename = "abcd_card_study/full_bin_study/asymptotic/fitDiagnostics.root"
  filename = "abcd_card_study/full_bin_study/fit_450/fitDiagnostics.root"
  fit_file = ROOT.TFile(filename)

  if model == "background_only":
    norm_fit_b = fit_file.Get("norm_fit_b")
    it_norm_fit_b = norm_fit_b.createIterator()
    var_norm_fit_b = it_norm_fit_b.Next()
    while var_norm_fit_b:
      variable_name = var_norm_fit_b.GetName()
      variable_value = var_norm_fit_b.getValV()
      variable_error = var_norm_fit_b.getError()
      var_norm_fit_b = it_norm_fit_b.Next()
      if ('/sig' in variable_name): continue
      #print(variable_name, parse_name(variable_name), variable_value, variable_error)
      iPlane, abcd_name = parse_name(variable_name)
      if (iPlane not in counts): counts[iPlane] = [[[],[],[],[]],[[],[],[],[]]]
      if abcd_name == "D": 
        counts[iPlane][0][3] = [variable_value, variable_error]
        counts[iPlane][1][3] = [variable_value, variable_error]
      elif abcd_name == "C":
        counts[iPlane][0][2] = [variable_value, variable_error]
        counts[iPlane][1][2] = [variable_value, variable_error]
      elif abcd_name == "B0":
        counts[iPlane][0][1] = [variable_value, variable_error]
      elif abcd_name == "A0":
        counts[iPlane][0][0] = [variable_value, variable_error]
      elif abcd_name == "B1":
        counts[iPlane][1][1] = [variable_value, variable_error]
      elif abcd_name == "A1":
        counts[iPlane][1][0] = [variable_value, variable_error]

    counts_string = ''
    for iPlane in counts:
      for iABCD in range(2):
        counts_string += "iplane: "+str(iPlane)+" ibin: "+str(iABCD)+"\n"
        counts_string += "fit value: "+str(counts[iPlane][iABCD][3][0])+" "+str(counts[iPlane][iABCD][1][0])+" "+str(counts[iPlane][iABCD][2][0])+" "+str(counts[iPlane][iABCD][0][0])+"\n"
        counts_string += "fit unc: "+str(counts[iPlane][iABCD][3][1])+" "+str(counts[iPlane][iABCD][1][1])+" "+str(counts[iPlane][iABCD][2][1])+" "+str(counts[iPlane][iABCD][0][1])+"\n"
    print(counts_string)

  if model == "signal":
    norm_fit_s = fit_file.Get("norm_fit_s")
    it_norm_fit_s = norm_fit_s.createIterator()
    var_norm_fit_s = it_norm_fit_s.Next()
    while var_norm_fit_s:
      variable_name = var_norm_fit_s.GetName()
      variable_value = var_norm_fit_s.getValV()
      variable_error = var_norm_fit_s.getError()
      var_norm_fit_s = it_norm_fit_s.Next()
      #print(variable_name, parse_name(variable_name), variable_value, variable_error)
      iPlane, abcd_name = parse_name(variable_name)
      if (iPlane not in counts): counts[iPlane] = [[[],[],[],[]],[[],[],[],[]]]
      if abcd_name == "D": 
        counts[iPlane][0][3] = addBin([variable_value, variable_error], counts[iPlane][0][3])
        counts[iPlane][1][3] = addBin([variable_value, variable_error], counts[iPlane][1][3])
      elif abcd_name == "C":
        counts[iPlane][0][2] = addBin([variable_value, variable_error], counts[iPlane][0][2])
        counts[iPlane][1][2] = addBin([variable_value, variable_error], counts[iPlane][1][2])
      elif abcd_name == "B0":
        counts[iPlane][0][1] = addBin([variable_value, variable_error], counts[iPlane][0][1])
      elif abcd_name == "A0":
        counts[iPlane][0][0] = addBin([variable_value, variable_error], counts[iPlane][0][0])
      elif abcd_name == "B1":
        counts[iPlane][1][1] = addBin([variable_value, variable_error], counts[iPlane][1][1])
      elif abcd_name == "A1":
        counts[iPlane][1][0] = addBin([variable_value, variable_error], counts[iPlane][1][0])

    counts_string = ''
    for iPlane in counts:
      for iABCD in range(2):
        counts_string += "iplane: "+str(iPlane)+" ibin: "+str(iABCD)+"\n"
        counts_string += "fit value: "+str(counts[iPlane][iABCD][3][0])+" "+str(counts[iPlane][iABCD][1][0])+" "+str(counts[iPlane][iABCD][2][0])+" "+str(counts[iPlane][iABCD][0][0])+"\n"
        counts_string += "fit unc: "+str(counts[iPlane][iABCD][3][1])+" "+str(counts[iPlane][iABCD][1][1])+" "+str(counts[iPlane][iABCD][2][1])+" "+str(counts[iPlane][iABCD][0][1])+"\n"
    print(counts_string)
