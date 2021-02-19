#!/bin/env python
import os,sys
from collections import OrderedDict
import re

def printHtWeights(htweight_filename, htweight_function_name):
  # Make htweight_dict from file
  # htweight_dict = {[ht_min, ht_max]: [data/mc, error]}
  htweight_dict = OrderedDict()
  with open(htweight_filename) as htweight_file:
    for line in htweight_file:
      binRange = [int(it) for it in re.findall('Bin.*?\d+,.*?\d+', line)[0][3:].split(',')]
      ratioAndError = [float(it) for it in re.findall('Ratio  =.*?\d+\.\d+.*?\d+\.\d+', line)[0][8:].split('+-')]
      htweight_dict[tuple(binRange)] = ratioAndError

  # Print htweight_dict
  nBin = len(htweight_dict)
  htweight_string = ''
  for iBin, (binRange, ratioAndError) in enumerate(htweight_dict.iteritems()):
    #print(iBin, binRange, ratioAndError)
    if iBin==0: 
      htweight_string += 'const NamedFunc '+htweight_function_name+'("'+htweight_function_name+'", [](const Baby &b) -> NamedFunc::ScalarType{\n'
      htweight_string += '  if (b.SampleType()<0) return 1.;\n'
      htweight_string += '  float weight = 1;\n'
      htweight_string += '  if (b.ht()>'+str(binRange[0])+' && b.ht()<='+str(binRange[1])+') weight = '+str(ratioAndError[0])+';\n'
    elif iBin==nBin-1:
      htweight_string += '  else if (b.ht()>'+str(binRange[0])+') weight = '+str(ratioAndError[0])+';\n'
    else:
      htweight_string += '  else if (b.ht()>'+str(binRange[0])+' && b.ht()<='+str(binRange[1])+') weight = '+str(ratioAndError[0])+';\n'
  htweight_string += '  return weight;\n'
  htweight_string += '});\n'

  print (htweight_string)

if __name__ == "__main__":
  year = '2017'

  #htweight_filename = "../txt/ht_"+year+".txt";
  #htweight_function_name = 'weight_ht_'+year
  #printHtWeights(htweight_filename, htweight_function_name)
  #htweight_filename = "../txt/ht_lowdrmax_"+year+".txt";
  #htweight_function_name = 'weight_ht_lowdrmax_'+year
  #printHtWeights(htweight_filename, htweight_function_name)
  #htweight_filename = "../txt/ht_highdrmax_"+year+".txt";
  #htweight_function_name = 'weight_ht_highdrmax_'+year
  #printHtWeights(htweight_filename, htweight_function_name)

  htweight_filename = "../txt/ht_sideband_2016.txt";
  htweight_function_name = 'weight_ht_sideband_2016'
  printHtWeights(htweight_filename, htweight_function_name)
  htweight_filename = "../txt/ht_sideband_2017.txt";
  htweight_function_name = 'weight_ht_sideband_2017'
  printHtWeights(htweight_filename, htweight_function_name)
  htweight_filename = "../txt/ht_sideband_2018.txt";
  htweight_function_name = 'weight_ht_sideband_2018'
  printHtWeights(htweight_filename, htweight_function_name)

