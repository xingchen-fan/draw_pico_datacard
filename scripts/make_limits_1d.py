#!/bin/env python
import os
import sys
import multiprocessing
import glob
import subprocess
from time import time

def runCommand(command, debug):
  print(command)
  if debug: raw_input("Press enter")
  #print(command)
  #os.system(command)
  #return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  out, err = process.communicate()
  return out.decode("utf-8")+'\n'+err.decode("utf-8")

if __name__ == '__main__':
  tag = sys.argv[1]
  debug = False
  if len(sys.argv) >= 3: debug = True

  outputFolder = 'cards/cards_1d_kappa_'+tag
  commandFilename = 'cards_cards_1d_kappa_'+tag+'_limits_cmds.py'
  resultsFilename = outputFolder+'/resolved_limits_kappa.txt'

  commands = [
  'scons && ./run/higgsino/write_kappa_datacards_run2.exe -o '+outputFolder+' -m CN -p "127_0,150_0,175_0,200_0,225_0,250_0,275_0,300_0,325_0,350_0,375_0,400_0,425_0,450_0,475_0,500_0,550_0,600_0,650_0,700_0,750_0,800_0,850_0,900_0,950_0,1000_0,1100_0,1200_0"',
  './scripts/write_combine_cmds.py -c '+outputFolder+' -m CN',
  'mv '+commandFilename+' '+commandFilename[:-2]+'json '+outputFolder,
  './scripts/run_commands.py '+outputFolder+'/'+commandFilename,
  'cat '+outputFolder+'/scan_point*mLSP-0/limit*txt | sort >> '+resultsFilename,
  'rp ./run/higgsino/plot_limit.exe -f '+resultsFilename+' --tag '+tag
  ]

  for command in commands:
    out= runCommand(command, debug)
    print (out)
