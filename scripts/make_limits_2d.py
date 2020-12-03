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

  outputFolder = 'cards/cards_2d_kappa_'+tag
  commandFilename = 'cards_cards_2d_kappa_'+tag+'_limits_cmds.py'
  resultsFilename = outputFolder+'/resolved_limits_kappa.txt'

  commands = [
  'scons && ./run/higgsino/write_kappa_datacards_run2.exe -o '+outputFolder+' -m N1N2',
  './scripts/write_combine_cmds.py -c '+outputFolder+' -m N1N2',
  'mv '+commandFilename+' '+commandFilename[:-2]+'json '+outputFolder,
  './scripts/run_commands.py '+outputFolder+'/'+commandFilename,
  'cat '+outputFolder+'/scan_point*/limit*txt | sort >> '+resultsFilename,
  'rp ./run/higgsino/limit_scan.exe -f '+resultsFilename+' -m N1N2 --tag '+tag
  ]

  for command in commands:
    out= runCommand(command, debug)
    print (out)
