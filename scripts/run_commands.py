#!/bin/env python
import os
import sys
import multiprocessing
import glob
import subprocess
from time import time

def runCommand(command):
  print(command)
  #print(command[1])
  os.system(command[1])

if __name__ == '__main__':
  t0 = time()

  cmsName = sys.argv[1]
  commands = subprocess.check_output(cmsName).decode('utf8').rstrip().split('\n')
  commands_info = [(index, commands[index]) for index in range(len(commands))]

  #pool = multiprocessing.Pool(processes=8)
  pool = multiprocessing.Pool()
  pool.map(runCommand, commands_info)

  print('Total commands: '+str(len(commands)))

  print('\nProgram took %.0fm %.0fs.' % ((time()-t0)/60,(time()-t0)%60))
