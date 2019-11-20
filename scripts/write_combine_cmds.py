#!/usr/bin/env python3

import os, argparse, sys
from glob import glob

parser = argparse.ArgumentParser(description="Submits batch jobs to calculate Higgsino limits for full scan.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c","--card_dir", default= os.path.join(os.getcwd(),"cards"),
                    help="Directory where the datacards are")
args = parser.parse_args()

tag = args.card_dir.strip('/').replace('/','_')
card_paths = glob(os.path.join(args.card_dir,'*.txt'))
print("Found {} datacards.\n".format(len(card_paths)))

cleanup = glob(os.path.join(args.card_dir,'scan_point*'))
if (len(cleanup)):
  sys.exit("Cards Directory already contains previous combine output. Please clean up or create new input folder")

out_cmd_file = tag+'_limits_cmds.py'
cmdfile = open(out_cmd_file,'w')
cmdfile.write('#!/bin/env python\n')
for icard_path in card_paths:
  in_dir = os.path.dirname(icard_path)
  file_name = os.path.basename(icard_path)
  cmd = '{}/run/higgsino/scan_point.exe -i {} -f {}'.format(os.getcwd(), in_dir, file_name)
  cmdfile.write('print(\"'+cmd+'\")\n')

cmdfile.close()
os.chmod(out_cmd_file, 0o755)
os.system('convert_cl_to_jobs_info.py '+out_cmd_file+' '+ out_cmd_file.replace('.py','.json'))

print('Last line in '+out_cmd_file+' is:')
os.system('cat '+out_cmd_file+' | tail -n 1\n')

print("To submit jobs do: ")
print('auto_submit_jobs.py '+out_cmd_file.replace('.py','.json'))
  
