#!/bin/env python
import glob
import ROOT
import os, sys
import subprocess
import json
import argparse

# fileGlob: *.root
# searchString: run==1&&event==1000
def find_files(folder, fileGlob, treeName, searchString):
  print("[Info] Finding file")
  filenames = glob.glob(folder+"/"+fileGlob)
  print("Finding "+folder+"/"+fileGlob)
  found = []
  for filename in filenames:
    rootFile = ROOT.TFile(filename)
    tree_rootFile = rootFile.Get(treeName)
    if tree_rootFile:
      if (tree_rootFile.GetEntries(searchString) != 0):
        print("  "+filename)
        found.append(filename)
  return found

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='''Finds root files that pass cut.''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-f','--search_folder', required=True, help='Folder that contains root files to search between')
  parser.add_argument('-g','--glob', default = '*.root', help='Glob string to select searched root files by filename. Default is *.root.')
  parser.add_argument('-t','--tree_name', required=True, help = 'Tree name in file')
  parser.add_argument('-c','--cut_string', required=True, help = 'ROOT cut string')
  args = parser.parse_args()

  print(find_files(args.search_folder, args.glob, args.tree_name, args.cut_string))
