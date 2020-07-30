#!/usr/bin/env python
from __future__ import print_function

###### Script that compares the yields in old and new ntuples
from ROOT import TChain, TH1D
import os, sys, subprocess
import glob
import json
import string
import time
import math

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def findBaseSampleNames(folder):
  infiles = set() # to remove duplicates
  #for file in glob.glob(folder+'/*TTJ*.root'):
  for file in glob.glob(folder+'/*.root'):
    dataset_tag = file.split('/')[-1]
    dataset_tag = dataset_tag.split('__RunIISummer16NanoAODv5__')[0]
    dataset_tag = dataset_tag.split('__RunIIFall17NanoAODv5__')[0]
    dataset_tag = dataset_tag.split('__RunIIAutumn18NanoAODv5__')[0]
    dataset_tag = dataset_tag.split('_ext')[0]
    dataset_tag = dataset_tag.split('_Tune')[0]
    dataset_tag.replace('pico_','')
    infiles.add(dataset_tag)
  sortedfiles = list()
  for file in infiles:
    sortedfiles.append(file)
  sortedfiles = sorted(sortedfiles)

  return sortedfiles

# Setting folders
oldfolder    = '/net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/unskimmed/'
newfolder    = '/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/unskimmed/'

oweight = "weight"
nweight = "weight"

## Finding tags for each dataset
sortedfiles = findBaseSampleNames(newfolder)
print ('Will compare below samples')
print(*sortedfiles,sep="\n")

print ('\nOLD FOLDER: '+oldfolder)
print ('NEW FOLDER: '+newfolder)
print ('OLD WEIGHT "'+oweight+'"  - NEW WEIGHT "'+nweight)

print ('\n{:>40}'.format(' Ntuple name          ')+'{:>16}'.format('Difference')+'{:>17}'.format('Old yield'),end='')
print ('{:>17}'.format('New yield')+'{:>17}'.format('Old entries')+'{:>17}'.format('New entries'))
print ('=' * 128)
not_in_old = list()
not_in_new = list()
rows = list()
line = 1
histo = TH1D("histo","",10,0,10)
for ifile in sortedfiles:
  ochain = TChain("tree")
  oldntuples = oldfolder+"/*"+ifile+'*root'
  no = ochain.Add(oldntuples)
  if no == 0:
    not_in_old.append(ifile)
    continue;
  nchain = TChain("tree")
  newntuples = newfolder+"/*"+ifile+'*root'
  nn = nchain.Add(newntuples)
  if nn == 0:
    not_in_new.append(ifile)
    continue;

  no = ochain.Draw("1>>histo",oweight,"goff")
  oldtot = histo.Integral()
  nn = nchain.Draw("1>>histo",nweight,"goff")
  newtot = histo.Integral()
    
  if oldtot != 0 : diff = (newtot-oldtot)*100/oldtot
  elif newtot == 0 : diff = 0
  else :  diff = 999

  pretag = ""
  posttag = ""
  ## Appending rows with significant differences for later printing
  if abs(diff) > 150/math.sqrt(no+1) and abs(diff) > 150/math.sqrt(nn+1): 
    rows.append([ifile, abs(diff), oldtot, newtot, no, nn])
    pretag = bcolors.FAIL
    posttag = bcolors.ENDC
  ## Printing all rows
  print (pretag+'{:>40}'.format(ifile)+'{:>14.2f}'.format(diff)+' %'+'{:>17.2f}'.format(oldtot),end='')
  print ('{:>17.2f}'.format(newtot)+'{:>17}'.format(no)+'{:>17}'.format(nn)+posttag)
  if line == 5 : 
    print
    line = 0
  line += 1

## Sorting rows by difference
if len(rows) > 0: 
  print (bcolors.FAIL + "\nSamples off by more than 1.5 sigma"+ bcolors.ENDC)
  print ('\n{:>40}'.format(' Ntuple name          ')+'{:>16}'.format('Difference')+'{:>17}'.format('Old yield'),end='')
  print ('{:>17}'.format('New yield')+'{:>17}'.format('Old entries')+'{:>17}'.format('New entries'))
  print ('=' * 128)
rows = sorted(rows, key=lambda rows: rows[1])
line = 1
for row in rows:
  ifile = row[0]
  diff = row[1]
  oldtot = row[2]
  newtot = row[3]
  no = row[4]
  nn = row[5]

  ## Printing rows with significant differences
  print ('{:>40}'.format(ifile)+'{:>14.2f}'.format(diff)+' %'+'{:>17.2f}'.format(oldtot),end='')
  print ('{:>17.2f}'.format(newtot)+'{:>17}'.format(no)+'{:>17}'.format(nn))
  if line == 5 : 
    print
    line = 0
    line += 1

if len(not_in_old) > 0:
  print (bcolors.BOLD + '\nNtuples not found in '+oldfolder+':'+ bcolors.ENDC)
  for ntu in not_in_old:
    print ('\t'+ntu)
if len(not_in_new) > 0:
  print (bcolors.BOLD + '\nNtuples not found in '+newfolder+':'+ bcolors.ENDC)
  for ntu in not_in_new:
    print ('\t'+ntu)

print ()
sys.exit(0)
