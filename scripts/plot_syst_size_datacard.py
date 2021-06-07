#!/usr/bin/env python
#utility for plotting systematic from MC statistics given datacards
import ROOT
import os.path
import sys
import math
import argparse

def get_average_systematic(prod_mass, lsp_mass, folder_name, syst_name):
  '''function that returns the average value of the systematic syst_name given 
  a datacard via produced sparticle mass, LSP mass, and a folder
  '''
  datacard_filename = folder_name+'/datacard-TChiHH_mChi-'+chi1_mass+'_mLSP-'+chi0_mass+'_Tune_2016,2017,2018_priority1_resolved.txt'
  if (not os.path.isfile(datacard_filename)):
    datacard_filename = folder_name+'/datacard-T5HH_mGluino-'+chi1_mass+'_mLSP-'+chi0_mass+'_Tune_2016,2017,2018_priority1_resolved.txt'
  if (not os.path.isfile(datacard_filename)):
    #skip this mass point
    return 0
  
  datacard_file = open(datacard_filename,'r')
  datacard_lines = datacard_file.read().split('\n')
  bin_yields = []
  syst_tuples = []
  total_events = 0
  weighted_stat_syst = 0
  
  for line in datacard_lines:
    entries = line.split()
    if (len(entries) > 0):
      if (entries[0] == 'rate'):
        #yield line
        entry_number = 1
        while entry_number < len(entries):
          bin_yields.append(float(entries[entry_number]))
          total_events += float(entries[entry_number])
          entry_number += 2
        if (total_events == 0):
          print('ERROR: no events')
          return

      elif (entries[1] == 'lnN'):

        if (entries[0][0:4] == 'stat'):
          #signal systematic from MC statistics
          entry_number = 2
          entry_index = 0
          while entry_number < len(entries):
            if entries[entry_number] == '-':
              entry_number += 2
              entry_index += 1
              continue
            syst_value = abs(float(entries[entry_number])-1)
            weighted_stat_syst += syst_value*bin_yields[entry_index]
            entry_number += 2
            entry_index += 1

        elif (entries[0][0:2] != 'cr'):
          #signal systematic uncertainty that is not MC statistics
          weighted_systematic = 0
          entry_number = 2
          entry_index = 0
          while entry_number < len(entries):
            syst_value = 0
            if '/' in entries[entry_number]:
              #asymmetric uncertainties
              down_syst = float(entries[entry_number].split('/')[0])
              up_syst = float(entries[entry_number].split('/')[1])
              down_syst = abs(down_syst-1.0)
              up_syst = abs(up_syst-1.0)
              syst_value = (down_syst+up_syst)/2.0
            else:
              syst_value = abs(float(entries[entry_number])-1)
            weighted_systematic += syst_value*bin_yields[entry_index]
            entry_number += 2
            entry_index += 1
          syst_tuples.append((entries[0], weighted_systematic/total_events*100))

  #sort systmatics by size
  syst_tuples.append(('stat', weighted_stat_syst/total_events*100))
  syst_tuples = sorted(syst_tuples, key=lambda syst: syst[1])
  return_value = 0
  for syst in syst_tuples:
    #print('Average systematic from {0} is {1}'.format(syst[0], syst[1]))
    if syst[0] == syst_name:
      return_value = syst[1]
  return return_value


if __name__ == '__main__':
  #parse arguments
  parser = argparse.ArgumentParser(description=('generates a histogram of average '
                                   'uncertainty from MC statistics given a '
                                   'collection of datacards'),
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-d","--directory", default= ".",
                      help="Directory where the datacards are")
  parser.add_argument("-m","--model", default="TChiHH2D",
                      help="TChiHH2D, TChiHH1D, T5HH2D, or T5HH1D")
  args = parser.parse_args()

  #fill histogram using datacards
  stat_hist = ROOT.TH2D('stat_hist',
      'Uncertainty from MC statistics [%]; m_{#chi_{1}} [GeV]; m_{#chi_{0}} [GeV]',
      33,0,825,29,0,725)

  if args.model == "TChiHH2D":
    chi1_mass_int = 127
    while chi1_mass_int <= 800:
      chi1_mass = str(chi1_mass_int)
      chi0_mass_int = 0
      while chi0_mass_int <= chi1_mass_int-125:
        #print('Analyzing signal point ({0},{1})'.format(chi1_mass_int, chi0_mass_int))
        chi0_mass = str(chi0_mass_int)
        stat_hist_size = get_average_systematic(chi1_mass, chi0_mass, args.directory, 'stat')
        #technically should use SetBinContent
        #modify for new stats
        if chi0_mass_int == 0:
            stat_hist_size *= math.sqrt(100000.0)
        else:
            stat_hist_size *= math.sqrt(44000.0)
        #
        if ((chi1_mass_int < 400) and (chi0_mass_int < 150)):
            stat_hist_size *= 1.0/math.sqrt(84000.0)
        elif (chi0_mass_int >= (chi1_mass_int-250)):
            stat_hist_size *= 1.0/math.sqrt(84000.0)
        else:
            stat_hist_size *= 1.0/math.sqrt(42000.0)
        stat_hist.Fill(chi1_mass_int, chi0_mass_int, stat_hist_size)
        chi0_mass_int += 25
      if chi1_mass_int == 127:
        chi1_mass_int = 150
      else:
        chi1_mass_int += 25
  
  #draw histogram
  stat_hist_smooth = stat_hist.Clone('stat_hist_smooth')
  for xbin in range(stat_hist.GetNbinsX()+1):
    for ybin in range(stat_hist.GetNbinsY()):
      if (xbin > (ybin+4)):
        if (stat_hist_smooth.GetBinContent(xbin,ybin)==0):
          #average bins around it
          tot_around = 0
          n_around = 0
          around_coords = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
          for around_coord in around_coords:
            single_around = stat_hist.GetBinContent(xbin+around_coord[0],ybin+around_coord[1])
            if single_around > 0:
              n_around += 1.0
            tot_around += single_around
          if n_around > 0:
            tot_around = tot_around/n_around
            stat_hist_smooth.SetBinContent(xbin,ybin,tot_around)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetPaintTextFormat('3.1f')
  canvas = ROOT.TCanvas('canvas','canvas',1000,1000)
  canvas.SetRightMargin(0.15)
  canvas.SetLeftMargin(0.15)
  canvas.SetTopMargin(0.15)
  canvas.SetBottomMargin(0.15)
  stat_hist_smooth.Draw('colz')
  stat_hist.SetMarkerSize(0.4)
  stat_hist.Draw('text same')
  canvas.Draw()
  canvas.SaveAs('stat_syst.pdf')
  
  #if (dim == '1d'):
  #  chi1_mass_int = 127
  #  while chi1_mass_int < 1500:
  #    chi1_mass = str(chi1_mass_int)
  #    print('Analyzing signal point ({0},0)'.format(chi1_mass_int))
  #    avr_ret = get_average_systematics(chi1_mass, '0', datacard_folder, syst_dict)
  #    syst_dict = avr_ret[0]
  #    total_entries += avr_ret[1]
  #    if chi1_mass_int == 127:
  #      chi1_mass_int = 150
  #    else:
  #      chi1_mass_int += 25
  #elif (dim == '2d'):
  #  chi1_mass_int = 127
  #  while chi1_mass_int <= 800:
  #    chi1_mass = str(chi1_mass_int)
  #    chi0_mass_int = 0
  #    while chi0_mass_int <= chi1_mass_int-125:
  #      print('Analyzing signal point ({0},{1})'.format(chi1_mass_int, chi0_mass_int))
  #      chi0_mass = str(chi0_mass_int)
  #      avr_ret = get_average_systematics(chi1_mass, chi0_mass, datacard_folder, syst_dict)
  #      syst_dict = avr_ret[0]
  #      total_entries += avr_ret[1]
  #      chi0_mass_int += 25
  #    if chi1_mass_int == 127:
  #      chi1_mass_int = 150
  #    else:
  #      chi1_mass_int += 25
  #elif (dim == 't5hh'):
  #  gluino_mass_int = 1000
  #  while gluino_mass_int <= 2600:
  #    gluino_mass = str(gluino_mass_int)
  #    chi0_mass_int = 1
  #    while chi0_mass_int <= gluino_mass_int:
  #      print('Analyzing signal point ({0},{1})'.format(gluino_mass_int, chi0_mass_int))
  #      chi0_mass = str(chi0_mass_int)
  #      avr_ret = get_average_systematics(gluino_mass, chi0_mass, datacard_folder, syst_dict, True)
  #      syst_dict = avr_ret[0]
  #      total_entries += avr_ret[1]
  #      if chi0_mass_int == 1:
  #        chi0_mass_int = 50
  #      else:
  #        chi0_mass_int += 50
  #    gluino_mass_int += 50

