#!/usr/bin/env python
import glob
import re
import ROOT
import json
import array

def get_mChi_mLSP(filename):
  mChi_string = re.findall("mChi-[0-9]*", filename)
  mChi = int(mChi_string[0][5:])
  mLSP_string = re.findall("mLSP-[0-9]*", filename)
  mLSP = int(mLSP_string[0][5:])
  return mChi, mLSP

def getEntries(filename, cut_string):
  chain = ROOT.TChain("tree")
  chain.Add(filename)
  return chain.GetEntries(cut_string)

def getSignalEntries(pico_filename):
  ''' returns a list with properties from the file and yields in the format
  [mChi, mLSP, all entries, baseline, 4b, 4b-lowDrmax-highMet, 4b-highDrmax-lowMET]
  '''
  baseline_cut = ("met/mht<2 && met/met_calo<2&&weight<1.5&&"
                 "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                 "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                 "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))")
  signal_cut = ("met/mht<2 && met/met_calo<2&&weight<1.5&&"
                 "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                 "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                 "(nbt>=2&&nbm>=3&&nbl>=4)")
  high_signal_cut = ("met/mht<2 && met/met_calo<2&&weight<1.5&&"
                 "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                 "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                 "(nbt>=2&&nbm>=3&&nbl>=4)&&met>400&&hig_cand_drmax[0]<1.1")
  low_signal_cut = ("met/mht<2 && met/met_calo<2&&weight<1.5&&"
                 "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                 "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                 "(nbt>=2&&nbm>=3&&nbl>=4)&&met<200&&hig_cand_drmax[0]>=1.1")
  mChi = get_mChi_mLSP(pico_filename)[0]
  mLSP = get_mChi_mLSP(pico_filename)[1]
  allEntries = getEntries(pico_filename,"1")
  baselineEntries = getEntries(pico_filename, baseline_cut)
  signal_4b = getEntries(pico_filename, signal_cut)
  signal_4b_lowdrmax_highMET = getEntries(pico_filename, high_signal_cut)
  signal_4b_highdrmax_lowMET = getEntries(pico_filename, low_signal_cut)
  return [mChi, mLSP, allEntries, baselineEntries, signal_4b, signal_4b_lowdrmax_highMET, signal_4b_highdrmax_lowMET]

def get_signal_entries_v2(pico_filename):
  ''' returns a list with properties from the file and yields in the format
  [mChi, mLSP, all entries, SR]
  '''
  baseline_cut = ("met/mht<2 && met/met_calo<2&&weight<1.5&&"
                 "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                 "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                 "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))")
  signal_cut = ("met/mht<2 && met/met_calo<2&&weight<1.5&&"
                 "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                 "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                 "(nbt>=2&&nbm>=3)&&hig_camd_am[0]>100&&hig_cand_am[0]<140")
  mChi = get_mChi_mLSP(pico_filename)[0]
  mLSP = get_mChi_mLSP(pico_filename)[1]
  chain = ROOT.TChain("tree")
  chain.Add(pico_filename)
  allEntries = chain.GetEntries("1")
  #baselineEntries = chain.GetEntries(baseline_cut)
  signalregionEntries = chain.GetEntries(signal_cut)
  return [mChi,mLSP,allEntries,signalregionEntries]

def printSignalInfo(signalInfo):
  print ("mChi: "+str(signalInfo[0])+" mLSP: "+str(signalInfo[1]))
  print ("All entries:          "+str(signalInfo[2]))
  print ("Baseline:             "+str(signalInfo[3]))
  print ("4b_signal:            "+str(signalInfo[4]))
  print ("4b_lowDrmax_largeMet: "+str(signalInfo[5]))

def makeData(pico_foldername, outputSignalInfoFilename):
  allSignalInfo = []

  pico_filenames = glob.glob(pico_foldername+"/*.root")
  #print(pico_files)

  # Making data
  for pico_filename in pico_filenames:
    signalInfo = get_signal_entries_v2(pico_filename)
    print('Processing mChi: '+str(signalInfo[0])+' mLSP: '+str(signalInfo[1]))
    #printSignalInfo(signalInfo)
    allSignalInfo.append(signalInfo)
  # Save data
  with open(outputSignalInfoFilename, 'w') as signalInfoFile:
    json.dump(allSignalInfo, signalInfoFile)
  print("Saved to "+outputSignalInfoFilename)

  # Laod data
  with open(outputSignalInfoFilename) as signalInfoFile:
    allSignalInfo = json.load(signalInfoFile)
  for signalInfo in allSignalInfo:
    printSignalInfo(signalInfo)

def makeLowEdgesFromZero(values):
  sorted_values = sorted(values)
  #print(sorted_values)
  lowEdges = []
  for index in range(len(sorted_values)-1):
    highDiff = sorted_values[index+1] - sorted_values[index]
    if index == 0:
      lowEdge = sorted_values[index]-highDiff*1./2
      if lowEdge<0: lowEdge = 0
      lowEdges.append(lowEdge)
    lowEdges.append(sorted_values[index]+highDiff*1./2)
    if index == len(sorted_values)-2:
      lowEdges.append(sorted_values[index+1]+highDiff*1./2)
  return array.array('d',lowEdges)

def getBinIndex(value, binValues):
  index = 1
  for binValue in sorted(binValues):
    if (value - binValue)==0: break
    else: index += 1
  return index

def drawHistogram(hist_title, hist_year, hist_signalInfoIndex, allSignalInfoPerYear):
  # Find all the mChi, mLSP values
  mChiSet = set()
  mLSPSet = set()
  for year in allSignalInfoPerYear:
    for signalInfo in allSignalInfoPerYear[year]:
      # [mChi, mLSP, all entries, baseline, 4b,  4b-lowDrmax-highMet]
      mChiSet.add(signalInfo[0])
      mLSPSet.add(signalInfo[1])
  print(sorted(mChiSet))
  print(sorted(mLSPSet))
   
  mChiLowEdges = makeLowEdgesFromZero(mChiSet)
  mLSPLowEdges = makeLowEdgesFromZero(mLSPSet)
  #print(mChiLowEdges)
  #print(mLSPLowEdges)

  # x axis: 100 to 1500 in 25 ns steps for mChi
  # y axis: 0 to 1100 in 25 ns steps for mLSP
  histTotalEntries = ROOT.TH2D(hist_title, hist_title, len(mChiLowEdges)-1, mChiLowEdges, len(mLSPLowEdges)-1, mLSPLowEdges)
  partialHistTotalEntries = ROOT.TH2D("partial_"+hist_title, "partial_"+hist_title, len(mChiLowEdges)-1, mChiLowEdges, len(mLSPLowEdges)-1, mLSPLowEdges)
  for signalInfo in allSignalInfoPerYear[hist_year]:
    # [mChi, mLSP, all entries, baseline, 4b,  4b-lowDrmax-highMet]
    imChi = getBinIndex(signalInfo[0], sorted(mChiSet))
    imLSP = getBinIndex(signalInfo[1], sorted(mLSPSet))
    entries = signalInfo[hist_signalInfoIndex]
    histTotalEntries.SetBinContent(imChi,imLSP,entries)
    # Ignore some points
    if signalInfo[0]%100!=0: continue
    if signalInfo[1]%50!=0: continue
    partialHistTotalEntries.SetBinContent(imChi,imLSP,entries)
  c1 = ROOT.TCanvas("c1","c1",500,500)
  c1.SetRightMargin(0.15)
  histTotalEntries.Draw("colz")
  partialHistTotalEntries.SetMarkerSize(0.8)
  partialHistTotalEntries.Draw("TEXT same")
  c1.SaveAs(histTotalEntries.GetName()+".pdf")
  return histTotalEntries

def make_efficiency_plot(denom_hist, num_hist, str_year):
  eff_hist = num_hist.Clone('eff_hist')
  eff_hist.Divide(denom_hist)
  eff_hist_text = eff_hist.Clone('eff_hist_text')
  for xbin in range(eff_hist.GetNbinsX()):
    for ybin in range(eff_hist.GetNbinsY()):
      if (xbin > ybin):
        if (eff_hist.GetBinContent(xbin,ybin)==0):
          #average bins around it
          tot_around = 0
          n_around = 0
          around_coords = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
          for around_coord in around_coords:
            single_around = eff_hist_text.GetBinContent(xbin+around_coord[0],ybin+around_coord[1])
            if single_around > 0:
              n_around += 1.0
            tot_around += single_around
          if n_around > 0:
            tot_around = tot_around/n_around
            eff_hist.SetBinContent(xbin,ybin,tot_around)
      if ((xbin % 2 == 1) or (ybin % 2 == 0)):
        eff_hist_text.SetBinContent(xbin,ybin,0)
  eff_hist_text
  ROOT.gStyle.SetPaintTextFormat('1.4f');
  eff_hist.SetTitle('Baseline Cut Efficiency '+str_year+'; m_{#chi_{1}}; m_{#chi_{0}}')
  c1 = ROOT.TCanvas("c1","c1",1000,1000)
  c1.SetLogz()
  c1.SetRightMargin(0.15)
  c1.SetLeftMargin(0.15)
  c1.SetTopMargin(0.15)
  c1.SetBottomMargin(0.15)
  eff_hist.Draw("colz")
  eff_hist_text.SetMarkerSize(0.4)
  eff_hist_text.Draw('text same')
  c1.SaveAs('TChiHH2D_eff_'+str_year+'.pdf')

if __name__ == "__main__":

  # Make data
  #makeData("/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath/2016/SMS-TChiHH_2D_fastSimJmeCorrection/raw_pico","signalInfo_2016.txt")
  #makeData("/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath/2017/SMS-TChiHH_2D_fastSimJmeCorrection/raw_pico","signalInfo_2017.txt")
  #makeData("/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath/2018/SMS-TChiHH_2D_fastSimJmeCorrection/raw_pico","signalInfo_2018.txt")

  # Collect all data
  allSignalInfoPerYear = {}
  for year in [2016, 2017, 2018]:
    outputSignalInfoFilename = "signalInfo_"+str(year)+".txt"
    with open(outputSignalInfoFilename) as signalInfoFile:
      # [mChi, mLSP, all entries, baseline, 4b, 4b-lowDrmax-highMe, 4b-highDrmax-lowMETt]
      # [mChi, mLSP, all entries, SR]
      allSignalInfo = json.load(signalInfoFile)
      index1D = 0;
      index2D = 0;
      for point in allSignalInfo:
        if point[1] == 0:
          print("index: "+str(index1D)+" Year: "+str(year)+" mChi: "+str(point[0])+" mLSP: "+str(point[1]))
          index1D += 1
        else:
          print("index: "+str(index2D)+" Year: "+str(year)+" mChi: "+str(point[0])+" mLSP: "+str(point[1]))
          index2D += 1
    allSignalInfoPerYear[year] = allSignalInfo

  # Draw data
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(False)
  # [mChi, mLSP, all entries, baseline, 4b, 4b-lowDrmax-highMet, 4b-highDrmax-lowMET]
  # [mChi, mLSP, all entries, SR]
  for year in [2016, 2017, 2018]:
    denom_hist = drawHistogram("total_Entries_"+str(year), year, 2, allSignalInfoPerYear)
    num_hist = drawHistogram("signal_Entries_"+str(year), year, 3, allSignalInfoPerYear)
    make_efficiency_plot(denom_hist,num_hist,str(year))
  
  #years = ['2016','2017','2018']
  #for year in years:
  #  pico_foldername = '/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath/'+year+'/SMS-TChiHH_2D_fastSimJmeCorrection/unskimmed'
  #  pico_filenames = glob.glob(pico_foldername+"/*.root")

  #  # Making data
  #  for pico_filename in pico_filenames:
  #    signal_info = get_signal_entries_v2(pico_filename)
  #    print('Processing '+year+', mChi: '+str(signal_info[0])+' mLSP: '+str(signal_info[1]))
  #    #printSignalInfo(signalInfo)
