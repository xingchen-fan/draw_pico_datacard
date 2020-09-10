#!/usr/bin/env python
import os, math
import ROOT
ROOT.gROOT.SetBatch(True)

def getProfilesFromUrl(profile_url, histNames):
  profile_file = ROOT.TFile.Open(profile_url)
  profiles = []
  for histName in histNames:
    profile = profile_file.Get(histName)
    profile.SetDirectory(0)
    profiles.append(profile)
  return profiles, profile_file

def printWeights(profile_mc_url, profile_data_url, tag=""):
  profiles_mc, profile_mc_file = getProfilesFromUrl(profile_mc_url, ['pu_mc'])
  # Should plus and minus be used?
  profile_data_names = ['pileup', 'pileup_plus', 'pileup_minus']
  profiles_data, profile_data_file = getProfilesFromUrl(profile_data_url, profile_data_names)
  print("profile_mc bins: "+str(profiles_mc[0].GetNbinsX())+" profile_data bins: "+str(profiles_data[0].GetNbinsX()))

  for iProfile in xrange(len(profile_data_names)):
    profile_data_name = profile_data_names[iProfile]
    profile_data = profiles_data[iProfile]

    # Draw the histogram
    canvas_name = profile_data_name+"_data_vs_mc"
    canvas = ROOT.TCanvas(canvas_name, canvas_name,500,500)
    profile_data.DrawNormalized("hist")
    profile_mc = profiles_mc[0].DrawNormalized("same")
    profile_mc.SetLineColor(2)
    #profile_mc.Draw("hist same")
    print("Saving comparison histogram as "+canvas_name+"_"+tag+".pdf")
    canvas.SaveAs(canvas_name+"_"+tag+".pdf")

    # Calculate weight
    profile_data.Scale(1./profile_data.Integral())
    profile_mc = profiles_mc[0].DrawNormalized()
    profile_data.Divide(profile_mc)
    wgt = [profile_data.GetBinContent(i+1) for i in range(profile_data.GetNbinsX())]
    if profile_data_name=="pileup": 
        print ("Nominal weights:")
        for j,iwgt in enumerate(wgt):
            print "NPV: "+'{:>3d}'.format(j+1),
            print (" Weight: "+'{:>10.3e}'.format(iwgt))

    print ("------> Vector for weight_tools:")
    print "w_pu_"+profile_data_name,
    print (" = vector<double>({"+', '.join('{:.3e}'.format(x) for x in wgt)+"});")

    # Draw weight histogram
    canvas_name = profile_data_name+"_weight_"+tag
    canvas = ROOT.TCanvas(canvas_name, canvas_name,500,500)
    if tag=="2018": profile_data.SetMaximum(10)
    profile_data.SetMarkerStyle(8)
    profile_data.SetTitle("Pileup Weight")
    profile_data.Draw("hist p")
    canvas.SaveAs(canvas_name+".pdf")
  

if __name__ == "__main__":
  # TODO Doesn't like SL7 ROOT for some reason... Can get file from url
  # ----- NEEDED INPUTS --------------------------------------
  # MC profile and data profile for each year.
  # Data profile can be found using https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
  # A reference for MC profile can be found at https://gitlab.cern.ch/cp3-cms/bamboo/-/blob/master/bamboo/scripts/makePUReWeightJSON.py
  # NanoAOD provides root files for data and MC profiles at https://github.com/cms-nanoAOD/nanoAOD-tools/tree/019649a008c2a0756becbdff21811af9b0c0593c/python/postprocessing/data/pileup
  
  profile_mc_2016_url = 'https://github.com/cms-nanoAOD/nanoAOD-tools/raw/019649a008c2a0756becbdff21811af9b0c0593c/python/postprocessing/data/pileup/pileup_profile_Summer16.root'
  profile_data_2016_url = 'https://github.com/cms-nanoAOD/nanoAOD-tools/raw/019649a008c2a0756becbdff21811af9b0c0593c/python/postprocessing/data/pileup/PileupData_GoldenJSON_Full2016.root'
  profile_mc_2017_url = 'https://github.com/cms-nanoAOD/nanoAOD-tools/raw/019649a008c2a0756becbdff21811af9b0c0593c/python/postprocessing/data/pileup/mcPileup2017.root'
  profile_data_2017_url = 'https://github.com/cms-nanoAOD/nanoAOD-tools/raw/019649a008c2a0756becbdff21811af9b0c0593c/python/postprocessing/data/pileup/PileupHistogram-goldenJSON-13tev-2017-99bins_withVar.root'
  profile_mc_2018_url = 'https://github.com/cms-nanoAOD/nanoAOD-tools/raw/019649a008c2a0756becbdff21811af9b0c0593c/python/postprocessing/data/pileup/mcPileup2018.root'
  profile_data_2018_url = 'https://github.com/cms-nanoAOD/nanoAOD-tools/raw/019649a008c2a0756becbdff21811af9b0c0593c/python/postprocessing/data/pileup/PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root'
  

  printWeights(profile_mc_2016_url, profile_data_2016_url, "2016")
  printWeights(profile_mc_2017_url, profile_data_2017_url, "2017")
  printWeights(profile_mc_2018_url, profile_data_2018_url, "2018")


## ----- NEEDED INPUTS --------------------------------------
## Instructions:
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
#
#cmssw = os.getenv("CMSSW_BASE")+"/src/"
## retrieve latest json with per-lumiblock luminosity information
## /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
## fyi, how it gets derived: https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/660/1.html
#lumi_json = cmssw + "babymaker/data/json/pileup_latest.txt"
#
## obtain certified lumis JSON for dataset we are trying to reweight to
## /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/
#json = cmssw + "babymaker/data/json/golden_Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16.json"
#
## Minbias cross-section, (corrected) announcement for ICHEP16:
## https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/613/2/1/1/1.html
#mbxsec, mbxsec_relunc = 69200, 0.046
#
## get hist of the MC PU profile
#gROOT.SetBatch(kTRUE)
#mcfile = TChain("tree")
#mcfile.Add("/net/cms2/cms2r0/babymaker/babies/2017_01_21/mc/unprocessed/fullbaby_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_*")
#hmc = TH1D("hmc","hmc",75,0,75)
#mcfile.Draw("ntrupv_mean>>hmc","","norm")
#
## To reweight both the in- and out-of-time pile up
## pileupCalc should be used in mode "true" and the variable to reweight is npvtru_mean:
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookExercisesHATS_Pileup_2013
#calc_mode = "true" 
#
## ----------------------------------------------------------
#
#mbxsec_dict = {
#     "nom": mbxsec,
#      "up": mbxsec*(1+mbxsec_relunc),
#    "down": mbxsec*(1-mbxsec_relunc)
#}
#
#for imb in mbxsec_dict:
#    cmd  = "pileupCalc.py"
#    cmd += " -i "+json
#    cmd += " --inputLumiJSON "+lumi_json
#    cmd += " --calcMode "+calc_mode
#    cmd += " --minBiasXsec "+str(mbxsec_dict[imb])
#    cmd += " --maxPileupBin "+str(hmc.GetNbinsX())
#    cmd += " --numPileupBins "+str(hmc.GetNbinsX())
#    cmd += " pileup_"+imb+".root"
#
#    print "Obtaining data pile up distribution for variation:", imb
#    print cmd
#    os.system(cmd)
#
#    fdata = TFile("pileup_"+imb+".root","READ")
#    htmp = fdata.Get("pileup").Clone("pu_"+imb)
#    htmp.Scale(1./htmp.Integral())
#    htmp.Divide(hmc)
#    htmp.SetDirectory(0)
#    fdata.Close()
#    wgt = [htmp.GetBinContent(i+1) for i in range(htmp.GetNbinsX())]
#    if imb=="nom": 
#        print "Nominal weights:"
#        for j,iwgt in enumerate(wgt):
#            print "NPV: "+'{:>3d}'.format(j+1),
#            print " Weight: "+'{:>10.3e}'.format(iwgt)
#
#    print "------> Vector for weight_tools:"
#    print "w_pu_"+imb,
#    print " = vector<double>({"+', '.join('{:.3e}'.format(x) for x in wgt)+"});"
