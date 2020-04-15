#!/usr/bin/env python

from ROOT import TChain, TH1D
from collections import OrderedDict 

ra2 = TChain("PreSelection")
ra2.Add("/net/cms25/cms25r5/jbkim/ra2b/ra2b_v18d/skim_met150/*root")

pico = TChain("tree")
pico.Add("/net/cms37/data1/pico/NanoAODv5/higgsino_ra2bsynch/2016/mc_skim150/raw_pico/*root")



cuts = OrderedDict()
cuts["MET>150"] = ["1","1"]
cuts["0 mus"] = ["nvmu==0", "Sum$(Muons.Pt()>10 && abs(Muons.Eta())<2.4 && Muons_mediumID && Muons_MiniIso<0.2)==0"]
cuts["0 els"] = ["nvel==0", "Sum$(Electrons.Pt()>10 && abs(Electrons.Eta())<2.5 && Electrons_MiniIso<0.1)==0"]
cuts["0 el trk"] = ["Sum$(abs(tk_pdgid)==11)==0", 
                    "Sum$(TAPElectronTracks.Pt()>5 && abs(TAPElectronTracks.Eta())<2.5 && abs(TAPElectronTracks_dxypv)<0.2 && abs(TAPElectronTracks_dzpv)<0.1 && ((TAPElectronTracks.Pt()<25 && TAPElectronTracks_pfRelIso03chg*TAPElectronTracks.Pt()<5) || TAPElectronTracks_pfRelIso03chg<0.2) && TAPElectronTracks_mT<100)==0"]
cuts["0 mu trk"] = ["Sum$(abs(tk_pdgid)==13)==0", 
                    "Sum$(TAPMuonTracks.Pt()>5 && abs(TAPMuonTracks.Eta())<2.5 && abs(TAPMuonTracks_dxypv)<0.2 && abs(TAPMuonTracks_dzpv)<0.1 && ((TAPMuonTracks.Pt()<25 && TAPMuonTracks_pfRelIso03chg*TAPMuonTracks.Pt()<5) || TAPMuonTracks_pfRelIso03chg<0.2) && TAPMuonTracks_mT<100)==0"]
cuts["0 pi trk"] = ["Sum$(abs(tk_pdgid)==211)==0", 
                    "Sum$(TAPPionTracks.Pt()>10 && abs(TAPPionTracks.Eta())<2.5 && abs(TAPPionTracks_dxypv)<0.2 && abs(TAPPionTracks_dzpv)<0.1 && ((TAPPionTracks.Pt()<25 && TAPPionTracks_pfRelIso03chg*TAPPionTracks.Pt()<5) || TAPPionTracks_pfRelIso03chg<0.1) && TAPPionTracks_mT<100)==0"]
cuts["4-5 jets"] = ["(njet==4 || njet==5)","(NJets==4 || NJets==5)"]
cuts["2b, tight"] = ["nbt>=2", "Sum$(Jets.Pt()>30 && abs(Jets.Eta())<=2.4 && Jets_bJetTagDeepCSVBvsAll>0.8953)>=2"]
cuts["low dPhi"] = ["!low_dphi_mht","DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3"]


print("{:<15} | {:>15} | {:>15} | {:>10}".format("Cut","Pico","RA2b", "Ratio"))
print('-'*70)

full_cut_pico, full_cut_ra2 = "1", "1"
for icut in cuts.keys():
  full_cut_pico = full_cut_pico + "&&" + cuts[icut][0]
  full_cut_ra2 = full_cut_ra2 + "&&" + cuts[icut][1]
  npico = pico.GetEntries(full_cut_pico)
  nra2 = ra2.GetEntries(full_cut_ra2)
  ratio = float(npico)/float(nra2)
  print("{:<15} | {:>15d} | {:>15d} | {:>10.3f}".format(icut, npico, nra2, ratio))

