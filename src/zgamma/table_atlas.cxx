#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TObject.h"
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

NamedFunc mWindow(float width) {
  double mh = 125.09;
  double max(mh+width/2), min(mh-width/2);
  TString cut = "llphoton_m[0] > " + to_string(min) + "&& llphoton_m[0] <" + to_string(max);
  NamedFunc out(cut);
  return out;
}

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  string bfolder("/net/cms29/cms29r0/pico/NanoAODv5/zgamma_channelIslandsv3/2016/");
  string mc_path( bfolder+"mc/merged_zgmc_llg/");
  string sig_path(bfolder+"HToZG/merged_zgmc_llg/");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  NamedFunc trigs(el_trigs || mu_trigs);
  auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma",       back, 
                       TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG*"}, trigs);
  auto proc_ewkzg = Process::MakeShared<Baby_pico>("EWK Z#gamma",      back, 
                       TColor::GetColor("#39a9df"),{mc_path+"*LLAJJ*"},   trigs);
  auto proc_dy    = Process::MakeShared<Baby_pico>("DY",               back, 
                       TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs && "stitch_dy");
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZ#gamma", sig, 
                       kRed     ,{sig_path+"*.root"},   trigs);
  auto proc_hzg_vbf = Process::MakeShared<Baby_pico>("HToZ#gamma VBF", sig, 
                       kMagenta                   ,{sig_path+"*VBF*.root"},   trigs);
  auto proc_hzg_gg  = Process::MakeShared<Baby_pico>("HToZ#gamma ggF", sig, 
                       TColor::GetColor("#ff0000"),{sig_path+"*GluGlu*.root"},   trigs);
  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_ewkzg, proc_hzg};
  vector<shared_ptr<Process>> theory_procs = {proc_dy, proc_smzg, proc_hzg};
  NamedFunc pTt("pTt",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 g = AssignGamma(b).Vect();
    TVector3 h = AssignH(b).Vect();
    TVector3 z = AssignZ(b).Vect();
    g.SetZ(0); h.SetZ(0); z.SetZ(0);
    return h.Cross((z-g).Unit()).Mag();
  });
  NamedFunc el_check("el_check",[](const Baby &b) -> NamedFunc::ScalarType{
    if(abs(b.el_eta()->at(b.ll_i1()->at(0)) - b.el_eta()->at(b.ll_i2()->at(0))) < 0.075 &&
       abs(b.el_phi()->at(b.ll_i1()->at(0)) - b.el_phi()->at(b.ll_i2()->at(0))) < 0.125)
       return false;
    return true;
  });
  NamedFunc gen_relpT("gen_relpT",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector h = AssignH(b,true);
    TLorentzVector g = AssignGamma(b,true);
    return g.Pt()/h.M();
  });
  NamedFunc relpT("relpT",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.nphoton() == 1) return b.photon_pt()->at(0)/b.llphoton_m()->at(0);
    double relp(-1);
    if(b.photon_pt()->at(1) > b.photon_pt()->at(0)) 
      if(b.photon_drmin()->at(1) > 0.3) 
        for(int i = 0; i < b.nllphoton(); i++)
          if(b.llphoton_iph()->at(i) == 1) {
            relp = b.photon_pt()->at(1)/b.llphoton_m()->at(i);
            break;
          }
    if(relp < 0) relp = b.photon_pt()->at(0)/b.llphoton_m()->at(0);
    return relp;
  });
  NamedFunc llp_m("llp_m",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.nphoton() == 1) return b.llphoton_m()->at(0);
    double m(-1);
    if(b.photon_pt()->at(1) > b.photon_pt()->at(0)) 
      if(b.photon_drmin()->at(1) > 0.3) 
        for(size_t i = 0; i < b.llphoton_pt()->size(); i++)
          if(b.llphoton_iph()->at(i) == 1) {
            m = b.llphoton_m()->at(i);
            break;
          }
    if(m < 0) m = b.llphoton_m()->at(0);
    return m;
  });
  NamedFunc baseline("nphoton > 0 && photon_id[0] && photon_elveto[0]");
  vector<NamedFunc> bline = {relpT > 0.12,
                             "ll_m[0] > 81.2 && ll_m[0] < 101.2",
                             "photon_drmin[0] > 0.3",
                             "llphoton_m[0] > 105 && llphoton_m[0] < 160"};
  for(size_t i = 0; i < bline.size(); i++)
    baseline = baseline && bline.at(i);
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 26 && el_pt[ll_i2[0]] > 17",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 26 && mu_pt[ll_i2[0]] > 8"};
  vector<float> windows = {3.7, 3.7, 3.8, 4.1, 3.9, 4.0, 4.0};
  vector<NamedFunc> mass;
  for(size_t im = 0; im < windows.size(); im++)
    mass.push_back(llp_m > 125.09 - windows.at(im)/2 && llp_m < 125.09 + windows.at(im)/2);
  lep.at(0) = lep.at(0) && el_check;
  NamedFunc vbf("vbf_mva > 0.925 && dijet_deta > 2");
//   NamedFunc relpT("photon_pt[0]/llphoton_m[0] > 0.4");
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.weight();
    if(b.type() == 6200)
      return weight/1.664;
    return        weight;
  });
  PlotMaker pm;
  pm.Push<Table>("ATLAS", vector<TableRow>{
      TableRow("VBF enriched",             baseline && mass.at(0) && vbf && (lep.at(0) || lep.at(1))                ,0,0,wgt),
      TableRow("High relative $p_T$",      baseline && mass.at(1) && relpT > 0.4 && !vbf && (lep.at(0) || lep.at(1)),0,0,wgt),
      TableRow("High $p_{Tt}$ ee",         baseline && mass.at(2) && pTt > 40 && lep.at(0) && !vbf && relpT < 0.40  ,0,0,wgt),
      TableRow("Low $p_{Tt}$ ee",          baseline && mass.at(3) && pTt < 40 && lep.at(0) && !vbf && relpT < 0.40  ,0,0,wgt),
      TableRow("High $p_{Tt}$ $\\mu\\mu$", baseline && mass.at(4) && pTt > 40 && lep.at(1) && !vbf && relpT < 0.40  ,0,0,wgt),
      TableRow("Low $p_{Tt}$ $\\mu\\mu$",  baseline && mass.at(5) && pTt < 40 && lep.at(1) && !vbf && relpT < 0.40  ,0,0,wgt),
      TableRow("Inclusive",                baseline && mass.at(6) && (lep.at(0) || lep.at(1))                       ,1,0,wgt)
	    },procs,false);
  pm.min_print_ = true;
  pm.MakePlots(139);
}


