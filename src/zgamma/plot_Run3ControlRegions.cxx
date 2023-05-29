#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Reader.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/slide_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

//------------------------------------GIST--------------------------------------//
//This plot file can be used for htozgamma_joshuatree_v0 or later versions of picos
//args: 
//2. basic/comprehenive/BDTs - provide a 0, 1, or 2 where 0 will plot only the sideband comparisons for mc and data
//   a 1 will plot various control regions. Both will have plots split by flavor and combined flavor.
//   An input of 2 will plot the various BDT regions. In jt_v1 picos there should be an int which gives
//   which untagged or dijet category the event belongs to.
//
//Last update 05/27/23 - Made to work only for basic or comprehenisive plot option
//------------------------------------------------------------------------------//


const NamedFunc signal_lead_electron_pt("signal_lead_electron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  for (unsigned iPart = 0; iPart<b.el_pt()->size(); iPart++) {
    if (b.el_sig()->at(iPart)) return b.el_pt()->at(iPart);
  }
  return -999;
});


const NamedFunc signal_sublead_electron_pt("signal_sublead_electron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  bool sublead = false;
  for (unsigned iPart = 0; iPart<b.el_pt()->size(); iPart++) {
    if (b.el_sig()->at(iPart)) {
      if (sublead == false) sublead = true;
      else return b.el_pt()->at(iPart);
      }
    }
  return -999;
});


const NamedFunc signal_lead_muon_pt("signal_lead_muon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  for (unsigned iPart = 0; iPart<b.mu_pt()->size(); iPart++) {
    if (b.mu_sig()->at(iPart)) return b.mu_pt()->at(iPart);
  }
  return -999;
});

const NamedFunc signal_sublead_muon_pt("signal_sublead_muon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  bool sublead = false;
    for (unsigned iPart = 0; iPart<b.mu_pt()->size(); iPart++) {
      if (b.mu_sig()->at(iPart)) {
        if (sublead == false) sublead = true;
      else return b.mu_pt()->at(iPart);
    }
  }
  return -999;
});


 NamedFunc photon_drmax("photon_drmax",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TVector3 photon = AssignGamma(b).Vect();
    TVector3 l1     = AssignL1(b).Vect();
    TVector3 l2     = AssignL2(b).Vect();
    return max(photon.DeltaR(l1),photon.DeltaR(l2));
  });


bool checkBit(int i, int n) {
  return((i%static_cast<int>(pow(2,n+1)))/static_cast<int>(pow(2,n)));
}

int main(int argc, char *argv[]){

  //These are the two input arguements, the first controls which input files to use the other is for what control regions to plot
  if(argc==1){cout<<"Please provide the option for basic/comprehensive control region plots and whether you want unskimmed (0) or merged_zg (1) data files"<<endl;return -1;}
  if(argc==2){cout<<"Please provide the option for basic/comprehensive control region plots"<<endl;return -1;}
  int fileSel = atoi(argv[1]);
  int plotSel = atoi(argv[2]);
 
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  NamedFunc type("rtype",  [](const Baby &b) -> NamedFunc::ScalarType{ return b.type(); });

  //This block selects whether or not to use unskimmed or merged_zgmc(data)_llg files
  string bfolder_dataEE;  string bfolder_mcEE;  string bfolder_data;  string bfolder_mc;
  if(fileSel==0){
    bfolder_dataEE = "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_joshuatree_v0/2022EE/data/unskimmed/";
    bfolder_mcEE = "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_joshuatree_v0/2022EE/mc/unskimmed/";
    bfolder_data = "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_joshuatree_v0/2022/data/unskimmed/";
    bfolder_mc = "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_joshuatree_v0/2022/mc/unskimmed/";
  } else if (fileSel==1){
    bfolder_dataEE = "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_joshuatree_v0/2022EE/data/merged_zgdata_llg/";
    bfolder_mcEE = "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_joshuatree_v0/2022EE/mc/merged_zgmc_llg/";
    bfolder_data = "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_joshuatree_v0/2022/data/merged_zgdata_llg/";
    bfolder_mc = "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_joshuatree_v0/2022/mc/merged_zgmc_llg/";
  }

  //This code defines the triggers and trigger pt selections. The muon triggers may be incorrect, but I tried to add all the triggers it could be
  NamedFunc trigs_2022_ee   = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"; 
  NamedFunc trigs_2022_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ||HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL||HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL";
  NamedFunc trigs_2022 =  (trigs_2022_ee && signal_lead_electron_pt > 25 && signal_sublead_electron_pt > 15) || (trigs_2022_mumu && signal_lead_muon_pt > 20 && signal_sublead_muon_pt > 10);  


  //These blocks define the processes that are used to make the plots.
  auto proc_data_2022EE  = Process::MakeShared<Baby_pico>("Data - 2022EE",     Process::Type::data, kBlack, {bfolder_dataEE +"*.root"}, trigs_2022);
  //{bfolder_dataEE + "*EGamma*.root",bfolder_dataEE + "*Muon*.root"}, trigs_2022);
  auto proc_DY_2022EE    = Process::MakeShared<Baby_pico>("DY + Fake - 2022EE", Process::Type::background, TColor::GetColor("#f7ad19"),{bfolder_mcEE + "*DYto2L*"},  trigs_2022);

  auto proc_data_2022  = Process::MakeShared<Baby_pico>("Data",     Process::Type::data, kBlack, {bfolder_data + "*.root"},trigs_2022);
  //{bfolder_data + "*EGamma*.root",bfolder_data + "*Muon*.root"}, trigs_2022);
  auto proc_DY_2022    = Process::MakeShared<Baby_pico>("DY + Fake", Process::Type::background, TColor::GetColor("#f7ad19"),{bfolder_mc + "*DYto2L*"},  trigs_2022);

  vector<string> procs_name = {"_2022", "_2022EE"};
  vector<vector<shared_ptr<Process>>> procs_w_data = {{proc_data_2022, proc_DY_2022},{proc_data_2022EE, proc_DY_2022EE}};



  //This code creates the plotting styles for the 1D and 2D histograms
  PlotOpt style("txt/plot_styles.txt","Eff2D");
  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::both)
          .YTitleOffset(1.)
          //.LinMaximum(20000)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(1077)
          .CanvasWidth(900)
          .FileExtensions({"pdf"});
          //.Bottom(BottomType::ratio); //Added this to find best cuts. remove this later 


  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .YTitleOffset(1.)
                                     .LogMinimum(0.1)
                                     //.LogMaximum(2000)
                                     .Overflow(OverflowType::both)
                                     .CanvasWidth(600)
                                     .LabelSize(0.04)
                                     .UseCMYK(true)  
                                     .YAxis(YAxisType::log)
                                     .Title(TitleType::info)};
 





  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt lin_dnorm = lin_lumi().Stack(StackType::data_norm);
  vector<PlotOpt> ops_dnorm_wrat = { lin_dnorm.Bottom(BottomType::ratio) };
  vector<PlotOpt> ops_2D = {lin_stack};


  //This bit of code defines the weights for each event. Data has a "-" in the sample type so it should not be given a weight. 
  //Right now the MC is not being weighted by any luminosity value
  NamedFunc wgt("wgt",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.SampleTypeString().Contains("-")) {
      return 1;
    }
/*    
    double w_year = 1;
    if (b.SampleTypeString()=="2022") {
        w_year = 13.5;
    } else if (b.SampleTypeString()=="2022EE") {
        w_year = 20.7;
    }
    return w_year*b.w_lumi()/34.2;
*/
    return b.weight();
 });

  //This defines the plotting object used to create plots
  PlotMaker pm;


  //These strings and namedfuncs are used to create the control regions and different plots based on lepton flavor
  vector<string> tag_cuts_flav = {"_EL","_MU","_LL"};
  vector<NamedFunc> flavorCut = {"ll_lepid[0]==11" , "ll_lepid[0]==13", "(ll_lepid[0]==11 || ll_lepid[0]==13)"};
  NamedFunc VBF_cuts = "njet>1 && dijet_m > 400";
  vector<string> CR_labels;
  vector<NamedFunc> ControlRegions;

  //These if statements use the second input argument to decide which control regions to plot
  if(plotSel==0){
    ControlRegions = {"nllphoton > 0 && photon_id80[llphoton_iph[0]] && ll_m[llphoton_ill[0]] > 80 && ll_m[llphoton_ill[0]] < 100 && (photon_pt[llphoton_iph[0]]/llphoton_m[0] > 15.0/110.0) && photon_drmin[llphoton_iph[0]] > 0.4 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && !(llphoton_m[0] < 128 && llphoton_m[0] > 122)"};
    CR_labels = {"_sideband"};

  }else if(plotSel==1){
    CR_labels  = {"_sideband","_mllyl100","_inv_mll","_inv_photonWP80","_inv_rat"};
    ControlRegions = {"nllphoton > 0 &&  photon_id80[llphoton_iph[0]] && ll_m[llphoton_ill[0]] > 80 && ll_m[llphoton_ill[0]] < 100 && (photon_pt[llphoton_iph[0]]/llphoton_m[0] > 15.0/110.0) && photon_drmin[llphoton_iph[0]] > 0.4 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && !(llphoton_m[0] < 128 && llphoton_m[0] > 122)", //Sideband
     "nllphoton > 0 &&  photon_id80[llphoton_iph[0]] && ll_m[llphoton_ill[0]] > 80 && ll_m[llphoton_ill[0]] < 100 && (photon_pt[llphoton_iph[0]]/llphoton_m[0] > 15.0/110.0) && photon_drmin[llphoton_iph[0]] > 0.4 && llphoton_m[0] < 100", //mlly < 100
     "nllphoton > 0 &&  photon_id80[llphoton_iph[0]] && !(ll_m[llphoton_ill[0]] > 80 && ll_m[llphoton_ill[0]] < 100) && (photon_pt[llphoton_iph[0]]/llphoton_m[0] > 15.0/110.0) && photon_drmin[llphoton_iph[0]] > 0.4 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && !(llphoton_m[0] < 128 && llphoton_m[0] > 122)", //!(80 < mll < 100)
"nllphoton > 0 &&  !photon_id80[llphoton_iph[0]] && ll_m[llphoton_ill[0]] > 80 && ll_m[llphoton_ill[0]] < 100 && (photon_pt[llphoton_iph[0]]/llphoton_m[0] > 15.0/110.0) && photon_drmin[llphoton_iph[0]] > 0.4 && llphoton_m[0] > 100 && llphoton_m[0] < 180", //!(photon_id80)
"nllphoton > 0 &&  photon_id80[llphoton_iph[0]] && ll_m[llphoton_ill[0]] > 80 && ll_m[llphoton_ill[0]] < 100 && (photon_pt[llphoton_iph[0]]/llphoton_m[0] < 15.0/110.0) && photon_drmin[llphoton_iph[0]] > 0.4 && llphoton_m[0] > 100 && llphoton_m[0] < 180" //inverted pt/mlly cut
     };
  }


//This is the start of the plotting. the plots are seperated by "year" (right now just 2022 and 2022EE), flavor, and control region
for(unsigned int k(0); k < flavorCut.size(); k++){
  for(unsigned int j(0); j < ControlRegions.size(); j++){
    for(unsigned int i(0); i < procs_w_data.size();  i++){
     pm.Push<Hist2D>(
       Axis(70,50 ,120,  "ll_m[0]",      "m_{ll} [GeV]",{50,80,100}),
       Axis(40,100,180,  "llphoton_m[0]","m_{ll#gamma} [GeV]",{100,180}),
       ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_2D).Tag("ShortName:plot_CR_2D_mll_mlly" + CR_labels[j] + procs_name[i]).Weight(wgt); 

      pm.Push<Hist1D>(Axis(40,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]",  {}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_Z_mlly" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
      pm.Push<Hist1D>(Axis(80,60,120, "ll_m[0]", "m_{ll} [GeV]",  {80,100}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_ZZ_mll" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);


    pm.Push<Hist1D>(Axis(8,0,8,  "nphoton", "N_{#gamma}", {}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_ROL1_nphoton_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(8,0,8,  "nel", "N_{e}", {}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_ROL1_nel_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(8,0,8,  "nmu", "N_{#mu}", {}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_ROL1_nmu_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);

    pm.Push<Hist1D>(Axis(8,0,8,  "njet", "N_{jets}", {}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_ROL1_njet_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);


    //-- Kinematic IDMVA Variables --//
    pm.Push<Hist1D>(Axis(40,-1,1,  "llphoton_costheta[0]", "cos(#theta)", {}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_6_ctheta_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,-1,1,  "llphoton_cosTheta[0]", "cos(#Theta)", {}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_4_cTHETA_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,-3.2,3.2,  "llphoton_psi[0]", "#phi", {}), ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_Z11_phi_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(30,0.05,0.65, "llphoton_pt[0]/llphoton_m[0]", "p_{T}(ll#gamma)/m_{ll#gamma}", {15.0/110}), ControlRegions[j] && flavorCut[k] , procs_w_data[i] , ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_3_ptlly_mlly_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);

    //Creating the correct eta plots
    if(k==0 || k==2){
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,  "el_eta[ll_i1[0]]", "#eta(e_{1})", {}), ControlRegions[j] && flavorCut[0], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_8a_eta_el1_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,  "el_eta[ll_i2[0]]", "#eta(e_{2})", {}), ControlRegions[j] && flavorCut[0], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_9a_eta_el2_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,5,65,  "el_pt[ll_i1[0]]", "p_{T}(e_{1})", {}), ControlRegions[j] && flavorCut[0], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_88a_pt_el1_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(45,10,100,  "el_pt[ll_i2[0]]", "p_{T}(e_{2})", {}), ControlRegions[j] && flavorCut[0], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_99a_pt_el2_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    }
    if(k==1 || k==2) {
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,  "mu_eta[ll_i1[0]]", "#eta(#mu_{1})", {}), ControlRegions[j] && flavorCut[1], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_8b_eta_mu1_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,  "mu_eta[ll_i2[0]]", "#eta(#mu_{2})", {}), ControlRegions[j] && flavorCut[1], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_9b_eta_mu2_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,0,50,  "mu_pt[ll_i1[0]]", "p_{T}(#mu_{1})", {}), ControlRegions[j] && flavorCut[1], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_88b_pt_mu1_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,0,50,  "mu_pt[ll_i2[0]]", "p_{T}(#mu_{2})", {}), ControlRegions[j] && flavorCut[1], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_99b_pt_mu2_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    }

    //Photon variable plots 
    pm.Push<Hist1D>(Axis(25,10,60,  "photon_pt[llphoton_iph[0]]", "p_{T}(#gamma)", {}), ControlRegions[j] && flavorCut[k], procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_88b_pt_mu1_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,  "photon_eta[llphoton_iph[0]]", "#eta(#gamma)", {}), ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_Z10_eta_gamma_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(25,0,3.5, "photon_drmin[llphoton_iph[0]]", "min( #Delta R(#gamma,l) )",  {}), ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_2_CoarseminDR_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(20,0.01,0.25, "photon_pterr[llphoton_iph[0]]/photon_pt[llphoton_iph[0]]", "#sigma(#gamma)/E(#gamma)", {}), ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_5_sigma_gamma_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(40,0.3,1,  "photon_idmva[llphoton_iph[0]]", "#gamma - IDMVA", {}), ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_1_IDMVA_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(45,0.5,5, photon_drmax, "max( #Delta R(#gamma,l) )",  {}), ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_7_maxDR_" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);



    //VBF Plots
    pm.Push<Hist1D>(Axis(90,30,300 , "jet_pt[0]"    ,"Leading jet p_{T} [GeV]"   ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_9_jet1_pt" + tag_cuts_flav[k]  + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(90,30,300 , "jet_pt[1]"    ,"Subleading jet p_{T} [GeV]",{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_8_jet2_pt" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);


    pm.Push<Hist1D>(Axis(12,400,800 , "dijet_m "    ,"m_{jj} [GeV]"   ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_9_jet1_pt" + tag_cuts_flav[k]  + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(12,0,9    ,"dijet_deta","#Delta#eta(j_{1},j_{2})"   ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_1_deta" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(10,0,3    ,"dijet_dphi","#Delta#phi(j_{1},j_{2})"   ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_6_dphi" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);

    pm.Push<Hist1D>(Axis(10, 0,  1 , "llphoton_dijet_balance"    ,"System balance"            ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_4_sysbal" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(10,0,3    , "llphoton_dijet_dphi[0]" ,"#Delta#phi(Z#gamma,jj)"    ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_3_jjlly_dphi" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis(20,0,140  , "llphoton_pTt2[0]"        ,"p_{Tt} [GeV]"              ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_5_pTt2" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);

    pm.Push<Hist1D>(Axis(15,0.4,4.4, "photon_jet_mindr[llphoton_iph[0]]"  ,"Min. #DeltaR(#gamma,j)"         ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_10_yjdr" + tag_cuts_flav[k] + CR_labels[j] + procs_name[i]);
    pm.Push<Hist1D>(Axis( 8,0,6    ,"photon_zeppenfeld[llphoton_iph[0]]"  ,"#gamma zeppenfeld"         ,{}), VBF_cuts && ControlRegions[j] && flavorCut[k] , procs_w_data[i], ops_dnorm_wrat).Weight(wgt).Tag("ShortName:plot_CR_dijet_7_zep" + tag_cuts_flav[k]  + CR_labels[j] + procs_name[i]);
    }
  }
}


  //This chunk of code is needed to actually create the plots. The 1.0 can be changed to a luminosity value. 
  //Right now it is set to 2022 luminosity value but will need to be changed for the luminosity of 2022 and 2022EE
  pm.min_print_ = true;
  pm.MakePlots(1.01);
}


