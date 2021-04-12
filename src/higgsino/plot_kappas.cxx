// Cheatsheet:
//    To obtain plots with no mismeasurement and kappa's with expected data stats:
//         ./run/hig/plot_kappa.exe --mm mc_as_data -l 36.2
//    There are 4 possibilities for the sample, requested with option -s. These are: search, zll, qcd, ttbar
//    Option -t plots the kappas with a tighter selection, see basecuts below, e.g.
//         ./run/hig/plot_kappa.exe --mm mc_as_data -t -s zll -l 36.2

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>
#include <string>
#include <algorithm>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "RooStats/NumberCountingUtils.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"


#include "TVector2.h"
#include "TMath.h"
#include "Math/Vector4D.h"

using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
using std::cout;
using std::setw;
using std::endl;
using std::map;
using std::set;
using std::shared_ptr;
using std::pair;
using std::tuple;
using std::get;

namespace{
  bool split_bkg = true;
  bool do_signal = true;
  bool do_zbi = true;
  bool do_incl_met = true;
  bool do_bin_drmax = true;
  bool do_old_bin = false;
  bool debug = false;
  bool do_highnb = false;
  bool do_midnb = false;
  bool unblind = false;
  int digits_table = 2;
  TString sample = "search";
  string alt_scen = "mc_as_data"; //e.g. "mc", "data", "mc_as_data" or any systematic defined in sys_weights.cfg or in the defining scenarios section below
  float lumi=1.;
  bool quick_test = false;
  // for office use only
  vector<TString> syst_names;
  vector<vector<float>> syst_values;
  string sys_wgts_file = "txt/sys_weights.cfg";
  string year_string = "2016,2017,2018";
  // string_options is split by comma. ex) option1,option2 
  // Use HigUtilities::is_in_string_options(string_options, "option2") to check if in string_options.
  // Options: use_old_trigger,split_ttbar_met,search_ttbar_same_bins,save_entries_weights_to_file,print_entries_weights
  string string_options = "";
}

struct abcd_def{
  TString scenario;
  vector<TString> planecuts;
  vector<TString> bincuts;
};

TString printTable(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds, 
                   vector<vector<float> > yieldsPlane, vector<shared_ptr<Process> > &proc_sigs, float total_luminosity);
void plotKappa(abcd_def &abcd, vector<vector<vector<float> > >  &kappas, 
               vector<vector<vector<float> > >  &kappas_mm, vector<vector<vector<float> > >  &kmcdat, float total_luminosity);
void plotTable(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds, 
                   vector<vector<float> > yieldsPlane, vector<shared_ptr<Process> > &proc_sigs, float total_luminosity);
vector<vector<float> > findPreds(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
                                 vector<vector<vector<float> > > &kappas, 
                                 vector<vector<vector<float> > > &kappas_mm, 
                                 vector<vector<vector<float> > > &kmcdat, 
                                 vector<vector<vector<float> > > &preds, bool saveData=false);


void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  std::chrono::high_resolution_clock::time_point begTime;
  begTime = std::chrono::high_resolution_clock::now();

  Palette colors("txt/colors.txt", "default");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// What ntuples to use?  //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<int> years; 
  HigUtilities::parseYears(year_string, years);
  //years = {2016,2017,2018};
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  string higgsino_version = "";
  // trigger_version: 0: old, 1: new
  int trigger_version = !HigUtilities::is_in_string_options(string_options, "use_old_trigger");

  string base_dir(string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/");
  string mc_skim_dir("mc/merged_higmc_higloose/"), data_skim_dir("data/merged_higdata_higloose/"), sig_skim_dir("SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/");
  if (sample=="ttbar")    {mc_skim_dir = "mc/merged_higmc_higlep1T/"; data_skim_dir = "data/merged_higdata_higlep1T/"; sig_skim_dir = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/";} 
  else if (sample=="zll") {mc_skim_dir = "mc/merged_higmc_higlep2T/"; data_skim_dir = "data/merged_higdata_higlep2T/"; sig_skim_dir = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep2T/";} 
  else if (sample=="qcd") {mc_skim_dir = "mc/merged_higmc_higqcd/";  data_skim_dir = "data/merged_higdata_higqcd/"; sig_skim_dir = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higqcd/";} 

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({
                                    "*TTJets_*Lept*", 
                                    "*_TTZ*.root", "*_TTW*.root",
                                    "*_TTGJets*.root", "*_ttHTobb*.root","*_TTTT*.root"
                                  });
  mctags["vjets"]   = set<string>({ 
                                    "*DYJetsToLL*.root","*_ZJet*.root", "*_WJetsToLNu*.root"
                                 });
  mctags["other"]   = set<string>({ // no events pass from QCD HT 100-200
                                   // "*QCD_HT100to200_Tune*", 
                                   "*QCD_HT200to300_Tune*", 
                                   "*QCD_HT300to500_Tune*",  // these have very low stats
                                   "*QCD_HT500to700_Tune*",
                                   "*QCD_HT700to1000_Tune*", 
                                   "*QCD_HT1000to1500_Tune*", 
                                   "*QCD_HT1500to2000_Tune*", 
                                   "*QCD_HT2000toInf_Tune*",
                                   "*_ST_*.root",
                                   "*_WH*.root", "*_ZH_HToBB*.root",
                                   "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"
                                 });

  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining baseline cuts ///////////////////////////////////////
  TString c_hig_trim = "hig_cand_drmax[0]<2.2  && hig_cand_dm[0]<40 && hig_cand_am[0]<200";
  string baseline_s = "njet>=4 && njet<=5";
  if (sample=="search") baseline_s += " && nvlep==0 && ntk==0 && !low_dphi_met &&"+c_hig_trim;
  else if (sample=="ttbar") baseline_s += " && nlep==1 && mt<100 &&"+c_hig_trim;
  else if (sample=="zll") baseline_s += " && nlep==2 && met<50 &&"+c_hig_trim;
  else if (sample=="qcd") baseline_s += " && nvlep==0 && ntk==0 && low_dphi_met &&"+c_hig_trim;

  NamedFunc baseline = baseline_s; 
  if (sample=="search") baseline = baseline && Higfuncs::final_pass_filters;
  else if (sample=="ttbar") baseline = baseline && Higfuncs::final_ttbar_pass_filters && Higfuncs::lead_signal_lepton_pt>30;
  else if (sample=="zll") baseline = baseline && Higfuncs::final_zll_pass_filters;
  else if (sample=="qcd") baseline = baseline && Higfuncs::final_qcd_pass_filters;

  if (sample=="search") do_incl_met = false;
  else do_incl_met = true;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////  
  // define signal processes
  vector<string> sigMasses;
  if (do_old_bin) sigMasses = {"225", "400", "700"};
  else sigMasses = {"175", "500", "950"};
  vector<shared_ptr<Process> > proc_sigs;
  for (unsigned isig(0); isig<sigMasses.size(); isig++)
    proc_sigs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigMasses[isig]+",1)", Process::Type::signal, 1, 
      attach_folder(base_dir, years, sig_skim_dir,{"*TChiHH_mChi-"+sigMasses[isig]+"*mLSP-0_*.root"}), baseline));

  // define background processes
  auto proc_ttx = Process::MakeShared<Baby_pico>("tt+X", Process::Type::background, colors("tt_1l"),
    attach_folder(base_dir, years, mc_skim_dir, mctags["ttx"]), baseline && "stitch");
  auto proc_vjets = Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
    attach_folder(base_dir, years, mc_skim_dir, mctags["vjets"]), baseline && "stitch");
  auto proc_other = Process::MakeShared<Baby_pico>("Other", Process::Type::background, colors("other"),
    attach_folder(base_dir, years, mc_skim_dir, mctags["other"]), baseline && "stitch");

  // define data or pseudo-data process
  set<string> names_data(attach_folder(base_dir, years, data_skim_dir, {"*.root"}));
  if(alt_scen != "data") { 
    if (quick_test) {
      names_data = {base_dir+"/2016/"+mc_skim_dir+"/*TTJets_SingleLeptFromT_Tune*"};
    } else {
      set<string> names_allmc;
      for (auto &iset: mctags) 
        names_allmc.insert(iset.second.begin(), iset.second.end());    
      names_data = attach_folder(base_dir, years, mc_skim_dir, names_allmc);
    }
  }
  NamedFunc base_data = baseline; 
  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc met_triggers = Higfuncs::met_trigger;
  if (trigger_version == 0) {
    lepton_triggers = Higfuncs::el_trigger_v0 || Higfuncs::mu_trigger_v0;
    met_triggers = Higfuncs::met_trigger_v0;
  }
  if (sample=="zll") base_data = base_data && lepton_triggers;
  else if (sample == "ttbar") base_data = base_data && (lepton_triggers || met_triggers);
  else if (sample == "qcd") base_data = base_data && met_triggers;
  else  base_data = base_data && met_triggers; // search
  if (alt_scen != "data") base_data = baseline && "stitch";
  auto proc_data = Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, names_data,base_data);

  // define combination of processes to use
  vector<shared_ptr<Process> > all_procs = vector<shared_ptr<Process> >{proc_ttx, proc_vjets, proc_other};
  if(quick_test) {
    auto proc_bkg = Process::MakeShared<Baby_pico>("All_bkg", Process::Type::background, colors("tt_1l"),
                    {base_dir+"/2016/"+mc_skim_dir+"/*TTJets_SingleLeptFromT*"}, baseline && "stitch");
    all_procs = vector<shared_ptr<Process> >{proc_bkg};
    split_bkg = false;
  }
  if (do_signal){
    for(size_t ind=0; ind<proc_sigs.size(); ind++)
      all_procs.push_back(proc_sigs[ind]);
  }
  all_procs.push_back(proc_data);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining scenarios  //////////////////////////////////////////
  NamedFunc nom_wgt = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years*Functions::w_pileup;
  NamedFunc wgt_syst_trg_eff = "weight"*Higfuncs::w_years*Functions::w_pileup;
  if (trigger_version == 0) nom_wgt = "weight"*Higfuncs::eff_higtrig_run2_v0*Higfuncs::w_years*Functions::w_pileup;

  vector<string> scenarios;
  map<string, NamedFunc> weights;
  weights.emplace("nominal", nom_wgt);
  if(alt_scen == "data"){
    scenarios = vector<string>{"data"};
  } else if(alt_scen == "bctag"){ 
    // example custom-defined systematic, can define multiple instead of just one, if desired
    // in general, can define any distortion of the MC with a named func and compare to default using this...
    scenarios = vector<string>();
    scenarios.push_back("syst_bctag");
    weights.emplace("syst_bctag", nom_wgt*"sys_bctag[0]");
  } else if(alt_scen == "syst_ttx_up"){
    scenarios = vector<string>();
    scenarios.push_back("syst_ttx_up");
    weights.emplace("syst_ttx_up", nom_wgt/(1+Higfuncs::wgt_syst_ttx));
  } else if(alt_scen == "syst_ttx_down"){
    scenarios = vector<string>();
    scenarios.push_back("syst_ttx_down");
    weights.emplace("syst_ttx_down", nom_wgt/(1-Higfuncs::wgt_syst_ttx));
  } else if(alt_scen == "syst_vjets_up"){
    scenarios = vector<string>();
    scenarios.push_back("syst_vjets_up");
    weights.emplace("syst_vjets_up", nom_wgt/(1+Higfuncs::wgt_syst_vjets));
  } else if(alt_scen == "syst_vjets_down"){
    scenarios = vector<string>();
    scenarios.push_back("syst_vjets_down");
    weights.emplace("syst_vjets_down", nom_wgt/(1-Higfuncs::wgt_syst_vjets));
  } else if(alt_scen == "syst_qcd_down"){
    scenarios = vector<string>();
    scenarios.push_back("syst_qcd_down");
    weights.emplace("syst_qcd_down", nom_wgt/(1-Higfuncs::wgt_syst_qcd));
  } else if(alt_scen == "syst_qcd_up"){
    scenarios = vector<string>();
    scenarios.push_back("syst_qcd_up");
    weights.emplace("syst_qcd_up", nom_wgt/(1+Higfuncs::wgt_syst_qcd));
  } else if(alt_scen == "syst_comp"){
    scenarios = vector<string>();
    scenarios.push_back("syst_comp");
    weights.emplace("syst_comp", nom_wgt*Higfuncs::wgt_syst_comp);
  } else if(alt_scen == "syst_trg_eff"){
    scenarios = vector<string>();
    scenarios.push_back("syst_trg_eff");
    weights.emplace("syst_trg_eff", wgt_syst_trg_eff);
  } else if(alt_scen == "syst_lowptleptrig"){
    scenarios = vector<string>();
    scenarios.push_back("syst_lowptleptrig");
    weights.emplace("syst_lowptleptrig", nom_wgt*Higfuncs::wgt_syst_lowptleptrig);
  } else if(alt_scen == "syst_ttpt"){
    scenarios = vector<string>();
    scenarios.push_back("syst_ttpt");
    weights.emplace("syst_ttpt", nom_wgt*Higfuncs::wgt_syst_ttbarpt);
  } else if(alt_scen == "mc_as_data"){
    scenarios = vector<string>{"mc_as_data"}; 
    weights.emplace("mc_as_data", nom_wgt);
  } else {
    // run on all scenarios from the sys_cfg file
    scenarios = vector<string>{alt_scen}; 
    for(const auto &scen: scenarios)
      weights.emplace(scen, nom_wgt*Functions::MismeasurementWeight(sys_wgts_file, scen));
  } 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining binning    //////////////////////////////////////////
  bool split_ttbar_met = HigUtilities::is_in_string_options(string_options, "split_ttbar_met");
  vector<TString> metcuts;
  string metdef = "met";
  if (sample=="zll") metdef = "ll_pt[0]";
  if (sample=="search" || sample=="qcd"){
    if (do_old_bin) {
      metcuts.push_back(metdef+">150&&"+metdef+"<=200");
      metcuts.push_back(metdef+">200&&"+metdef+"<=300");
      metcuts.push_back(metdef+">300&&"+metdef+"<=450");
      metcuts.push_back(metdef+">450");
    } else if (HigUtilities::is_in_string_options(string_options, "search_ttbar_same_bins")) {
      metcuts.push_back(metdef+">150&&"+metdef+"<=200");
      metcuts.push_back(metdef+">200&&"+metdef+"<=300");
      metcuts.push_back(metdef+">300");
    } else {
      metcuts.push_back(metdef+">150&&"+metdef+"<=200");
      metcuts.push_back(metdef+">200&&"+metdef+"<=300");
      metcuts.push_back(metdef+">300&&"+metdef+"<=400");
      metcuts.push_back(metdef+">400");
    }
  } else if (sample=="zll") {
    metcuts.push_back(metdef+">0&&"+metdef+"<=75");
    metcuts.push_back(metdef+">75&&"+metdef+"<=150");
    metcuts.push_back(metdef+">150&&"+metdef+"<=200");
    metcuts.push_back(metdef+">200&&"+metdef+"<=300");
    metcuts.push_back(metdef+">300");
  } else if (sample=="ttbar") {
    if (split_ttbar_met) {
      metcuts.push_back(metdef+">0&&"+metdef+"<=75");
      metcuts.push_back(metdef+">75&&"+metdef+"<=150");
      metcuts.push_back(metdef+">150&&"+metdef+"<=200");
      metcuts.push_back(metdef+">200&&"+metdef+"<=300");
      metcuts.push_back(metdef+">300");
    } else if (HigUtilities::is_in_string_options(string_options, "search_ttbar_same_bins")) {
      metcuts.push_back(metdef+">150&&"+metdef+"<=200");
      metcuts.push_back(metdef+">200&&"+metdef+"<=300");
      metcuts.push_back(metdef+">300");
    }else {
      metcuts.push_back(metdef+">0&&"+metdef+"<=75");
      metcuts.push_back(metdef+">75&&"+metdef+"<=150");
      metcuts.push_back(metdef+">150");
    }
  }
  if (do_incl_met) {
    if (sample=="qcd" || sample=="search") { // add an inclusive bin
      metcuts.push_back(metdef+">150");
    } else if (sample=="ttbar" ) {
      if (!split_ttbar_met || !HigUtilities::is_in_string_options(string_options, "search_ttbar_same_bins")) {
        metcuts.push_back(metdef+">0");
      }
    } else if (sample=="zll"){
      metcuts.push_back(metdef+">0");
    }
  }

  vector<TString> nbcuts;
  if (sample=="ttbar" || sample=="search" || do_highnb){
    nbcuts.push_back("nbt==2&&nbm==2");
    nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
    nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
  } else if (do_midnb){
    nbcuts.push_back("nbm==1");
    nbcuts.push_back("nbm==2");
  } else {
    nbcuts.push_back("nbm==0");
    nbcuts.push_back("nbm==1");
  }

  vector<TString> planecuts;
  for (unsigned imet(0); imet<metcuts.size(); imet++){
    if (!do_bin_drmax || do_old_bin) {
      planecuts.push_back(metcuts[imet]);
    } else {
      planecuts.push_back(metcuts[imet] + "&& hig_cand_drmax[0]<=1.1");
      planecuts.push_back(metcuts[imet] + "&& hig_cand_drmax[0]>1.1");
    }
  }

  ////// ABCD cuts
  vector<TString> abcd_regions = {"!(hig_cand_am[0]>100 && hig_cand_am[0]<=140) && nb_cr",
                                  "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140) && nb_sr",
                                  "hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nb_cr",
                                  "hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nb_sr"};

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining ABCD scenarios //////////////////////////////////////////
  PlotMaker pm;

  vector<abcd_def> abcds;
  for(auto &iscen: scenarios) {
    // the bins are the bins in the signal region, so skip the first Nb bin, which is the control
    vector<TString> bincuts = vector<TString>(nbcuts.begin()+1, nbcuts.end());

    abcd_def iabcd{iscen, planecuts, bincuts};
    abcds.push_back(iabcd);

    vector<TableRow> table_cuts, table_cuts_mm;
    for(auto &iplane: planecuts) { // N.B. loop order determines order of yields, so do not change...
      for(auto &inb: bincuts){
        for(auto &ireg: abcd_regions){
          TString totcut = iplane+"&&"+ireg;
          totcut.ReplaceAll("nb_cr", nbcuts[0]).ReplaceAll("nb_sr",inb);
          table_cuts.push_back(TableRow(totcut.Data(), totcut.Data(),0,0,weights.at("nominal")));
          if(alt_scen != "data") 
            table_cuts_mm.push_back(TableRow(totcut.Data(), totcut.Data(),0,0,weights.at(iscen)));
        } // Loop over ABCD cuts
      } // Loop over bin cuts
    } // Loop over plane cuts
    TString tname = "preds"; tname += iscen;
    // all_procs is a list of processes
    pm.Push<Table>(tname.Data(),  table_cuts, all_procs, true, false);
    tname += iscen;
    //pm.Push<Table>(tname.Data(),  table_cuts_mm, all_procs, true, false);
    pm.Push<Table>(tname.Data(),  table_cuts_mm, all_procs, true, false);
  } // Loop over ABCD scenarios

  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  ////////////////////////// Save the yields //////////////////////////////////////
  vector<TString> tablenames;
  for(size_t iscen=0; iscen<abcds.size(); iscen++) {
    // allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
    // if split_bkg [2/4] Other, [3/5] ttx, [4/6] vjets
    // allyields[dataType][bin] = GammaParams
    vector<vector<GammaParams> > allyields;
    Table * yield_table;
    if(alt_scen != "data"){
      yield_table = static_cast<Table*>(pm.Figures()[iscen*2].get());
      Table * yield_table_mm = static_cast<Table*>(pm.Figures()[iscen*2+1].get());
      // Yield is used to get info from specific process. 
      // DataYield can't be used because weight is assumed to be 1.
      allyields.push_back(yield_table_mm->Yield(proc_data.get(), lumi));
    } else {
      yield_table = static_cast<Table*>(pm.Figures()[iscen].get());
      allyields.push_back(yield_table->DataYield());
    }
    allyields.push_back(yield_table->BackgroundYield(lumi));
    if(do_signal){
      for(size_t ind=0; ind<proc_sigs.size(); ind++)
        allyields.push_back(yield_table->Yield(proc_sigs[ind].get(), lumi));
    }
    if(split_bkg){
      allyields.push_back(yield_table->Yield(proc_other.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_vjets.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_ttx.get(), lumi));
    }
    // Print out all the yields
    if (debug) {
      cout<<"Total number of cuts per process:"<<allyields[0].size()<<endl;
      cout<<"----------------------- Yields table -----------------------"<<endl;
      cout<<setw(10)<<"Cut"<<setw(10)<<"Data"<<setw(10)<<"SM bkg.";
      cout<<setw(10)<<"Signal 1"<<setw(10)<<"Signal 2";
      if (split_bkg) cout<<setw(10)<<"tt+X"<<setw(10)<<"Z+jets"<<setw(10)<<"Other";
      cout<<endl;
      for (unsigned j=0; j<allyields[0].size(); j++){
        if (j%(abcd_regions.size()*2)==0)
          cout<<"Plane: "<<planecuts[j/(abcd_regions.size()*2)]<<endl; // planecuts+1 because of the inclusive MET bin
        cout<<setw(10)<<CodeToRootTex(abcd_regions[j%abcd_regions.size()].Data());
        for (unsigned i=0; i<allyields.size(); i++){
          cout<<setw(10)<<RoundNumber(allyields[i][j].Yield(),1);
        }
        cout<<endl;
      }
      cout<<"-----------------------    End.    -----------------------"<<endl;
    }

    //// Calculating kappa and Total bkg prediction
    vector<vector<vector<float> > > kappas, kappas_mm, kmcdat, preds;
    if (debug) cout<<"Finding predictions"<<endl;
    bool save_entries_weights_to_file = HigUtilities::is_in_string_options(string_options, "save_entries_weights_to_file");
    vector<vector<float> > yieldsPlane = findPreds(abcds[iscen], allyields, kappas, kappas_mm, kmcdat, preds, /*saveData*/ save_entries_weights_to_file);

    if (debug) cout<<"Making tables."<<endl;
    TString fullname = printTable(abcds[iscen], allyields, kappas, preds, yieldsPlane, proc_sigs, total_luminosity);
    tablenames.push_back(fullname);
    
    //// Plotting kappa
    if (debug) cout<<"Making kappa plots."<<endl;
    plotKappa(abcds[iscen], kappas, kappas_mm, kmcdat, total_luminosity);

    if (sample=="search" && alt_scen=="data") plotTable(abcds[iscen], allyields, kappas, preds, yieldsPlane, proc_sigs, total_luminosity);

  } // Loop over ABCD methods

  if(alt_scen=="data" || alt_scen=="mc_as_data" || alt_scen=="mc"){
    //// Printing names of ouput files
    cout<<endl<<"===== Tables that can be compiled"<<endl;
    for(size_t ind=0; ind<tablenames.size(); ind++)
      cout<<" pdflatex "<<tablenames[ind]<<" &> /dev/null"<<endl;
    cout<<endl;
  }

  double seconds = (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding "<<abcds.size()<<" tables took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// Prints table with results
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
// if split_bkg: [2/4] Other, [3/5] tt1l, [4/6] tt2l
TString printTable(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds, 
                   vector<vector<float> > yieldsPlane, vector<shared_ptr<Process> > &proc_sigs, float total_luminosity){
  //cout<<endl<<"Printing table (significance estimation can take a bit)"<<endl;
  //// Table general parameters
  yieldsPlane.clear();
  TString ump = " & ";

  //// Setting output file name
  TString lumi_s = RoundNumber(total_luminosity, 0);
  if (lumi_s=="1") lumi_s = "137";
  TString outname = "tables/table_pred_"+sample+"_lumi"+lumi_s+(do_highnb?"_highnb":"")+(do_midnb?"_midnb":""); 
  outname += "_"+abcd.scenario+".tex";
  ofstream out(outname);

  //// Printing main table preamble
  out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\\begin{tabular}[tbp!]{ l ";
  size_t Nsig = proc_sigs.size(); // Number of signal points (for now it cannot be changed)
  size_t Ncol = 1;
  if(split_bkg) {out << "|ccc"; Ncol +=3;}
  out << "|cc"; Ncol +=2; //for kappa and tot bkg
  out << "|c"; Ncol += 1;
  if (unblind || alt_scen == "data") {
    out << "c "<<(do_zbi?"c":""); // for predicted and observed
    Ncol += 1 + (do_zbi?1:0);
  }
  if(do_signal) {
    for(size_t ind=0; ind<Nsig; ind++){
      out<<"|c"<<(do_zbi?"c":"");
      Ncol += 1 + (do_zbi?1:0);
    } 
  }

  out<<"}\\hline\\hline\n";
  out<<"${\\cal L}="<<lumi_s<<"$ fb$^{-1}$ ";
  if(split_bkg) out << " & Other & V$+$jets & $t\\bar{t}$ ";
  out << "& $\\kappa$ & MC bkg.";
  out << " & Pred. ";
  if (unblind || alt_scen == "data") out << "& Obs. "<<(do_zbi?"& Signi.":"");
  if(do_signal) {
    for(size_t ind=0; ind<Nsig; ind++) {
      TString signame = proc_sigs[ind]->name_.c_str();
      if(do_zbi) out << "& \\multicolumn{2}{c"<<(ind<Nsig-1?"|":"")<<"}{" << signame <<"}";
      else  out << "& " << signame;
    }
  }
  out << " \\\\ \\hline\\hline\n";

  vector<TString> binNames({"SBD, crb", "SBD, xb", "HIG, crb", "HIG, xb"});
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// Printing results////////////////////////////////////////////////
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    out<<endl<< "\\multicolumn{"<<Ncol<<"}{c}{$"<<CodeToLatex(abcd.planecuts[iplane].Data())
       <<"$ }  \\\\ \\hline\n";
    for(size_t ibin=0; ibin < abcd.bincuts.size(); ibin++){
      for(auto iabcd:{0,2,1,3}){
        if (ibin>0 && iabcd%2==0) continue; // don't print the 2b regions again
        size_t index = iplane*abcd.bincuts.size()*4+ibin*4+iabcd;
        if(iabcd==1) out << "\\hline" << endl;
        //// Printing bin name
        TString binName = binNames[iabcd];
        if (!do_highnb && (sample=="zll" || sample=="qcd")) {
          if(do_midnb) binName.ReplaceAll("xb", "2b").ReplaceAll("crb", "1b");
          else binName.ReplaceAll("xb", "1b").ReplaceAll("crb", "0b");
        } else {
          if(ibin==0) binName.ReplaceAll("xb", "3b");
          else binName.ReplaceAll("xb", "4b");
          binName.ReplaceAll("crb", "2b");
        }
        out << binName;
        //// Printing Other, tt1l, tt2l
        if(split_bkg){
          size_t offset = (do_signal?Nsig:0);
          out << ump <<RoundNumber(allyields[offset+2][index].Yield(), digits_table) << "$\\pm$" << RoundNumber(allyields[offset+2][index].Uncertainty(), digits_table)
              << ump <<RoundNumber(allyields[offset+3][index].Yield(), digits_table) << "$\\pm$" << RoundNumber(allyields[offset+3][index].Uncertainty(), digits_table)
              << ump <<RoundNumber(allyields[offset+4][index].Yield(), digits_table) << "$\\pm$" << RoundNumber(allyields[offset+4][index].Uncertainty(), digits_table);
        }
        //// Printing kappa
        out<<ump;
        if(iabcd==3) out  << "$"    << RoundNumber(kappas[iplane][ibin][0], digits_table)
                          << "^{+"  << RoundNumber(kappas[iplane][ibin][1], digits_table)
                          << "}_{-" << RoundNumber(kappas[iplane][ibin][2], digits_table) <<"}$ ";
        //// Printing MC Bkg yields
        out << ump << RoundNumber(allyields[1][index].Yield(), digits_table) << "$\\pm$" << RoundNumber(allyields[1][index].Uncertainty(), digits_table);
        //// Printing background predictions
        out << ump;
        if(iabcd==3) out << "$"    << RoundNumber(preds[iplane][ibin][0], digits_table)
                         << "^{+"  << RoundNumber(preds[iplane][ibin][1], digits_table)
                         << "}_{-" << RoundNumber(preds[iplane][ibin][2], digits_table) <<"}$ ";
        if (unblind || alt_scen == "data") {
          //// Printing observed events in data and Obs/MC ratio
          out << ump;
          if(iabcd==3) { // signal region
            if (unblind) out << RoundNumber(allyields[0][index].Yield(), 0);
          }
          else out << RoundNumber(allyields[0][index].Yield(), 0);
          //// Printing Zbi significance
          if(do_zbi) { // "$"+RoundNumber(Significance( ),1)+"\\sigma$"
            out << ump;
            //if(iabcd==3) out << "$"+RoundNumber(Significance(allyields[0][index].Yield(), preds[iplane][ibin][0], 
            //                                      preds[iplane][ibin][1], preds[iplane][ibin][2]),1)+"\\sigma$";
            if(iabcd==3) out << "$"+RoundNumber(utilities::to_significance(utilities::calculate_pvalue(allyields[0][index].Yield(), preds[iplane][ibin][0], 
                                                  preds[iplane][ibin][1], preds[iplane][ibin][2])),1)+"\\sigma$";
          }
        }
        //// Printing signal yields
        if(do_signal){
          for(size_t ind=0; ind<Nsig; ind++) {
            out<<ump<<RoundNumber(allyields[2+ind][index].Yield(), digits_table);
            if(do_zbi){
              out << ump;
              if(iabcd==3) {
                //float signif = Significance(preds[iplane][ibin][0]+allyields[2+ind][index].Yield(),preds[iplane][ibin][0], preds[iplane][ibin][1], preds[iplane][ibin][2]);
                float signif = utilities::to_significance(utilities::calculate_pvalue(preds[iplane][ibin][0]+allyields[2+ind][index].Yield(),preds[iplane][ibin][0], preds[iplane][ibin][1], preds[iplane][ibin][2]));
                if (allyields[2+ind][index].Yield()>0.45 && signif>0.45) {
                  out<<"$"+RoundNumber(signif,1)+"\\sigma$";
                } else {
                  out<<"$ - $";
                }
              }
            } // if do_zbi
          } // Loop over signals
        } // if do_signal
        out << "\\\\ \n";
      } // Loop over bin cuts
    } // Loop over ABCD cuts
    out << "\\hline\\hline\n";
  } // Loop over plane cuts
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //// Printing footer and closing file
  out<< "\\end{tabular}"<<endl;
  out << "}\n"; // For resizebox
  out.close();

  //// Copying header and table to the compilable file
  TString fullname = outname; fullname.ReplaceAll("table_","fulltable_");
  ofstream full(fullname);
  ifstream header("txt/header.tex");
  full<<header.rdbuf();
  header.close();
  if(!abcd.scenario.Contains("signal")) full << "\\usepackage[landscape]{geometry}\n\n";
  full << "\\begin{document}\n\n";
  full << "\\begin{preview}\n";
  // full << "\\caption{" << abcd.caption <<".}\\vspace{0.1in}\n\\label{tab:"<<abcd.scenario<<"}\n";
  ifstream outtab(outname);
  full << outtab.rdbuf();
  outtab.close();
  full << "\\end{preview}\n";
  full << "\\end{document}\n";
  full.close();

  return fullname;
} // printTable

//// Makes kappa plots
void plotKappa(abcd_def &abcd, vector<vector<vector<float> > > &kappas, 
               vector<vector<vector<float> > > &kappas_mm, vector<vector<vector<float> > > &kmcdat, float total_luminosity){

  double markerSize = 1.1;

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Kappa");
  if(kappas.size() >= 1) { // Used to be 4
    opts.CanvasWidth(2000);
    markerSize = 1.5;
  }
  setPlotStyle(opts);

  //// k_ordered has all the kappas, k[plane][bin][value_index]
  // value_index: 0 = central value, 1 = up value, 2 = down value
  vector<vector<vector<float> > > k_ordered, kmd_ordered, k_ordered_mm;
   
  int nbins = 0; // Total number of njet bins (used in the base histo)
  for(size_t iplane=0; iplane < kappas.size(); iplane++) {
    k_ordered.push_back(vector<vector<float> >());
    kmd_ordered.push_back(vector<vector<float> >());
    k_ordered_mm.push_back(vector<vector<float> >());
    for(size_t ibin=0; ibin < kappas[iplane].size(); ibin++){
      TString bincut = abcd.bincuts[ibin];
      bincut.ReplaceAll(" ","");
      bincut.ReplaceAll("mm_","");
      int index;
      do{
        index = bincut.First('[');
        bincut.Remove(index, bincut.First(']')-index+1);
      }while(index>=0);

      k_ordered[iplane].push_back(kappas[iplane][ibin]);
      kmd_ordered[iplane].push_back(kmcdat[iplane][ibin]);
      k_ordered_mm[iplane].push_back(kappas_mm[iplane][ibin]);
      nbins++;
    } // Loop over bin cuts
  } // Loop over plane cuts

  //// Plotting kappas
  TCanvas can("can","");
  can.SetFillStyle(4000);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex klabel; klabel.SetTextFont(42); klabel.SetTextAlign(23);


  float minx = 0.5, maxx = nbins+0.5, miny = 0;
  float maxy = 3;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  TString ytitle = "#kappa";
  if(alt_scen!="data" && alt_scen!="mc") ytitle += " (Scen. = "+alt_scen+")";
  histo.SetTitleOffset(0.55,"y");
  histo.SetTitleSize(0.07,"y");
  histo.SetYTitle(ytitle);
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  int bin = 0;
  vector<double> vx, vexh, vexl, vy, veyh, veyl;
  vector<double> vx_mm, vexh_mm, vexl_mm, vy_mm, veyh_mm, veyl_mm;
  vector<double> vx_kmd, vexh_kmd, vexl_kmd, vy_kmd, veyh_kmd, veyl_kmd;
  for(size_t iplane=0; iplane < k_ordered.size(); iplane++) {
    for(size_t ibin=0; ibin < k_ordered[iplane].size(); ibin++){
      bin++;
      histo.GetXaxis()->SetBinLabel(bin, "");
      vx.push_back(bin);
      vexl.push_back(0);
      vexh.push_back(0);
      vy.push_back(k_ordered[iplane][ibin][0]);
      veyh.push_back(k_ordered[iplane][ibin][1]);
      veyl.push_back(k_ordered[iplane][ibin][2]);

      //// MC kappas with data uncertainties
      vx_kmd.push_back(bin);
      vexl_kmd.push_back(0);
      vexh_kmd.push_back(0);
      vy_kmd.push_back(kmd_ordered[iplane][ibin][0]);
      float ekmdUp = sqrt(pow(k_ordered[iplane][ibin][1],2) +
                          pow(kmd_ordered[iplane][ibin][1],2));
      float ekmdDown = sqrt(pow(k_ordered[iplane][ibin][2],2) +
                            pow(kmd_ordered[iplane][ibin][2],2));
      veyh_kmd.push_back(ekmdUp);            
      veyl_kmd.push_back(ekmdDown);

      //// Data/pseudodata kappas
      vx_mm.push_back(bin+0.1);
      vexl_mm.push_back(0);
      vexh_mm.push_back(0);
      vy_mm.push_back(k_ordered_mm[iplane][ibin][0]);
      if(alt_scen=="data" || alt_scen=="mc_as_data") {
        veyh_mm.push_back(k_ordered_mm[iplane][ibin][1]);
        veyl_mm.push_back(k_ordered_mm[iplane][ibin][2]);         
      } else {     
        veyh_mm.push_back(0);
        veyl_mm.push_back(0);
      }

      // Printing difference between kappa and kappa_mm
      float kap = k_ordered[iplane][ibin][0], kap_mm = k_ordered_mm[iplane][ibin][0];
      TString text = ""; 
      //if(alt_scen=="data") text = "#Delta_{#kappa} = "+RoundNumber((kap_mm-1)*100,0,1)+"%"; // with respect to kappa being 1
      if(alt_scen=="data") text = "#Delta_{#kappa} = "+RoundNumber((kap_mm-kap)*100,0,kap)+"%"; // with respect to kap because mc is value that will be used.
      else if (alt_scen=="mc_as_data" || alt_scen=="mc") text = "#Delta_{#kappa}="+RoundNumber((kap-1)*100,0,1)+"%";
      else {
        /*if fake mismeasure*/ text = "#Delta_{#kappa}="+RoundNumber((kap_mm-kap)*100,0,kap)+"%";
        cout<<"bin["<<ibin<<"] plane["<<iplane<<"] Delta kappa: "<<RoundNumber((kap_mm-kap)*100,0,kap)+"%"<<endl;
      }
      klabel.SetTextSize(abcd.planecuts.size()>=10 ? 0.025 : 0.035);
      klabel.DrawLatex(bin, 0.85*maxy, text);
      //// Printing stat uncertainty of kappa_mm/kappa
      float kapUp = k_ordered[iplane][ibin][1], kapDown = k_ordered[iplane][ibin][2];
      float kap_mmUp = k_ordered_mm[iplane][ibin][1];
      float kap_mmDown = k_ordered_mm[iplane][ibin][2];
      if(alt_scen=="data") {
        //TString unc_ = RoundNumber(kap_mmUp*100,0, 1)>RoundNumber(kap_mmDown*100,0, 1) ? RoundNumber(kap_mmUp*100,0, 1) : RoundNumber(kap_mmDown*100,0, 1);
        //TString unc_ = RoundNumber(kap_mmUp*100,0, kap_mm)>RoundNumber(kap_mmDown*100,0, kap_mm) ? RoundNumber(kap_mmUp*100,0, kap_mm) : RoundNumber(kap_mmDown*100,0, kap_mm);
        float kapUp_sum = sqrt(pow(kapUp,2)+pow(kap_mmUp,2));
        float kapDown_sum = sqrt(pow(kapDown,2)+pow(kap_mmDown,2));
        TString unc_ = RoundNumber(kapUp_sum*100,0, kap_mm)>RoundNumber(kapDown_sum*100,0, kap_mm) ? RoundNumber(kapUp_sum*100,0, kap_mm) : RoundNumber(kapDown_sum*100,0, kap_mm);
        // text = "#sigma_{stat}=^{+"+RoundNumber(kap_mmUp*100,0, 1)+"%}_{-"+RoundNumber(kap_mmDown*100,0, 1)+"%}";
        text = "#sigma_{st}^{#Delta}="+unc_+"%";
      } else if (alt_scen == "mc_as_data") {
        // For Delta uncertainty
        //float kapUp_sum = sqrt(pow(kapUp,2)+pow(kap_mmUp,2));
        //float kapDown_sum = sqrt(pow(kapDown,2)+pow(kap_mmDown,2));
        //TString unc_ = RoundNumber(kapUp_sum*100,0, kap_mm)>RoundNumber(kapDown_sum*100,0, kap_mm) ? RoundNumber(kapUp_sum*100,0, kap_mm) : RoundNumber(kapDown_sum*100,0, kap_mm);
        // For MC uncertainty
        //TString unc_ = RoundNumber(kapUp*100,0, kap)>RoundNumber(kapDown*100,0, kap) ? RoundNumber(kapUp*100,0, kap) : RoundNumber(kapDown*100,0, kap);
        // For data uncertainty
        TString unc_ = RoundNumber(kap_mmUp*100,0, kap_mm)>RoundNumber(kap_mmDown*100,0, kap_mm) ? RoundNumber(kap_mmUp*100,0, kap_mm) : RoundNumber(kap_mmDown*100,0, kap_mm);
        text = "#sigma_{st}^{dat}="+unc_+"%";
        cout<<"bin["<<ibin<<"] plane["<<iplane<<"] unc kappa: "<<unc_+"%"<<endl;
      } else {
        TString unc_ = RoundNumber(kapUp*100,0, kap)>RoundNumber(kapDown*100,0, kap) ? RoundNumber(kapUp*100,0, kap) : RoundNumber(kapDown*100,0, kap);
        // text = "#sigma_{stat}=^{+"+RoundNumber(kapUp*100,0, 1)+"%}_{-"+RoundNumber(kapDown*100,0, 1)+"%}";
        text = "#sigma_{st}^{mc}="+unc_+"%";
      }
      klabel.DrawLatex(bin, 0.78*maxy, text);
      // adding label to indicate the ABCD corresponding to each kappa value
      klabel.SetTextSize(abcd.planecuts.size()>=10 ? 0.025 : 0.035);
      
      if (sample=="search" || sample=="ttbar") text = ibin%2==0 ? "3b/2b" : "4b/2b"; 
      else if (sample=="zll") text = do_midnb ? "2b/1b" : "1b/0b"; 
      else if (sample=="qcd") {
        if (do_highnb) text =  ibin%2==0 ? "3b/2b" : "4b/2b"; 
        else text = do_midnb ? "2b/1b" : "1b/0b"; 
      }
      klabel.DrawLatex(bin, 0.92*maxy, text);
    } // Loop over bin cuts

    // Drawing lines separating MET planes
    line.SetLineStyle(iplane%2==0 ? 2:1); line.SetLineWidth(1); line.SetLineColor(kBlack);
    if (iplane<k_ordered.size()-1) line.DrawLine(bin+0.5, miny, bin+0.5, maxy);

    // All the labels on the X-axis...
    TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);
    label.SetTextSize(abcd.planecuts.size()>=10 ? 0.025 : 0.045);
    string metdef = "met";
    if (sample=="zll") metdef = "ll_pt";
    double lmargin(opts.LeftMargin()), rmargin(opts.RightMargin()), bmargin(opts.BottomMargin());
    if (abcd.planecuts[iplane].Contains("drmax")) {
      //separate the MET and dRmax cuts
      string plabel = abcd.planecuts[iplane].Data();
      string drlabel = plabel.substr(plabel.find("hig_cand_drmax"), plabel.length());
      drlabel = CodeToRootTex(drlabel);
      label.DrawLatexNDC(lmargin+(1-rmargin-lmargin)/abcd.planecuts.size()*(iplane+0.5), 0.25, drlabel.c_str());

      if (iplane%2==0) { // write the MET bin only once per each pair of drmax bins
        ReplaceAll(plabel," ", "");
        string metlabel = plabel.substr(0, plabel.find("&&hig_cand_drmax"));
        metlabel = CodeToRootTex(metlabel); 
        ReplaceAll(metlabel, "<"+metdef+"#leq","#minus"); 
        if(metlabel == "> 0") metlabel = "Inclusive";
        label.DrawLatexNDC(lmargin+(1-rmargin-lmargin)/abcd.planecuts.size()*(iplane+1), 0.13, metlabel.c_str());
      }
    } else { // if not binning in dRmax
      string plabel = CodeToRootTex(abcd.planecuts[iplane].Data());
      ReplaceAll(plabel, "<"+metdef+"#leq","#minus"); 
      if(plabel == "> 0") plabel = "Inclusive";
      label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., -0.03*maxy, plabel.c_str());
    }
    // write pT miss at the bottom of the plot...
    if (iplane==0) {
      if (sample=="zll") label.DrawLatexNDC(lmargin+(1-rmargin-lmargin)/2., bmargin-0.08, "p_{T}^{Z} [GeV]");
      else label.DrawLatexNDC(lmargin+(1-rmargin-lmargin)/2., bmargin-0.08, "p_{T}^{miss} [GeV]");
    }
  } // Loop over plane cuts

  //// Drawing legend and TGraphs
  if (debug) cout<<"Building up TGraphs"<<endl;
  double legX(0.595), legY(1-0.035), legH(0.05), legW(0.15);
  if(alt_scen != "data") legX = 0.65;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()*1.15); leg.SetFillColor(0);leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(2);
  TGraphAsymmErrors graph_kmd = TGraphAsymmErrors(vx_kmd.size(), &(vx_kmd[0]), &(vy_kmd[0]),
                                      &(vexl_kmd[0]), &(vexh_kmd[0]), &(veyl_kmd[0]), 
                                      &(veyh_kmd[0]));
  graph_kmd.SetMarkerStyle(21); graph_kmd.SetMarkerSize(markerSize); graph_kmd.SetMarkerColor(2);
  graph_kmd.SetLineColor(1); graph_kmd.SetLineWidth(2);
  if(alt_scen=="mc_as_data") graph_kmd.Draw("p0 same");

  TGraphAsymmErrors graph = TGraphAsymmErrors(vx.size(), &(vx[0]), &(vy[0]),
                                  &(vexl[0]), &(vexh[0]), &(veyl[0]), &(veyh[0]));
  graph.SetMarkerStyle(21); graph.SetMarkerSize(markerSize); graph.SetMarkerColor(2);
  graph.SetLineColor(2); graph.SetLineWidth(2);
  graph.Draw("p0 same");

  TGraphAsymmErrors graph_mm = TGraphAsymmErrors(vx_mm.size(), &(vx_mm[0]), &(vy_mm[0]),
                                     &(vexl_mm[0]), &(vexh_mm[0]), &(veyl_mm[0]), 
                                     &(veyh_mm[0]));
  graph_mm.SetMarkerStyle(20); graph_mm.SetMarkerSize(markerSize*1.2); graph_mm.SetMarkerColor(1);
  graph_mm.SetLineColor(1); graph_mm.SetLineWidth(2);
  if(alt_scen!="mc_as_data" && alt_scen!="mc") graph_mm.Draw("p0 same");

  leg.AddEntry(&graph, "MC", "ep");
  TString data_s = (alt_scen=="data"?"Data":"Pseudodata");
  if(alt_scen!="mc_as_data" && alt_scen!="mc") 
    leg.AddEntry(&graph_mm, data_s, "ep");

  leg.Draw();

  //// Drawing CMS labels and line at 1
  TString cmsPrel = "#font[62]{CMS} #scale[0.8]{#font[52]{}}";
  TString cmsSim = "#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}";
  TLatex cmslabel;
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  if(alt_scen != "data") cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,cmsSim);
  else cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,cmsPrel);
  cmslabel.SetTextAlign(31);
  //cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");
  cmslabel.SetTextSize(0.053);

  ///// Luminosity and energy
  TString abcd_title;
  if(sample.Contains("zll")) abcd_title = "Dilepton control region";
  else if(sample.Contains("qcd")) abcd_title = "Low #Delta#phi control region";
  else if(sample.Contains("ttbar")) abcd_title = "Single-lepton control region";
  else abcd_title = "Search region";

  TString title = "";
  TString fontstyle = RoundNumber(opts.Font()+10,0);
  if(alt_scen!="mc_as_data") title = "#font[42]{"+RoundNumber(total_luminosity, 0)+" fb^{-1} (13 TeV)}";
  if(sample=="search" && alt_scen != "data") title = "#font["+fontstyle+"]{"+abcd_title+"}";
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, title);

  ///// sample name
  cmslabel.SetTextAlign(11);
  title = "#font["+fontstyle+"]{"+abcd_title+"}";
  if(!(sample=="search" && alt_scen != "data")) cmslabel.DrawLatex(opts.LeftMargin()+0.2, 1-opts.TopMargin()+0.015, title);

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname="plots/kappa_"+sample+"_" +abcd.scenario + (do_highnb?"_highnb":"")+(do_midnb?"_midnb":"");
  fname += "_lumi"+RoundNumber(lumi, 0)+".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl; 

}

TH1D * getTH1DFromAllyields(size_t yieldIndex, vector<vector<GammaParams> > &allyields, abcd_def &abcd, int planeType = -1, bool blindSignalRegion=false) {
  // planeType = -1: all, 1: low-drmax, 2: high-drmax
  int nUnrolledBins;
  if (planeType==-1) nUnrolledBins = abcd.planecuts.size()*(2+2*abcd.bincuts.size());
  else nUnrolledBins = abcd.planecuts.size()*(2+2*abcd.bincuts.size())/2;
  string unrolledHistName = "unrolledHist"+string(RoundNumber(yieldIndex,1).Data());
  TH1D * unrolledHist = new TH1D(unrolledHistName.c_str(), unrolledHistName.c_str(), nUnrolledBins, -0.5, nUnrolledBins-0.5);
  vector<TString> binNames({"SBD, crb", "SBD, xb", "HIG, crb", "HIG, xb"});
  int iUnrollBin = 1;
  // plane is drmax-met binning
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    if (planeType == 1 && iplane%2 == 1) continue;
    if (planeType == 2 && iplane%2 == 0) continue;
    //cout<<"iplane: "<<iplane<<endl;
    // ibin is 3b, 4b binning
    for(size_t ibin=0; ibin < abcd.bincuts.size(); ibin++){
      // 0: 2b_SDB, 1: 2b_SDB, 2: 2b_HIG, 3: xb_HIG
      for(auto iabcd:{0,2,1,3}){
        if (ibin>0 && iabcd%2==0) continue; // don't loop the 2b regions again
        size_t index = iplane*abcd.bincuts.size()*4+ibin*4+iabcd;
        // Get binName
        TString binName = binNames[iabcd];
        if (!do_highnb && (sample=="zll" || sample=="qcd")) {
          if(do_midnb) binName.ReplaceAll("xb", "2b").ReplaceAll("crb", "1b");
          else binName.ReplaceAll("xb", "1b").ReplaceAll("crb", "0b");
        } else {
          if(ibin==0) binName.ReplaceAll("xb", "3b");
          else binName.ReplaceAll("xb", "4b");
          binName.ReplaceAll("crb", "2b");
        }
        // Get Yields
        //cout<<"binName: "<<binName<<" "<<RoundNumber(allyields[yieldIndex][index].Yield(),1)<<endl;
        if (blindSignalRegion && iabcd==3) {
          // Do not unblind
        } else {
          unrolledHist->SetBinContent(iUnrollBin, allyields[yieldIndex][index].Yield());
          unrolledHist->SetBinError(iUnrollBin, allyields[yieldIndex][index].Uncertainty());
        }
        iUnrollBin++;
      }
    }
  }
  //TCanvas cTest("cTest", "cTest", 500, 500);
  //unrolledHist->Draw();
  //cTest.SaveAs("cTest.pdf");
  return unrolledHist;
}

TGraphAsymmErrors * getTGraphAsymmErrorsFromPreds(vector<vector<vector<float> > > & preds, abcd_def &abcd, int planeType = -1) {
  string unrolledPredHistName = "unrolledPredHist";
  //TH1D * unrolledPredHist = new TH1D(unrolledPredHistName.c_str(), unrolledPredHistName.c_str(), nUnrolledBins, -0.5, nUnrolledBins-0.5);
  int iUnrollBin = 0;
  vector<Double_t> xpoints;
  vector<Double_t> ypoints;
  vector<Double_t> xerror_low;
  vector<Double_t> xerror_high;
  vector<Double_t> yerror_low;
  vector<Double_t> yerror_high;
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    if (planeType == 1 && iplane%2 == 1) continue;
    if (planeType == 2 && iplane%2 == 0) continue;
    for(size_t ibin=0; ibin < abcd.bincuts.size(); ibin++){
      // 0: 2b_SDB, 1: 2b_SDB, 2: 2b_HIG, 3: xb_HIG
      for(auto iabcd:{0,2,1,3}){
        if (ibin>0 && iabcd%2==0) continue; // don't loop the 2b regions again
        if (iabcd==3) {
          //cout<<"iPlane: "<<iplane<<" iBin: "<<ibin<<" iabcd: "<<iabcd<<" iUnrollBin: "<<iUnrollBin<<" preds: "<<preds[iplane][ibin][0]<<" "<<preds[iplane][ibin][1]<<" "<<preds[iplane][ibin][2]<<endl;
          double errorHigh = preds[iplane][ibin][1];
          double errorLow = preds[iplane][ibin][2];
          double prediction = preds[iplane][ibin][0];
          xpoints.push_back(iUnrollBin);
          ypoints.push_back(prediction);
          xerror_low.push_back(0);
          xerror_high.push_back(0);
          yerror_low.push_back(errorLow);
          yerror_high.push_back(errorHigh);
        }
        iUnrollBin++;
      }
    }
  }
  TGraphAsymmErrors * unrolledPredHist = new TGraphAsymmErrors(xpoints.size(), &xpoints[0], &ypoints[0], &xerror_low[0], &xerror_high[0], &yerror_low[0], &yerror_high[0]);
  return unrolledPredHist;
}

void changeToPoissonUncertainty(TH1D * hist) {
  int nBins = hist->GetNbinsX();
  for (int iBin = 1; iBin < nBins+1; iBin++) {
    //double oldUncertainty = hist->GetBinError(iBin);
    double yield = hist->GetBinContent(iBin);
    hist->SetBinError(iBin, sqrt(yield));
  }
}

void changeHistToMinimum(TH1D * hist, double min) {
  int nBins = hist->GetNbinsX();
  for (int iBin = 1; iBin < nBins+1; iBin++) {
    //double oldUncertainty = hist->GetBinError(iBin);
    double yield = hist->GetBinContent(iBin);
    if (yield < min) hist->SetBinContent(iBin, min);
  }
}

void changePredToMinimum(TGraphAsymmErrors * graph, double min) {
  int nPoints = graph->GetN();
  for (int iPoint = 0; iPoint < nPoints; iPoint++) {
    //double oldUncertainty = hist->GetBinError(iBin);
    Double_t xPoint, yPoint;
    graph->GetPoint(iPoint, xPoint, yPoint);
    Double_t yerror_low = graph->GetErrorYlow(iPoint);
    Double_t yerror_high = graph->GetErrorYhigh(iPoint);
    if (yPoint < min) {
      graph->SetPoint(iPoint, xPoint, min);
      Double_t new_yerror_low = yPoint-yerror_low - min;
      if (new_yerror_low<0) new_yerror_low = 0;
      Double_t new_yerror_high = yPoint+yerror_high - min;
      if (new_yerror_high<0) new_yerror_high = 0;
      graph->SetPointError(iPoint, 0, 0, new_yerror_low, new_yerror_high);
    }
  }
}

void plotTable(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas __attribute__((unused)), vector<vector<vector<float> > > &preds, 
                   vector<vector<float> > __attribute__((unused)) yieldsPlane, vector<shared_ptr<Process> > __attribute__((unused)) &proc_sigs, float __attribute__((unused)) total_luminosity){
  // allyields [0] is data, [1] is bkg, [2] is sig1, [3] is sig2, [4] is sig3, [5] is other, [6] is vjet, [7] is ttx
  TH1D * bkgHistLow = getTH1DFromAllyields(/*yieldIndex*/1, allyields, abcd, /*planeType*/1);
  changeToPoissonUncertainty(bkgHistLow);

  TH1D * otherHistLow = getTH1DFromAllyields(/*yieldIndex*/5, allyields, abcd, /*planeType*/1);
  TH1D * vjetHistLow = getTH1DFromAllyields(/*yieldIndex*/6, allyields, abcd, /*planeType*/1);
  TH1D * ttxHistLow = getTH1DFromAllyields(/*yieldIndex*/7, allyields, abcd, /*planeType*/1);
  THStack * bkgStackHistLow = new THStack("bkgStackLow", "Low #Delta R_{max};Bin;Yield");
  otherHistLow->SetFillColor(kGray);
  bkgStackHistLow->Add(otherHistLow);
  vjetHistLow->SetFillColor(kOrange-4);
  bkgStackHistLow->Add(vjetHistLow);
  ttxHistLow->SetFillColor(kBlue-4);
  bkgStackHistLow->Add(ttxHistLow);
  bkgStackHistLow->SetMinimum(0.1);

  TH1D * sig2HistLow = getTH1DFromAllyields(/*yieldIndex*/3, allyields, abcd, /*planeType*/1);
  sig2HistLow->SetLineColor(kPink-4);
  sig2HistLow->SetLineWidth(2);

  TGraphAsymmErrors * predHistLow = getTGraphAsymmErrorsFromPreds(preds, abcd, /*planeType*/1);
  predHistLow->SetMarkerStyle(21);
  //predHistLow->SetMarkerSize(2);
  predHistLow->SetMarkerColor(kRed);
  predHistLow->SetLineColor(kRed);
  predHistLow->SetLineWidth(3);
  changePredToMinimum(predHistLow, 0.03);

  TH1D * dataHistLow = getTH1DFromAllyields(/*yieldIndex*/0, allyields, abcd, /*planeType*/1, /*blindSignalRegion*/!unblind);
  dataHistLow->SetMarkerStyle(kFullCircle);
  changeHistToMinimum(dataHistLow, 0.03);

  ////
  TH1D * bkgHistHigh = getTH1DFromAllyields(/*yieldIndex*/1, allyields, abcd, /*planeType*/2);
  changeToPoissonUncertainty(bkgHistHigh);

  TH1D * otherHistHigh = getTH1DFromAllyields(/*yieldIndex*/5, allyields, abcd, /*planeType*/2);
  TH1D * vjetHistHigh = getTH1DFromAllyields(/*yieldIndex*/6, allyields, abcd, /*planeType*/2);
  TH1D * ttxHistHigh = getTH1DFromAllyields(/*yieldIndex*/7, allyields, abcd, /*planeType*/2);
  THStack * bkgStackHistHigh = new THStack("bkgStackHigh", "High #Delta R_{max};Bin;Yield");
  otherHistHigh->SetFillColor(kGray);
  bkgStackHistHigh->Add(otherHistHigh);
  vjetHistHigh->SetFillColor(kOrange-4);
  bkgStackHistHigh->Add(vjetHistHigh);
  ttxHistHigh->SetFillColor(kBlue-4);
  bkgStackHistHigh->Add(ttxHistHigh);
  bkgStackHistHigh->SetMinimum(0.1);

  TH1D * sig2HistHigh = getTH1DFromAllyields(/*yieldIndex*/3, allyields, abcd, /*planeType*/2);
  sig2HistHigh->SetLineColor(kPink-4);
  sig2HistHigh->SetLineWidth(2);

  TGraphAsymmErrors * predHistHigh = getTGraphAsymmErrorsFromPreds(preds, abcd, /*planeType*/2);
  predHistHigh->SetMarkerStyle(21);
  //predHistHigh->SetMarkerSize(2);
  predHistHigh->SetMarkerColor(kRed);
  predHistHigh->SetLineColor(kRed);
  predHistHigh->SetLineWidth(3);
  changePredToMinimum(predHistHigh, 0.03);

  TH1D * dataHistHigh = getTH1DFromAllyields(/*yieldIndex*/0, allyields, abcd, /*planeType*/2, /*blindSignalRegion*/!unblind);
  dataHistHigh->SetMarkerStyle(kFullCircle);
  changeHistToMinimum(dataHistHigh, 0.04);


  //TH1D * sig1Hist = getTH1DFromAllyields(/*yieldIndex*/2, allyields, abcd, /*planeType*/1);
  //sig1Hist->SetLineColor(kRed);
  //sig1Hist->SetLineWidth(2);
  //TH1D * sig2Hist = getTH1DFromAllyields(/*yieldIndex*/3, allyields, abcd, /*planeType*/1);
  //sig2Hist->SetLineColor(kGreen);
  //sig2Hist->SetLineWidth(2);
  //TH1D * sig3Hist = getTH1DFromAllyields(/*yieldIndex*/4, allyields, abcd, /*planeType*/1);
  //sig3Hist->SetLineColor(kBlue);
  //sig3Hist->SetLineWidth(2);
  //gStyle->SetErrorX(0);

  TCanvas cTest("cTest", "cTest", 500, 500);
  TPad *pad1 = new TPad("top", "top", 0.0, 0.5, 1.0, 1.0);
  TPad *pad2 = new TPad("bot", "bot", 0.0, 0.0, 1.0, 0.5);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  bkgStackHistLow->Draw("hist");
  bkgHistLow->SetFillColor(kBlack);
  bkgHistLow->SetFillStyle(3013);
  bkgHistLow->Draw("same E2");
  cTest.Update();
  for (unsigned iPlane = 0; iPlane < 4; iPlane++) {
    float binIndex = iPlane*6;
    TLine * line = new TLine(binIndex-0.5,0,binIndex-0.5,pad1->GetUymax());
    line->SetLineWidth(3);
    line->Draw("same");
  }
  gPad->SetLogy();
  //sig1Hist->Draw("same hist");
  //sig2HistLow->Draw("same hist");
  //sig3Hist->Draw("same hist");
  //predHistLow->Draw("same e0 x0");
  predHistLow->Draw("same p");
  dataHistLow->Draw("same hist p");

  pad2->cd();
  bkgStackHistHigh->Draw("hist");
  bkgHistHigh->SetFillColor(kBlack);
  bkgHistHigh->SetFillStyle(3013);
  bkgHistHigh->Draw("same E2");
  cTest.Update();
  for (unsigned iPlane = 0; iPlane < 4; iPlane++) {
    float binIndex = iPlane*6;
    TLine * line = new TLine(binIndex-0.5,0,binIndex-0.5,pad2->GetUymax());
    line->SetLineWidth(3);
    line->Draw("same");
  }
  gPad->SetLogy();
  //sig2HistHigh->Draw("same hist");
  //predHistHigh->Draw("same e0 x0");
  predHistHigh->Draw("same p");
  dataHistHigh->Draw("same hist p");
  //TH1D * ratioHist = static_cast<TH1D*>(dataHist->Clone("data_over_mc"));
  //ratioHist->Divide(bkgHist);
  //ratioHist->SetMaximum(2);
  //ratioHist->SetMinimum(0);
  //ratioHist->Draw("pe");
  cTest.SaveAs("plots/unrolled_regions.pdf");
  cout<<"open plots/unrolled_regions.pdf"<<endl;
  return;
}

//// Calculating kappa and Total bkg prediction
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
vector<vector<float> > findPreds(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
               vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &kappas_mm, 
               vector<vector<vector<float> > > &kmcdat, vector<vector<vector<float> > > &preds, bool saveData){

  stringstream dataOutput;

  // Powers for kappa:   ({R1, R2, D3, R4})
  vector<float> pow_kappa({ 1, -1, -1,  1});
  // Powers for TotBkg pred:({R1, R2, D3,  R1, R2, D3, D4})
  vector<float> pow_totpred( {-1,  1,  1,   1, -1, -1,  1});
  vector<float> pow_datapred( {-1,  1,  1});

  //// total normalization
  //GammaParams Ndata_all, Nmc_all;
  //for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
  //  //// Counting yields in plane so that we can normalize the total yield in the ABCD in MC to that in the data
  //  for(size_t ibin=0; ibin < abcd.bincuts.size(); ibin++){
  //    for(size_t iabcd=0; iabcd < 4; iabcd++){
  //      size_t index = iplane*abcd.bincuts.size()*4 +ibin*4 + iabcd;
  //      Ndata_all += allyields[0][index];
  //      Nmc_all += allyields[1][index];
  //    } // Loop over ABCD cuts
  //  } // Loop over bin cuts
  //}
  //float dataMC_all = Ndata_all.Yield()/Nmc_all.Yield();

  float val(1.), valup(1.), valdown(1.);
  vector<vector<float> > yieldsPlane;
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    //// Counting yields in plane so that we can normalize the total yield in the ABCD in MC to that in the data
    GammaParams NdataPlane, NmcPlane;
    for(size_t ibin=0; ibin < abcd.bincuts.size(); ibin++){
      for(size_t iabcd=0; iabcd < 4; iabcd++){
        size_t index = iplane*abcd.bincuts.size()*4 +ibin*4 + iabcd;
        NdataPlane += allyields[0][index];
        NmcPlane += allyields[1][index];
      } // Loop over ABCD cuts
    } // Loop over bin cuts
    //cout<<"Plane "<<iplane<<": MC is "<<NmcPlane<<", data is "<<NdataPlane<<endl;
    float Nobs = NdataPlane.Yield(), Nmc = NmcPlane.Yield();
    float dataMC = Nobs/Nmc; // Different for each bin
    float edataMC = sqrt(pow(sqrt(Nobs)/Nmc,2) + pow(Nobs*NmcPlane.Uncertainty()/Nmc/Nmc,2));
    yieldsPlane.push_back({dataMC, edataMC});

    kappas.push_back(vector<vector<float> >());
    kmcdat.push_back(vector<vector<float> >());
    kappas_mm.push_back(vector<vector<float> >());
    preds.push_back(vector<vector<float> >());
    for(size_t ibin=0; ibin < abcd.bincuts.size(); ibin++){
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      //// Pushing data yields for predictions
      for(size_t iabcd=0; iabcd < 3; iabcd++){
        size_t index = iplane*abcd.bincuts.size()*4 +ibin*4 + iabcd;
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[0][index].Yield());
        weights.back().push_back(1.);
      } // Loop over ABCD cuts

      // Throwing toys to find predictions with no kappa -> used for data stat. unc. in systematics table
      val = calcKappa(entries, weights, pow_datapred, valdown, valup);
      if(valdown<0) valdown = 0;

      vector<vector<float> > kentries;
      vector<vector<float> > kweights;
      vector<vector<float> > kkentries;
      vector<vector<float> > kkweights;
      vector<vector<float> > kentries_mm;
      vector<vector<float> > kweights_mm;
      //// Pushing MC yields for predictions and kappas
      for(size_t iabcd=0; iabcd < 4; iabcd++){
        size_t index = iplane*abcd.bincuts.size()*4 +ibin*4 + iabcd;
        //// Renormalizing MC to data
        //allyields[1][index] *= dataMC;
        //if(split_bkg) {
        //  size_t offset = (do_signal?Nsig:0);
        //  allyields[offset+1][index] *= dataMC;
        //  allyields[offset+2][index] *= dataMC;
        //  allyields[offset+3][index] *= dataMC;
        //}
        //cout<<"dataMC: "<<dataMC<<" allyields[1]["<<index<<"]: "<<allyields[1][index]<<endl;
        //// Renormalizing MC to data using global normalization
        //allyields[1][index] *= dataMC_all;
        //if(split_bkg) {
        //  size_t offset = (do_signal?Nsig:0);
        //  allyields[offset+1][index] *= dataMC_all;
        //  allyields[offset+2][index] *= dataMC_all;
        //  allyields[offset+3][index] *= dataMC_all;
        //}

        // Yields for predictions
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[1][index].NEffective());
        weights.back().push_back(allyields[1][index].Weight());
        // Yields for kappas
        kentries.push_back(vector<float>());
        kweights.push_back(vector<float>());
        kentries.back().push_back(allyields[1][index].NEffective());
        kweights.back().push_back(allyields[1][index].Weight());
        // Yields for kappas on pseudodata
        kentries_mm.push_back(vector<float>());
        kweights_mm.push_back(vector<float>());
        kentries_mm.back().push_back(allyields[0][index].Yield());
        kweights_mm.back().push_back(1.);
        // Yields for kappas_mc normalized to data
        kkentries.push_back(vector<float>());
        kkweights.push_back(vector<float>());
        kkentries.back().push_back(allyields[1][index].Yield());
        kkweights.back().push_back(1.);

      } // Loop over ABCD cuts
      bool print_entries_weights = HigUtilities::is_in_string_options(string_options, "print_entries_weights");
      if (print_entries_weights) cout<<"iplane: "<<iplane<<endl;
      dataOutput << "iplane: "<<iplane<<" ibin: "<<ibin<<endl;
      // Throwing toys to find predictions and uncertainties
      val = calcKappa(entries, weights, pow_totpred, valdown, valup);
      if (print_entries_weights) {
        cout<<"Prediction"<<endl;
        for (unsigned index = 0; index < 7; ++index) cout<<"entries["<<index<<"]: "<<entries[index][0]<<" ";
        cout<<endl;
        for (unsigned index = 0; index < 7; ++index) cout<<"weights["<<index<<"]: "<<weights[index][0]<<" ";
        cout<<endl;
        cout<<"val: "<<val<<" "<<valdown<<" "<<valup<<endl;
      }
      //bool do_data=false; bool verbose=false; double syst=-1; bool do_plot=true; int nrep=100000; int nSigma=1;
      //val = calcKappa(entries, weights, pow_totpred, valdown, valup, do_data, verbose, syst, do_plot, nrep, nSigma);
      if(valdown<0) valdown = 0;
      preds[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries, kweights, pow_kappa, valdown, valup);
      //bool do_data=false; bool verbose=false; double syst=-1; bool do_plot=true; int nrep=100000; int nSigma=1;
      //val = calcKappa(kentries, kweights, pow_kappa, valdown, valup, do_data, verbose, syst, do_plot, nrep, nSigma);
      if(valdown<0) valdown = 0;
      kappas[iplane].push_back(vector<float>({val, valup, valdown}));
      if (print_entries_weights) {
        cout<<"MC        : "<<kentries[0][0]<<" "<<kentries[1][0]<<" "<<kentries[2][0]<<" "<<kentries[3][0]<<endl;
        cout<<"MC weight : "<<kweights[0][0]<<" "<<kweights[1][0]<<" "<<kweights[2][0]<<" "<<kweights[3][0]<<endl;
        cout<<"MC kappa  : "<<val<<" "<<valup<<" "<<valdown<<endl;
      }
      dataOutput << "MC eff. entries: "<<kentries[0][0]<<" "<<kentries[1][0]<<" "<<kentries[2][0]<<" "<<kentries[3][0]<<endl;
      dataOutput << "MC eff. weights: "<<kweights[0][0]<<" "<<kweights[1][0]<<" "<<kweights[2][0]<<" "<<kweights[3][0]<<endl;
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries_mm, kweights_mm, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kappas_mm[iplane].push_back(vector<float>({val, valup, valdown}));
      if (print_entries_weights) {
        cout<<"Data        : "<<kentries_mm[0][0]<<" "<<kentries_mm[1][0]<<" "<<kentries_mm[2][0]<<" "<<kentries_mm[3][0]<<endl;
        cout<<"Data weight : "<<kweights_mm[0][0]<<" "<<kweights_mm[1][0]<<" "<<kweights_mm[2][0]<<" "<<kweights_mm[3][0]<<endl;
        cout<<"Data kappa  : "<<val<<" "<<valup<<" "<<valdown<<endl;
      }
      dataOutput << "Data entries: "<<kentries_mm[0][0]<<" "<<kentries_mm[1][0]<<" "<<kentries_mm[2][0]<<" "<<kentries_mm[3][0]<<endl;
      dataOutput << "Data weights: "<<kweights_mm[0][0]<<" "<<kweights_mm[1][0]<<" "<<kweights_mm[2][0]<<" "<<kweights_mm[3][0]<<endl;
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kkentries, kkweights, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kmcdat[iplane].push_back(vector<float>({val, valup, valdown}));
    } // Loop over bin cuts
  } // Loop over plane cuts

  if (saveData) {
    string data_filename = "counts_in_bins.txt";
    ofstream data_file(data_filename);
    data_file<<dataOutput.rdbuf()<<endl;
    data_file.close();
  }

  return yieldsPlane;
} // findPreds

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 's'},
      {"debug", no_argument, 0, 'd'},         // Debug: prints yields and cuts used
      {"quick", no_argument, 0, 'q'},         // Used inclusive ttbar for quick testing
      {"unblind", no_argument, 0, 'u'},         // Opens data or pseudodata
      {"scen", required_argument, 0, 0},       // mc, data, mc_as_data, ...
      {"midnb", no_argument, 0, 0},           // Check zll and qcd CRs for 2b
      {"highnb", no_argument, 0, 0},          // Check QCD CR at 3b and 4b
      {"incl_drmax", no_argument, 0, 0},          
      {"old_bin", no_argument, 0, 0},          
      {"year", required_argument, 0, 'y'},
      {"string_options", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:kdquo:y:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      sample = optarg;
      break;
    case 'd':
      debug = true;
      break;
    case 'q':
      quick_test = true;
      break;
    case 'u':
      unblind = true;
      break;
    case 'y':
      year_string = optarg;
      break;
    case 'o':
      string_options = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "scen"){
        alt_scen = optarg;
      }else if(optname == "highnb"){
        do_highnb = true;
      }else if(optname == "midnb"){
        do_midnb = true;
      }else if(optname == "incl_drmax"){
        do_bin_drmax = false;
      }else if(optname == "old_bin"){
        do_old_bin = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
