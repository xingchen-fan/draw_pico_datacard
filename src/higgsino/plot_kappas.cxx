// Cheatsheet:
//    To obtain plots with no mismeasurement and kappa's with expected data stats:
//         ./run/hig/plot_kappa.exe --mm mc_as_data -l 36.2
//    There are 4 possibilities for the sample, requested with option -s. These are: search, zll, qcd, ttbar
//    Option -t plots the kappas with a tighter selection, see basecuts below, e.g.
//         ./run/hig/plot_kappa.exe --mm mc_as_data -t -s zll -l 36.2

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>
#include <string>

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

using namespace std;

namespace{
  bool only_kappa = false;
  bool split_bkg = true;
  bool do_signal = true;
  bool debug = false;
  bool do_highnb = false;
  bool do_midnb = false;
  int digits_table = 1;
  TString sample = "search";
  string alt_scen = "mc_as_data";
  bool only_mc = (alt_scen != "data");
  float lumi=1.;
  bool quick_test = false;
  bool print_mc = true;
  bool do_zbi = false;
  // for office use only
  vector<TString> syst_names;
  vector<vector<float>> syst_values;
  string sys_wgts_file = "txt/sys_weights.cfg";
}

struct abcd_def{
  TString scenario;
  vector<TString> planecuts;
  vector<TString> bincuts;
};

TString printTable(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds, 
                   vector<vector<float> > yieldsPlane, vector<shared_ptr<Process> > &proc_sigs);
void plotKappa(abcd_def &abcd, vector<vector<vector<float> > >  &kappas, 
               vector<vector<vector<float> > >  &kappas_mm, vector<vector<vector<float> > >  &kmcdat);
vector<vector<float> > findPreds(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
                                 vector<vector<vector<float> > > &kappas, 
                                 vector<vector<vector<float> > > &kappas_mm, 
                                 vector<vector<vector<float> > > &kmcdat, 
                                 vector<vector<vector<float> > > &datapreds,
                                 vector<vector<vector<float> > > &preds);

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  Palette colors("txt/colors.txt", "default");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// What ntuples to use?  //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<int> years = {2016};
  // years = {2016, 2017, 2018};

  string base_dir(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/");
  string mc_skim_dir("mc/merged_higmc_higtight/"), data_skim_dir("mc/merged_higdata_higloose/");
  if (sample=="ttbar")    {mc_skim_dir = "mc/merged_higmc_higlep1/"; data_skim_dir = "merged_higdata_higlep1/";} 
  else if (sample=="zll") {mc_skim_dir = "mc/merged_higmc_higlep2/"; data_skim_dir = "merged_higdata_higlep2/";} 
  else if (sample=="qcd") {mc_skim_dir = "mc/merged_higmc_higqcd/";  data_skim_dir = "merged_higdata_higqcd/";} 
  string sig_skim_dir("SMS-TChiHH_2D/merged_higmc_higtight/");

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lept*", 
                                    "*_TTZ*.root", "*_TTW*.root",
                                    "*_TTGJets*.root", "*_ttHTobb*.root","*_TTTT*.root"
                                  });
  mctags["vjets"]   = set<string>({ 
                                   "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"
                                 });
  mctags["other"]   = set<string>({
                                   // "*QCD_HT100to200_Tune*", "*QCD_HT200to300_Tune*", // these have too low weights
                                   // "*QCD_HT300to500_Tune*", 
                                   // "*QCD_HT500to700_Tune*",
                                   // "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                   // "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*",
                                   "*_ST_*.root",
                                   "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                   "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"
                                 });

  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining baseline cuts ///////////////////////////////////////
  TString c_hig_trim = "hig_cand_drmax[0]<=2.2  && hig_cand_dm[0]<=40 && hig_cand_am[0]<200";
  string baseline_s = "njet>=4 && njet<=5";
  if (sample=="search") baseline_s += " && nvlep==0 && ntk==0 && !low_dphi_met &&"+c_hig_trim;
  else if (sample=="ttbar") baseline_s += " && nlep==1 && mt<100 &&"+c_hig_trim;
  else if (sample=="zll") baseline_s += " && nlep==2 && met<50 &&"+c_hig_trim;
  else if (sample=="qcd") baseline_s += " && nvlep==0 && ntk==0 && low_dphi_met &&"+c_hig_trim;
  NamedFunc baseline = baseline_s;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////  
  // define signal processes
  vector<string> sigMasses({"400", "900"});
  vector<shared_ptr<Process> > proc_sigs;
  for (unsigned isig(0); isig<sigMasses.size(); isig++)
    proc_sigs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigMasses[isig]+",1)", 
      Process::Type::signal, 1, attach_folder(base_dir, years, sig_skim_dir,{"*TChiHH_mChi-"+sigMasses[isig]+"*.root"}), 
      baseline && "pass"));

  // define background processes
  auto proc_ttx = Process::MakeShared<Baby_pico>("tt+X", Process::Type::background, colors("tt_1l"),
    attach_folder(base_dir, years, mc_skim_dir, mctags["ttx"]), baseline && "stitch && pass");
  auto proc_vjets = Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
    attach_folder(base_dir, years, mc_skim_dir, mctags["vjets"]), baseline && "stitch && pass");
  auto proc_other = Process::MakeShared<Baby_pico>("Other", Process::Type::background, colors("other"),
    attach_folder(base_dir, years, mc_skim_dir, mctags["other"]), baseline && "stitch && pass");

  // define data or pseudo-data process
  set<string> names_data(attach_folder(base_dir, years, data_skim_dir, {"*.root"}));
  if(only_mc) { // define pseudodata if only using MC
    if (quick_test) {
      names_data = {base_dir+"/2016/"+mc_skim_dir+"/*TTJets_SingleLeptFromT_Tune*"};
    } else {
      set<string> names_allmc;
      for (auto &iset: mctags) 
        names_allmc.insert(iset.second.begin(), iset.second.end());    
      names_data = attach_folder(base_dir, years, mc_skim_dir, names_allmc);
    }
  }
  NamedFunc base_data = baseline && "1" && "pass"; //INSERT_TRIGGERS_HERE
  if (only_mc) base_data = baseline && "stitch && pass";
  auto proc_data = Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, names_data,base_data);

  // define combination of processes to use
  vector<shared_ptr<Process> > all_procs = vector<shared_ptr<Process> >{proc_ttx, proc_vjets, proc_other};
  if(quick_test) {
    auto proc_bkg = Process::MakeShared<Baby_pico>("All_bkg", Process::Type::background, colors("tt_1l"),
                    {base_dir+"/2016/"+mc_skim_dir+"/*TTJets_SingleLeptFromT_Tune*"}, baseline && "stitch && pass");
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
  NamedFunc nom_wgt = "w_lumi*w_isr"*HigUtilities::w_years*Higfuncs::eff_higtrig;//Higfuncs::weight_higd * Higfuncs::eff_higtrig;
  vector<string> scenarios;
  map<string, NamedFunc> weights;
  weights.emplace("no_mismeasurement", nom_wgt);
  if(alt_scen == "data"){
    scenarios = vector<string>{"data"};
  } else if(alt_scen == "bctag"){ 
    only_kappa = true;
    // example custom-defined systematic, can define multiple instead of just one, if desired
    // in general, can define any distortion of the MC with a named func and compare to default using this...
    scenarios = vector<string>();
    scenarios.push_back("syst_bctag");
    weights.emplace("syst_bctag", nom_wgt*"sys_bctag[0]");
  } else if(alt_scen == "mc_as_data"){
    scenarios = vector<string>{"mc_as_data"}; 
    weights.emplace("mc_as_data", nom_wgt);
  } else {
    only_kappa = true;
    // run on all scenarios from the sys_cfg file
    scenarios = vector<string>{alt_scen}; 
    for(const auto &scen: scenarios)
      weights.emplace(scen, nom_wgt*Functions::MismeasurementWeight(sys_wgts_file, scen));
  } 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining binning    //////////////////////////////////////////
  vector<TString> metcuts;
  string metdef = "met";
  if (sample=="zll") metdef = "ll_pt[0]";
  if (sample=="ttbar" || sample=="zll"){
    metcuts.push_back(metdef+">0&&"+metdef+"<=75");
    metcuts.push_back(metdef+">75&&"+metdef+"<=150");
  }
  metcuts.push_back(metdef+">150&&"+metdef+"<=200");
  metcuts.push_back(metdef+">200&&"+metdef+"<=300");
  metcuts.push_back(metdef+">300&&"+metdef+"<=400");
  metcuts.push_back(metdef+">400");
  if (sample=="qcd") { // add an inclusive bin
    metcuts.push_back(metdef+">150");
  } else if (sample=="ttbar" || sample=="zll"){
    metcuts.push_back(metdef+">0");
  }

  vector<TString> nbcuts;
  if (sample=="ttbar" || sample=="search" || do_highnb){
    nbcuts.push_back("nbt==2&&nbm==2");
    nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
    nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
  } else if (do_midnb){
    nbcuts.push_back("nbm==1");
    nbcuts.push_back("nbt==2&&nbm==2");
  } else {
    nbcuts.push_back("nbm==0");
    nbcuts.push_back("nbm==1");
  }

  vector<TString> planecuts;
  for (auto imet: metcuts){
    planecuts.push_back(imet + "&& hig_cand_drmax[0]<=1.1");
    planecuts.push_back(imet + "&& hig_cand_drmax[0]>1.1");
    // planecuts.push_back(imet);
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
          table_cuts.push_back(TableRow(totcut.Data(), totcut.Data(),0,0,weights.at("no_mismeasurement")));
          if(only_mc) 
            table_cuts_mm.push_back(TableRow(totcut.Data(), totcut.Data(),0,0,weights.at(iscen)));
        } // Loop over ABCD cuts
      } // Loop over bin cuts
    } // Loop over plane cuts
    TString tname = "preds"; tname += iscen;
    pm.Push<Table>(tname.Data(),  table_cuts, all_procs, true, false);
    tname += iscen;
    pm.Push<Table>(tname.Data(),  table_cuts_mm, all_procs, true, false);
  } // Loop over ABCD scenarios

  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////
  vector<TString> tablenames;
  for(size_t iscen=0; iscen<abcds.size(); iscen++) {
    // allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
    // if split_bkg [2/4] Other, [3/5] tt1l, [4/6] tt2l
    vector<vector<GammaParams> > allyields;
    Table * yield_table;
    if(only_mc){
      yield_table = static_cast<Table*>(pm.Figures()[iscen*2].get());
      Table * yield_table_mm = static_cast<Table*>(pm.Figures()[iscen*2+1].get());
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
      allyields.push_back(yield_table->Yield(proc_ttx.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_vjets.get(), lumi));
    }
    if (debug) {
      cout<<"only_mc = "<<(only_mc? "True":"False")<<endl;
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
    vector<vector<vector<float> > > kappas, kappas_mm, kmcdat, datapreds, preds;
    if (debug) cout<<"Finding predictions"<<endl;
    vector<vector<float> > yieldsPlane = findPreds(abcds[iscen], allyields, kappas, kappas_mm, kmcdat, datapreds, preds);

    if (debug) cout<<"Making tables."<<endl;
    TString fullname = printTable(abcds[iscen], allyields, kappas, datapreds, yieldsPlane, proc_sigs);
    tablenames.push_back(fullname);
    
    //// Plotting kappa
    if (debug) cout<<"Making kappa plots."<<endl;
    plotKappa(abcds[iscen], kappas, kappas_mm, kmcdat);
  } // Loop over ABCD methods

  if(!only_kappa){
    //// Printing names of ouput files
    cout<<endl<<"===== Tables that can be compiled"<<endl;
    for(size_t ind=0; ind<tablenames.size(); ind++)
      cout<<" pdflatex "<<tablenames[ind]<<"  > /dev/null"<<endl;
    cout<<endl;
  }

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
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
                   vector<vector<float> > yieldsPlane, vector<shared_ptr<Process> > &proc_sigs){
  //cout<<endl<<"Printing table (significance estimation can take a bit)"<<endl;
  //// Table general parameters
  yieldsPlane.clear();
  TString ump = " & ";



  //// Setting output file name
  TString lumi_s = RoundNumber(lumi, 0);
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
  if(print_mc) {out << "|cc"; Ncol +=2;}
  if(!only_mc || abcd.scenario.Contains("mc_as_data")) {
    out << "|cc "<<(do_zbi?"c":"");
    Ncol += 2 + (do_zbi?1:0);
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
  if(print_mc) out << "& $\\kappa$ & MC bkg.";
  if(!only_mc || abcd.scenario.Contains("mc_as_data")) out << " & Pred.& Obs. "<<(do_zbi?"& Signi.":"");
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
          out << ump <<RoundNumber(allyields[offset+2][index].Yield(), digits_table)
              << ump <<RoundNumber(allyields[offset+3][index].Yield(), digits_table)
              << ump <<RoundNumber(allyields[offset+4][index].Yield(), digits_table);
        }
        if(print_mc) {
          //// Printing kappa
          out<<ump;
          if(iabcd==3) out  << "$"    << RoundNumber(kappas[iplane][ibin][0], digits_table)
                            << "^{+"  << RoundNumber(kappas[iplane][ibin][1], digits_table)
                            << "}_{-" << RoundNumber(kappas[iplane][ibin][2], digits_table) <<"}$ ";
          //// Printing MC Bkg yields
          out << ump << RoundNumber(allyields[1][index].Yield(), digits_table);
        } // print_mc
        //// Printing background predictions
        out << ump;
        if(iabcd==3) out << "$"    << RoundNumber(preds[iplane][ibin][0], digits_table)
                         << "^{+"  << RoundNumber(preds[iplane][ibin][1], digits_table)
                         << "}_{-" << RoundNumber(preds[iplane][ibin][2], digits_table) <<"}$ ";
        //// Printing observed events in data and Obs/MC ratio
        out << ump;
        if(iabcd==3) out << RoundNumber(allyields[0][index].Yield(), 0);
        else out << RoundNumber(allyields[0][index].Yield(), 0);
        //// Printing Zbi significance
        if(do_zbi) { // "$"+RoundNumber(Significance( ),1)+"\\sigma$"
          out << ump;
          if(iabcd==3) out << "$"+RoundNumber(Significance(allyields[0][index].Yield(), preds[iplane][ibin][0], 
                                                preds[iplane][ibin][1], preds[iplane][ibin][2]),1)+"\\sigma$";
        }
        //// Printing signal yields
        if(do_signal){
          for(size_t ind=0; ind<Nsig; ind++) {
            out<<ump<<RoundNumber(allyields[2+ind][index].Yield(), digits_table);
            if(do_zbi){
              out << ump;
              if(iabcd==3) 
                out<<"$"+RoundNumber(Significance(preds[iplane][ibin][0]+allyields[2+ind][index].Yield(),
                         preds[iplane][ibin][0], preds[iplane][ibin][1], preds[iplane][ibin][2]),1)+"\\sigma$";
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
               vector<vector<vector<float> > > &kappas_mm, vector<vector<vector<float> > > &kmcdat){

  bool label_up = false; //// Putting the MET labels at the bottom
  double markerSize = 1.1;

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Kappa");
  if(label_up) opts.BottomMargin(0.11);
  if(kappas.size() >= 1) { // Used to be 4
    opts.CanvasWidth(1600);
    markerSize = 1.5;
  }
  setPlotStyle(opts);

  struct kmarker{
    TString cut;
    int color;
    int style;
    vector<float> kappa;
  };
  //// k_ordered has all the kappas grouped in sets of nb cuts (typically, in bins of njet)
  vector<vector<vector<kmarker> > > k_ordered, kmd_ordered, k_ordered_mm;
  vector<kmarker> ind_bcuts; // nb cuts actually used in the plot
  vector<float> zz; // Zero length vector for the kmarker constructor
  // vector<kmarker> bcuts({{"nbm==1",2,21,zz}, {"nbm==2",4,20,zz}, {"nbm>=3",kGreen+3,22,zz}, 
  //                                                               {"nbm==0",kMagenta+2,23,zz}, 
  //                                                               {"nbl==0",kMagenta+2,23,zz}});
  vector<kmarker> bcuts({{"none",2,21,zz}});
   
  int cSignal = kBlue;
  float maxy = 2.4, fYaxis = 1.3;
  int nbins = 0; // Total number of njet bins (used in the base histo)
  for(size_t iplane=0; iplane < kappas.size(); iplane++) {
    k_ordered.push_back(vector<vector<kmarker> >());
    kmd_ordered.push_back(vector<vector<kmarker> >());
    k_ordered_mm.push_back(vector<vector<kmarker> >());
    for(size_t ibin=0; ibin < kappas[iplane].size(); ibin++){
      if(maxy < fYaxis*(kappas[iplane][ibin][0]+kappas[iplane][ibin][1])) 
        maxy = fYaxis*(kappas[iplane][ibin][0]+kappas[iplane][ibin][1]);
      if(maxy < fYaxis*(kmcdat[iplane][ibin][0]+kmcdat[iplane][ibin][1])) 
        maxy = fYaxis*(kmcdat[iplane][ibin][0]+kmcdat[iplane][ibin][1]);
      if(maxy < fYaxis*(kappas_mm[iplane][ibin][0])) 
        maxy = fYaxis*(kappas_mm[iplane][ibin][0]);
      TString bincut = abcd.bincuts[ibin];
      bincut.ReplaceAll(" ","");
      bincut.ReplaceAll("mm_","");
      int index;
      do{
        index = bincut.First('[');
        bincut.Remove(index, bincut.First(']')-index+1);
      }while(index>=0);

      k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[0].color, bcuts[0].style, kappas[iplane][ibin]}}));
      kmd_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[0].color, bcuts[0].style, kmcdat[iplane][ibin]}}));
      k_ordered_mm[iplane].push_back(vector<kmarker>({{bincut, 1, bcuts[0].style, kappas_mm[iplane][ibin]}}));
      nbins++;
      if(ind_bcuts.size()==0) ind_bcuts.push_back(bcuts[0]);
      
    } // Loop over bin cuts
  } // Loop over plane cuts

  //// Plotting kappas
  TCanvas can("can","");
  can.SetFillStyle(4000);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);
  if(k_ordered.size()>3) label.SetTextSize(0.04);
  TLatex klab; klab.SetTextFont(42); klab.SetTextAlign(23);


  float minx = 0.5, maxx = nbins+0.5, miny = 0;
  if(label_up) maxy = 2.6;
  maxy = 3;
  if(alt_scen=="syst_mcstat") maxy = 3;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  TString ytitle = "#kappa";
  if(alt_scen!="data" && alt_scen!="syst_mcstat") ytitle += " (Scen. = "+alt_scen+")";
  histo.SetTitleOffset(0.57,"y");
  histo.SetTitleSize(0.07,"y");
  histo.SetYTitle(ytitle);
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  int bin = 0;
  vector<vector<double> > vx(ind_bcuts.size()), vexh(ind_bcuts.size()), vexl(ind_bcuts.size());
  vector<vector<double> > vy(ind_bcuts.size()), veyh(ind_bcuts.size()), veyl(ind_bcuts.size());
  vector<vector<double> > vx_mm(ind_bcuts.size()), vexh_mm(ind_bcuts.size()), vexl_mm(ind_bcuts.size());
  vector<vector<double> > vy_mm(ind_bcuts.size()), veyh_mm(ind_bcuts.size()), veyl_mm(ind_bcuts.size());
  vector<vector<double> > vx_kmd(ind_bcuts.size()), vexh_kmd(ind_bcuts.size()), vexl_kmd(ind_bcuts.size());
  vector<vector<double> > vy_kmd(ind_bcuts.size()), veyh_kmd(ind_bcuts.size()), veyl_kmd(ind_bcuts.size());
  for(size_t iplane=0; iplane < k_ordered.size(); iplane++) {
    for(size_t ibin=0; ibin < k_ordered[iplane].size(); ibin++){
      bin++;
      histo.GetXaxis()->SetBinLabel(bin, "");
      // xval is the x position of the first marker in the group
      double xval = bin, nbs = k_ordered[iplane][ibin].size(), minxb = 0.15, binw = 0;
      // If there is more than one point in the group, it starts minxb to the left of the center of the bin
      // binw is the distance between points in the njet group
      if(nbs>1) {
        xval -= minxb;
        binw = 2*minxb/(nbs-1);
      }
      for(size_t ib=0; ib<k_ordered[iplane][ibin].size(); ib++){
        //// Finding which TGraph this point goes into by comparing the color of the TGraph and the point
        for(size_t indb=0; indb<ind_bcuts.size(); indb++){
          if(ind_bcuts[indb].color == k_ordered[iplane][ibin][ib].color){
            vx[indb].push_back(xval);
            vexl[indb].push_back(0);
            vexh[indb].push_back(0);
            vy[indb].push_back(k_ordered[iplane][ibin][ib].kappa[0]);
            veyh[indb].push_back(k_ordered[iplane][ibin][ib].kappa[1]);
            veyl[indb].push_back(k_ordered[iplane][ibin][ib].kappa[2]);

            //// MC kappas with data uncertainties
            vx_kmd[indb].push_back(xval);
            vexl_kmd[indb].push_back(0);
            vexh_kmd[indb].push_back(0);
            vy_kmd[indb].push_back(kmd_ordered[iplane][ibin][ib].kappa[0]);
            float ekmdUp = sqrt(pow(k_ordered[iplane][ibin][ib].kappa[1],2) +
                                pow(kmd_ordered[iplane][ibin][ib].kappa[1],2));
            float ekmdDown = sqrt(pow(k_ordered[iplane][ibin][ib].kappa[2],2) +
                                  pow(kmd_ordered[iplane][ibin][ib].kappa[2],2));
            veyh_kmd[indb].push_back(ekmdUp);            
            veyl_kmd[indb].push_back(ekmdDown);

            //// Data/pseudodata kappas
            vx_mm[indb].push_back(xval+0.1);
            vexl_mm[indb].push_back(0);
            vexh_mm[indb].push_back(0);
            vy_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[0]);
            if(alt_scen=="data" || alt_scen=="mc_as_data") {
              veyh_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[1]);
              veyl_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[2]);         
            } else {     
              veyh_mm[indb].push_back(0);
              veyl_mm[indb].push_back(0);
            }

            // Printing difference between kappa and kappa_mm
            float kap = k_ordered[iplane][ibin][ib].kappa[0], kap_mm = k_ordered_mm[iplane][ibin][ib].kappa[0];
            TString text = "#Delta_{#kappa} = "+RoundNumber((kap_mm-kap)*100,0,kap)+"%";
            if(alt_scen=="data") text = "#Delta_{#kappa} = "+RoundNumber((kap_mm-1)*100,0,1)+"%";
            else text = "#Delta_{#kappa} = "+RoundNumber((kap-1)*100,0,1)+"%";
            klab.SetTextColor(1);
            klab.SetTextSize(0.045);
            if (k_ordered.size()*k_ordered[iplane].size()>8) klab.SetTextSize(0.03);
            if(alt_scen!="mc_as_data") klab.DrawLatex(xval, 0.952*maxy, text);
            //// Printing stat uncertainty of kappa_mm/kappa
            float kapUp = k_ordered[iplane][ibin][ib].kappa[1], kapDown = k_ordered[iplane][ibin][ib].kappa[2];
            float kap_mmUp = k_ordered_mm[iplane][ibin][ib].kappa[1];
            float kap_mmDown = k_ordered_mm[iplane][ibin][ib].kappa[2];
            if(alt_scen=="data" || alt_scen=="mc_as_data") {              
              text = "#sigma_{stat} = ^{+"+RoundNumber(kap_mmUp*100,0, 1)+"%}_{-"+RoundNumber(kap_mmDown*100,0, 1)+"%}";
            } else {
              text = "#sigma_{stat} = ^{+"+RoundNumber(kapUp*100,0, 1)+"%}_{-"+RoundNumber(kapDown*100,0, 1)+"%}";
            }
            // adding label to indicate the ABCD corresponding to each kappa value
            klab.SetTextSize(0.05);
            if (k_ordered.size()*k_ordered[iplane].size()>8) klab.SetTextSize(0.035);
            if (sample=="search" || sample=="ttbar") text = ibin%2==0 ? "3b/2b" : "4b/2b"; 
            else if (sample=="zll") text = do_midnb ? "2b/1b" : "1b/0b"; 
            else if (sample=="qcd") {
              if (do_highnb) text =  ibin%2==0 ? "3b/2b" : "4b/2b"; 
              else text = do_midnb ? "2b/1b" : "1b/0b"; 
            }
            klab.DrawLatex(xval, 0.92*maxy, text);
            xval += binw;
          }
        } 
      } 
    } // Loop over bin cuts

    // Drawing lines separating MET planes
    line.SetLineStyle(iplane%2==0 ? 2:1); line.SetLineWidth(1); line.SetLineColor(kBlack);
    if (iplane<k_ordered.size()-1) line.DrawLine(bin+0.5, miny, bin+0.5, maxy);

    // All the labels on the X-axis...
    string metdef = "met";
    if (sample=="zll") metdef = "ll_pt";
    double lmargin(opts.LeftMargin()), rmargin(opts.RightMargin()), bmargin(opts.BottomMargin());
    label.SetTextSize(0.05);
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
        label.SetTextAlign(23);
        label.DrawLatexNDC(lmargin+(1-rmargin-lmargin)/abcd.planecuts.size()*(iplane+1), 0.13, metlabel.c_str());
      }
    } else { // if not binning in dRmax
      string plabel = CodeToRootTex(abcd.planecuts[iplane].Data());
      ReplaceAll(plabel, "<"+metdef+"#leq","#minus"); 
      if(plabel == "> 0") plabel = "Inclusive";
      label.SetTextSize(0.05);
      label.SetTextAlign(23);
      label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., -0.03*maxy, plabel.c_str());
    }
    // write pT miss at the bottom of the plot...
    if (iplane==0) {
      if (sample=="zll") label.DrawLatexNDC(lmargin+(1-rmargin-lmargin)/2., bmargin-0.08, 
                                          "p_{T}(l#kern[0.11]{#lower[-0.25]{^{+}}}#kern[0.3]{l#kern[0.05]{#lower[0.36]{^{#scale[1.3]{-}}}}}) [GeV]");
      else label.DrawLatexNDC(lmargin+(1-rmargin-lmargin)/2., bmargin-0.08, "p_{T}^{miss} [GeV]");
    }
  } // Loop over plane cuts

  //// Drawing legend and TGraphs
  if (debug) cout<<"Building up TGraphs"<<endl;
  double legX(opts.LeftMargin()+0.005), legY(1-0.035), legSingle = 0.05;
  legX = 0.595;
  if(only_mc) legX = 0.65;
  if(label_up) legY = 0.8;
  double legW = 0.15, legH = legSingle*(ind_bcuts.size()+1)/2;
  if(ind_bcuts.size()>3) legH = legSingle*((ind_bcuts.size()+1)/2);
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()*1.15); leg.SetFillColor(0);leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(2);
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  TGraphAsymmErrors graph_kmd[20]; 
  TGraphAsymmErrors graph_mm[20]; 
  for(size_t indb=0; indb<ind_bcuts.size(); indb++){
    graph_kmd[indb] = TGraphAsymmErrors(vx_kmd[indb].size(), &(vx_kmd[indb][0]), &(vy_kmd[indb][0]),
                                        &(vexl_kmd[indb][0]), &(vexh_kmd[indb][0]), &(veyl_kmd[indb][0]), 
                                        &(veyh_kmd[indb][0]));
    graph_kmd[indb].SetMarkerStyle(ind_bcuts[indb].style); graph_kmd[indb].SetMarkerSize(markerSize);
    graph_kmd[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph_kmd[indb].SetLineColor(1); graph_kmd[indb].SetLineWidth(2);
    if(alt_scen=="mc_as_data") graph_kmd[indb].Draw("p0 same");

    graph[indb] = TGraphAsymmErrors(vx[indb].size(), &(vx[indb][0]), &(vy[indb][0]),
                                    &(vexl[indb][0]), &(vexh[indb][0]), &(veyl[indb][0]), &(veyh[indb][0]));
    graph[indb].SetMarkerStyle(ind_bcuts[indb].style); graph[indb].SetMarkerSize(markerSize);
    graph[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph[indb].SetLineColor(ind_bcuts[indb].color); graph[indb].SetLineWidth(2);
    graph[indb].Draw("p0 same");

    graph_mm[indb] = TGraphAsymmErrors(vx_mm[indb].size(), &(vx_mm[indb][0]), &(vy_mm[indb][0]),
                                       &(vexl_mm[indb][0]), &(vexh_mm[indb][0]), &(veyl_mm[indb][0]), 
                                       &(veyh_mm[indb][0]));
    graph_mm[indb].SetMarkerStyle(20); graph_mm[indb].SetMarkerSize(markerSize*1.2);
    graph_mm[indb].SetMarkerColor(1);
    graph_mm[indb].SetLineColor(1); graph_mm[indb].SetLineWidth(2);
    if(alt_scen!="mc_as_data" && alt_scen!="mc_only") graph_mm[indb].Draw("p0 same");

    leg.AddEntry(&graph[indb], "MC", "ep");
    TString data_s = (alt_scen=="data"||alt_scen=="no_mismeasurement"?"Data":"Pseudodata");
    if(alt_scen!="mc_as_data") 
      leg.AddEntry(&graph_mm[indb], data_s, "ep");

  } // Loop over TGraphs
  leg.Draw();

  //// Drawing CMS labels and line at 1
  TString cmsPrel = "#font[62]{CMS} #scale[0.8]{#font[52]{}}";
  TString cmsSim = "#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}";
  TLatex cmslabel;
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  if(only_mc) cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,cmsSim);
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
  if(alt_scen!="mc_as_data" && alt_scen!="syst_mcstat") title = "#font[42]{"+RoundNumber(lumi, 0)+" fb^{-1} (13 TeV)}";
  if(abcd_title.Contains("Search") && only_mc) title = "#font["+fontstyle+"]{"+abcd_title+"}";
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, title);

  ///// Sample name
  cmslabel.SetTextAlign(11);
  title = "#font["+fontstyle+"]{"+abcd_title+"}";
  TString newSignal = "#color["; newSignal += cSignal; newSignal += "]{Signal}";
  title.ReplaceAll("Signal", newSignal);
  if(!(abcd_title.Contains("Search") && only_mc)) cmslabel.DrawLatex(opts.LeftMargin()+0.14, 1-opts.TopMargin()+0.015, title);

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname="plots/kappa_"+sample+"_tight_" +abcd.scenario + (do_highnb?"_highnb":"")+(do_midnb?"_highnb":"");
  fname += "_lumi"+RoundNumber(lumi, 0)+".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl; 

}

//// Calculating kappa and Total bkg prediction
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
vector<vector<float> > findPreds(abcd_def &abcd, vector<vector<GammaParams> > &allyields,
               vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &kappas_mm, 
               vector<vector<vector<float> > > &kmcdat, vector<vector<vector<float> > > &datapreds,
               vector<vector<vector<float> > > &preds){
  // Powers for kappa:   ({R1, R2, D3, R4})
  vector<float> pow_kappa({ 1, -1, -1,  1});
  // Powers for TotBkg pred:({R1, R2, D3,  R1, R2, D3, D4})
  vector<float> pow_totpred( {-1,  1,  1,   1, -1, -1,  1});
  vector<float> pow_datapred( {-1,  1,  1});

  float val(1.), valup(1.), valdown(1.);
  vector<vector<float> > yieldsPlane;
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    //// Counting yields in plane without double-counting R1/R3 yields when integrated
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
    float dataMC = Nobs/Nmc;
    float edataMC = sqrt(pow(sqrt(Nobs)/Nmc,2) + pow(Nobs*NmcPlane.Uncertainty()/Nmc/Nmc,2));
    yieldsPlane.push_back({dataMC, edataMC});
    

    kappas.push_back(vector<vector<float> >());
    kmcdat.push_back(vector<vector<float> >());
    kappas_mm.push_back(vector<vector<float> >());
    preds.push_back(vector<vector<float> >());
    datapreds.push_back(vector<vector<float> >());
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
      datapreds[iplane].push_back(vector<float>({val, valup, valdown}));

      vector<vector<float> > kentries;
      vector<vector<float> > kweights;
      vector<vector<float> > kkentries;
      vector<vector<float> > kkweights;
      vector<vector<float> > kentries_mm;
      vector<vector<float> > kweights_mm;
      //// Pushing MC yields for predictions and kappas
      for(size_t iabcd=0; iabcd < 4; iabcd++){
        size_t index = iplane*abcd.bincuts.size()*4 +ibin*4 + iabcd;
        // Renormalizing MC to data
        // allyields[1][index] *= dataMC;

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

      // Throwing toys to find predictions and uncertainties
      val = calcKappa(entries, weights, pow_totpred, valdown, valup);
      if(valdown<0) valdown = 0;
      preds[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries, kweights, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kappas[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries_mm, kweights_mm, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kappas_mm[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kkentries, kkweights, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kmcdat[iplane].push_back(vector<float>({val, valup, valdown}));
    } // Loop over bin cuts
  } // Loop over plane cuts

  return yieldsPlane;
} // findPreds

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 's'},    // Which sample to use: standard, 2015 data
      {"only_kappa", no_argument, 0, 'k'},    // Only plots kappa (no table)
      {"debug", no_argument, 0, 'd'},         // Debug: prints yields and cuts used
      {"quick", no_argument, 0, 'q'},           // Used inclusive ttbar for quick testing
      {"scen", required_argument, 0, 0},        // Mismeasurment scenario, 0 for data
      {"midnb", no_argument, 0, 0},           // Check zll and qcd CRs for 2b
      {"highnb", no_argument, 0, 0},          // Check QCD CR at 3b and 4b
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:kdq", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      sample = optarg;
      break;
    case 'k':
      only_kappa = true;
      break;
    case 'd':
      debug = true;
      break;
    case 'q':
      quick_test = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "scen"){
        alt_scen = optarg;
      }else if(optname == "highnb"){
        do_highnb = true;
      }else if(optname == "midnb"){
        do_midnb = true;
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
