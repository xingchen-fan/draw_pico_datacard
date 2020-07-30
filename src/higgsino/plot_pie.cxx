///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/slide_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  //fixme:simplify options
  float lumi = 1;
  // options "zll", "qcd", "ttbar", "search"
  string sample_name = "search";
  bool resolved = true;
}
  
int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  PlotOpt lin_shapes_info("txt/plot_styles.txt", "CMSPaper");
  lin_shapes_info.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::shapes);
  PlotOpt log_shapes_info = lin_shapes_info;
  log_shapes_info.YAxis(YAxisType::log);
  vector<PlotOpt> plt_types = {lin_shapes_info, log_shapes_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  set<int> years; 
  years = {2016, 2017, 2018};

  string base_dir(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/");
  string mc_skim_dir("mc/merged_higmc_preselect/");
  if (sample_name=="ttbar")    mc_skim_dir = "mc/merged_higmc_higlep1T/"; 
  else if (sample_name=="zll") mc_skim_dir = "mc/merged_higmc_higlep2T/"; 
  else if (sample_name=="qcd") mc_skim_dir = "mc/merged_higmc_higqcd/";  

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                                  "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"
                                  });
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({
                                   // "*QCD_HT100to200_Tune*", 
                                   "*QCD_HT200to300_Tune*", 
                                   "*QCD_HT300to500_Tune*",  // these have very low stats
                                   "*QCD_HT500to700_Tune*",
                                 "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                 "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                   "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  set<string> allmctags;
  for (auto &iset: mctags) {
      allmctags.insert(iset.second.begin(), iset.second.end());
  }

  // Baseline definitions
  NamedFunc wgt = "w_lumi*w_isr"*HigUtilities::w_years*Higfuncs::eff_higtrig;//Higfuncs::weight_higd * Higfuncs::eff_higtrig;

  TString c_hig_trim = "hig_cand_drmax[0]<=2.2  && hig_cand_dm[0]<=40 && hig_cand_am[0]<200";
  string baseline = "njet>=4 && njet<=5 && nbt>=2";
  if (sample_name=="search") baseline += " && nvlep==0 && ntk==0 && !low_dphi_met &&"+c_hig_trim;
  else if (sample_name=="ttbar") baseline += " && nlep==1 && mt<100 &&"+c_hig_trim;
  else if (sample_name=="zll") baseline += " && nlep==2 && met<50 &&"+c_hig_trim;
  else if (sample_name=="qcd") baseline += " && nvlep==0 && ntk==0 && low_dphi_met &&"+c_hig_trim;

  NamedFunc base_func = baseline && Functions::hem_veto && "stitch && pass";// && met/mht<2 && met/met_calo<2";

  map<string, vector<shared_ptr<Process> >> proc_sets;
  proc_sets["base"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background, kBlue-7,
    attach_folder(base_dir, years, mc_skim_dir, mctags["ttx"]),     base_func && "ntrutauh>0"));
  proc_sets["base"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background, colors("tt_1l"),
    attach_folder(base_dir, years, mc_skim_dir, mctags["ttx"]),     base_func && "ntrutauh==0"));
  proc_sets["base"].push_back(Process::MakeShared<Baby_pico>("V+jets",     Process::Type::background, kOrange+1,
    attach_folder(base_dir, years, mc_skim_dir, mctags["vjets"]),   base_func));
  proc_sets["base"].push_back(Process::MakeShared<Baby_pico>("Single t",   Process::Type::background, colors("single_t"),
    attach_folder(base_dir, years, mc_skim_dir, mctags["singlet"]), base_func));
  proc_sets["base"].push_back(Process::MakeShared<Baby_pico>("QCD",        Process::Type::background, colors("other"),
    attach_folder(base_dir, years, mc_skim_dir, mctags["qcd"]),     base_func)); 
  proc_sets["base"].push_back(Process::MakeShared<Baby_pico>("Other",      Process::Type::background, kGreen+1,
    attach_folder(base_dir, years, mc_skim_dir, mctags["other"]),   base_func));      

  set<string> allfiles = attach_folder(base_dir, years, mc_skim_dir, allmctags);
  proc_sets["ntrub"].push_back(Process::MakeShared<Baby_pico>
    ("0 B-hadron",       Process::Type::background, kAzure-4,   allfiles, base_func && Functions::ntrub<1));
  proc_sets["ntrub"].push_back(Process::MakeShared<Baby_pico>
    ("1 B-hadron",       Process::Type::background, kTeal-8,    allfiles, base_func && Functions::ntrub==1));
  proc_sets["ntrub"].push_back(Process::MakeShared<Baby_pico>
    ("2 B-hadrons",      Process::Type::background, kOrange-4,  allfiles, base_func && Functions::ntrub==2));
  proc_sets["ntrub"].push_back(Process::MakeShared<Baby_pico>
    ("3 B-hadrons",      Process::Type::background, kPink+2,    allfiles, base_func && Functions::ntrub==3));
  proc_sets["ntrub"].push_back(Process::MakeShared<Baby_pico>
    ("4 B-hadrons", Process::Type::background, kMagenta-1, allfiles, base_func && Functions::ntrub==4));

  vector<string> xcuts; // useful additional cut definitions
  xcuts.push_back("hig_cand_drmax[0]<=1.1");
  xcuts.push_back("hig_cand_drmax[0]>1.1");
  
  vector<string> nbcuts;
  if (sample_name=="zll" || sample_name=="qcd") {
    nbcuts.push_back("nbm==0");
    nbcuts.push_back("nbm==1");
    nbcuts.push_back("nbm==2");
  }
  if (sample_name=="ttbar" || sample_name=="search") {
    nbcuts.push_back("nbt==2&&nbm==2");
    nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
    nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
  }

  vector<string> metcuts;
  string metdef = "met";
  if (sample_name=="zll") metdef = "ll_pt[0]";
  if (sample_name=="search" || sample_name=="qcd"){
    metcuts.push_back(metdef+">150&&"+metdef+"<=200");
    metcuts.push_back(metdef+">200&&"+metdef+"<=300");
    metcuts.push_back(metdef+">300&&"+metdef+"<=400");
    metcuts.push_back(metdef+">400");
  } else if (sample_name=="ttbar" || sample_name=="zll") {
    metcuts.push_back(metdef+">0&&"+metdef+"<=75");
    metcuts.push_back(metdef+">75&&"+metdef+"<=150");
    metcuts.push_back(metdef+">150&&"+metdef+"<=200");
    metcuts.push_back(metdef+">200&&"+metdef+"<=300");
    metcuts.push_back(metdef+">300");
  }

  PlotMaker pm;
  SlideMaker sm("slide_pies_"+sample_name+".tex","1610");
  vector<string> pnames;
  vector<string> row_labels; vector<string> col_labels; 

  //  Make pie table with Nb as rows and additional cuts as columns
  //-----------------------------------------------------------------
  for(auto &ixcut: xcuts) 
    row_labels.push_back("$"+CodeToLatex(ixcut)+"$");
  for(auto &inb: nbcuts) 
    col_labels.push_back("$"+CodeToLatex(inb)+"$");

  vector<TableRow> table_cuts;
  for(auto &proc_set: proc_sets) {
    string tabname = sample_name+"_"+proc_set.first;
    pnames.clear();
    for(auto &ixcut: xcuts) {
      for(auto &inb: nbcuts) {
        string icut = inb+"&&"+ixcut;
        table_cuts.push_back(TableRow("", icut, 0, 0, wgt));  
        pnames.push_back("pie_"+tabname+"_"+CodeToPlainText(icut)+"_perc_lumi"+RoundNumber(lumi,0).Data()+".pdf");
      }
    }
    pm.Push<Table>(tabname, table_cuts, proc_set.second, true, true, true);
    sm.AddSlide(pnames, nbcuts.size(), "X-axis: Number of b-tags, Y-axis: additional cuts", col_labels, row_labels);  
  }
  
  //  Make pie table with Nb as rows and additional cuts as columns,
  //  one slide per set of additional cuts
  //--------------------------------------------------------------------------------
  table_cuts.clear(); row_labels.clear(); col_labels.clear(); 
  for(auto &imet: metcuts) 
    col_labels.push_back("$"+CodeToLatex(imet)+"$");
  for(auto &inb: nbcuts) 
    row_labels.push_back("$"+CodeToLatex(inb)+"$");

  for(auto &proc_set: proc_sets) {
    string tabname = sample_name+"_"+proc_set.first;
    for(auto &ixcut: xcuts) {
      pnames.clear();
      string slide_ttl = "Additional cuts: $"+CodeToLatex(ixcut)+"$";
      if (ixcut=="1") slide_ttl = "No additional cuts";
      for(auto &inb: nbcuts) {
        for(auto &imet: metcuts) {
          string icut = inb+"&&"+imet+"&&"+ixcut;
          table_cuts.push_back(TableRow("", icut, 0, 0, wgt));  
          pnames.push_back("pie_"+tabname+"_"+CodeToPlainText(icut)+"_perc_lumi"+RoundNumber(lumi,0).Data()+".pdf");
        }
      }
      sm.AddSlide(pnames, metcuts.size(), slide_ttl, col_labels, row_labels);
    }
    pm.Push<Table>(tabname,  table_cuts, proc_set.second, true, true, true);
  }

  pm.min_print_ = true;
  pm.multithreaded_ = true;

  pm.MakePlots(lumi);
  sm.Close();

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"bcats", required_argument, 0, 'b'},
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"sample", required_argument, 0, 's'},    // Which sample_name to use: standard, met150, 2015 data
      {"boo", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:l:gn:d:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 's':
      sample_name = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "boo"){
        resolved = false;
      } else {
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
