///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>
#include <regex>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "TLegend.h" // Controls error level reporting
#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/ordered_dict.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  float lumi = 1;
  bool doSignal = true;
  string tag = "";
}

float to_float(string const & inString) {
    std::smatch match;
    std::regex pattern("-?[0-9]+(\\.[0-9]*)?");
    std::regex_search (inString, match, pattern);
    //cout<<"to_float: "<<inString<<" :: "<<match.str()<<endl;
    return stof(match.str());
}

void regSearchAll(string inString, string const & regExPattern, vector<string> & stringList, vector<int> & positionList) {
    std::smatch match;
    std::regex pattern(regExPattern);
    //cout<<inString<<endl;
    int position_cursor = 0;
    while (std::regex_search (inString, match, pattern)) {
      //cout<<match.str()<<" "<<position_cursor + match.position()<<endl;
      stringList.push_back(match.str());
      positionList.push_back(position_cursor + match.position());
      inString = inString.substr(match.position()+1);
      position_cursor += match.position()+1;
    }
}

string regReplace(string const & inString, string const & regExPattern, string const & replace) {
  std::regex pattern(regExPattern);
  return std::regex_replace(inString, pattern, replace);
}

void getValuesFromTexTable(string const & filename, vector<vector<float> > & tableValues, vector<vector<float> > & tableUncertainties) {
  string line;
  ifstream texTable (filename);
  // tableValues[iRow][iColumn] = value
  while (getline(texTable, line)) {
    //cout<<line<<endl;
    vector<string> stringList;
    vector<int> positionList;
    regSearchAll(line, "(\\&[ ]*-?[0-9]+\\.[0-9]+[ ]*(\\&|\\\\\\\\))|(\\&[ ]*-?[0-9]+\\.[0-9]+\\$\\\\pm\\$[0-9]+\\.[0-9]+[ ]*(\\&|\\\\\\\\))", stringList, positionList);
    //cout<<"StringList size: "<<stringList.size()<<endl;
    //for (auto const & valueString : stringList) cout<<valueString<<endl;;
    vector<float> rowValues;
    vector<float> rowUncertainties;
    // Split into values and uncertainties
    for (auto const & valueString : stringList) {
      //cout<<valueString<<endl;
      if (Contains(valueString,"pm")) {
        vector<string> splitString;
        HigUtilities::stringToVectorString(valueString, splitString, "$\\pm$");
        rowValues.push_back(to_float(splitString[0]));
        rowUncertainties.push_back(to_float(splitString[1]));
        //for(auto const & item : splitString) cout<<item<<endl;
      } else {
        rowValues.push_back(to_float(valueString));
        rowUncertainties.push_back(0);
      }
    }
    tableValues.size(); tableUncertainties.size();
    //for (auto const & valueString : stringList) rowValues.push_back(to_float(valueString));
    if(rowValues.size()!=0) tableValues.push_back(rowValues);
    if(rowUncertainties.size()!=0) tableUncertainties.push_back(rowUncertainties);
  }
  texTable.close();
}

void replaceValuesInTexTablename(string const & filename, string const & newFilename, vector<vector<float> > const & tableValues, vector<vector<float> > const & tableUncertainties, int digits) {
  string newTexTable;
  string line;
  ifstream texTable (filename);
  int iRow = 0;
  while (getline(texTable, line)) {
    //cout<<"[org] "<<line<<endl;
    vector<string> stringList;
    vector<int> positionList;
    regSearchAll(line, "(\\&[ ]*-?[0-9]+\\.[0-9]+[ ]*(\\&|\\\\\\\\))|(\\&[ ]*-?[0-9]+\\.[0-9]+\\$\\\\pm\\$[0-9]+\\.[0-9]+[ ]*(\\&|\\\\\\\\))", stringList, positionList);
    if (stringList.size()!=0) {
      if (tableValues[iRow].size() != stringList.size()) cout<<"[Error] replaceValuesInTexTablename:: tableValue[iRow] size:"+to_string(tableValues[iRow].size())+" and stringList size:"+to_string(stringList.size())+" do not match"<<endl;
      int position_cursor = 0;
      for (unsigned iColumn=0; iColumn<tableValues[iRow].size(); ++iColumn) {
        //cout<<"[org] "<<stringList[iColumn]<<" tableValue: "<<tableValues[iRow][iColumn]<<" tableUncertainties: "<<tableUncertainties[iRow][iColumn]<<" "<<RoundNumber(tableUncertainties[iRow][iColumn], digits, 1)<<endl;
        string newColumn;
        if (Contains(stringList[iColumn], "pm")) {
          // $$ is to escape from $#
          newColumn = regReplace(stringList[iColumn],"-?[0-9]+(\\.[0-9]*)?\\$", string(RoundNumber(tableValues[iRow][iColumn], digits, 1).Data())+"$$");
          newColumn = regReplace(newColumn,"\\$[0-9]+\\.[0-9]+", "$$"+string(RoundNumber(tableUncertainties[iRow][iColumn], digits, 1).Data()));
        } else {
          newColumn = regReplace(stringList[iColumn],"-?[0-9]+(\\.[0-9]*)?", RoundNumber(tableValues[iRow][iColumn], digits, 1).Data());
        }
        //cout<<"[new] "<<newColumn<<endl;
        // Need to replace line
        line.replace(positionList[iColumn] + position_cursor, stringList[iColumn].size(), newColumn);
        position_cursor += newColumn.size() - stringList[iColumn].size();
        //cout<<"[new] "<<line<<endl;
      }
      iRow++;
    }
    newTexTable += line + "\n";
  }
  texTable.close();
 
  ofstream outFile(newFilename);
  outFile<<newTexTable;
  outFile.close();
  //cout<<"[Info] replaceValuesInTexTablename:: Wrote new tex table: "<<newFilename<<endl;
  cout<<"pdflatex "<<newFilename<<" &> /dev/null; #./python/texify.py tables"<<endl;
}

const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double weight = 1;
  //if (b.type()==106000) {
  //  return 35.9;
  //}
  if (b.SampleType()==2016){
    return weight*35.9;
  } else if (b.SampleType()==2017){
    return weight*41.5;
  } else {
    return weight*59.6;
  }
});

const NamedFunc w_years_scaleup("w_years_scaleup", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double weight = 1;
  //if (b.type()==106000) {
  //  return 35.9;
  //}
  return weight*137;
});

string getLuminosityString(string const & year_string) {
  set<int> years;
  if (year_string == "run2") years = {2016, 2017, 2018};
  else HigUtilities::parseYears(year_string, years);
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();
  return total_luminosity_string;
}

// nanoAodFolder: /net/cms29/cms29r0/pico/NanoAODv5/
// production: higgsino_humboldt
// dataType: mc/data/signal
// fileTag: (tt, single_t, zjets, wjets, qcd, other, all), (all), (200_1, 600_1, 950_1, ...)
// sample: search/ttbar/zll/qcd
// year_string: 2016/2017/2018/run2
void addProcess(string const & processName, Process::Type type, int color, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & dataType, string const & fileTag, string const & sample, string const & year_string, 
  vector<shared_ptr<Process> > & procs, bool use_met150_skim) {
  set<int> years;
  //years_a = {2016};
  if (year_string == "run2") years = {2016, 2017, 2018};
  else HigUtilities::parseYears(year_string, years);

  // Set base tags
  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
  //mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  // Combine all tags
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  set<string> fileNames;
  if (dataType == "data") fileNames = {"*.root"};
  else if (dataType == "signal") {
    vector<string> massPoints;
    HigUtilities::stringToVectorString(fileTag, massPoints, "_");
    // Special case for angeles production
    if (production == "higgsino_angeles" && massPoints[1] == "0") massPoints[1] = "1";
    fileNames = {"*TChiHH_mChi-"+massPoints[0]+"_mLSP-"+massPoints[1]+"_*.root"};
  }
  else fileNames = mctags[fileTag];

  // Set folders
  torch::OrderedDict<string, string> folderDict;
  string mcProductionFolder = nanoAodFolder+"/"+production;
  string dataProductionFolder = nanoAodFolder+"/"+production;
  string signalProductionFolder = nanoAodFolder+"/"+production;
  folderDict.insert("mc_production_folder", mcProductionFolder);
  // preselect:
  //   ((nbt>=2 && njet>=4 && njet<=5)||(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1))
  //   nvlep==0 && ntk==0 && !low_dphi_met && met>150 && 
  //folderDict.insert("search_mc_skim_folder", "mc/merged_higmc_higloose/");
  // higloose: 
  //   (nbt>=2 || nbdft>=2 || Sum$(fjet_pt>300 && fjet_msoftdrop>50)>0)&&
  //   met>150 && nvlep==0
  if (use_met150_skim) folderDict.insert("search_mc_skim_folder", "mc/merged_higmc_met150/");
  else folderDict.insert("search_mc_skim_folder", "mc/merged_higmc_higloose/");
  // higlep1T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || ((nbt>=2 || nbdft>=2) && njet>=4 && njet<=5)) &&
  //   nlep==1 && 
  //   (Max$(el_pt*el_sig)>40 || Max$(mu_pt*mu_sig)>40) // pass_1l_trig40 (sig is signal lepton)
  folderDict.insert("ttbar_mc_skim_folder", "mc/merged_higmc_higlep1T/");
  // higlep2T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nlep==2 && 
  //   @ll_m.size()>=1 && Sum$(ll_m>80 && ll_m<100)>=1
  //   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) // pass_2l_trig30
  folderDict.insert("zll_mc_skim_folder", "mc/merged_higmc_higlep2T/");
  // higqcd with met150:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nvlep==0 && ntk==0 && low_dphi_met &&
  //   met>150  // Since applied to met150 skim
  folderDict.insert("qcd_mc_skim_folder", "mc/merged_higmc_higqcd/");

  folderDict.insert("data_production_folder", dataProductionFolder);
  if (use_met150_skim) folderDict.insert("search_data_skim_folder","data/merged_higmc_met150/");
  else folderDict.insert("search_data_skim_folder","data/merged_higmc_higloose/");
  folderDict.insert("ttbar_data_skim_folder","data/merged_higmc_higlep1T/");
  folderDict.insert("zll_data_skim_folder","data/merged_higmc_higlep2T/");
  folderDict.insert("qcd_data_skim_folder","data/merged_higmc_higqcd/");

  folderDict.insert("signal_production_folder", signalProductionFolder);
  //folderDict.insert("signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higloose/");
  //folderDict.insert("signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_preselect/");
  if (use_met150_skim) folderDict.insert("search_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_met150/");
  else folderDict.insert("search_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higloose/");
  folderDict.insert("ttbar_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higlep1T/");
  folderDict.insert("zll_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higlep2T/");
  folderDict.insert("qcd_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higqcd/");

  // Set paths
  set<string> pathNames;
  pathNames = attach_folder(folderDict[dataType+"_production_folder"], years, folderDict[sample+"_"+dataType+"_skim_folder"], fileNames);

  procs.push_back(Process::MakeShared<Baby_pico>(processName, type, color, pathNames,"stitch"&&additionalCut));
}

void addAllMcProcesses(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & sample, string const & year_string, 
  vector<shared_ptr<Process> > & procs, bool use_met150_skim) {
  Palette colors("txt/colors.txt", "default");
  addProcess("t#bar{t}+X (#tau_{had}>0) "+processName_postfix, Process::Type::background, colors("tt_htau"), "ntrutauh>0"&&additionalCut,
    nanoAodFolder, production, "mc", "tt", sample, year_string,
    procs, use_met150_skim);
  addProcess("t#bar{t}+X (#tau_{had}=0) "+processName_postfix, Process::Type::background, colors("tt_1l"), "ntrutauh==0"&&additionalCut,
    nanoAodFolder, production, "mc", "tt", sample, year_string,
    procs, use_met150_skim);
  addProcess("Z+Jets "+processName_postfix, Process::Type::background, kOrange+1, additionalCut,
    nanoAodFolder, production, "mc", "zjets", sample, year_string,
    procs, use_met150_skim);
  addProcess("W+Jets "+processName_postfix, Process::Type::background, kGreen+1, additionalCut,
    nanoAodFolder, production, "mc", "wjets", sample, year_string,
    procs, use_met150_skim);
  addProcess("Single t "+processName_postfix, Process::Type::background, colors("single_t"), additionalCut,
    nanoAodFolder, production, "mc", "single_t", sample, year_string,
    procs, use_met150_skim);
  addProcess("QCD "+processName_postfix, Process::Type::background, colors("other"), additionalCut,
    nanoAodFolder, production, "mc", "qcd", sample, year_string,
    procs, use_met150_skim);
  addProcess("Other "+processName_postfix, Process::Type::background, kGray+2, additionalCut,
    nanoAodFolder, production, "mc", "other", sample, year_string,
    procs, use_met150_skim);
}

void addMultipleSignalProcesses(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & sample, string const & year_string, 
  vector<shared_ptr<Process> > & procs, bool use_met150_skim) 
{
  vector<string> sigm = {"175", "550", "950"};
  for (unsigned isig(0); isig<sigm.size(); isig++){
    addProcess("TChiHH("+sigm[isig]+",0) "+processName_postfix, Process::Type::signal, 1, additionalCut,
      nanoAodFolder, production, "signal", sigm[isig]+"_0", sample, year_string,
      procs, use_met150_skim);
  }
}

void constructRegionTable(vector<TableRow> & tablerows, 
  string const & sample, NamedFunc weight,
  bool const & bin_met, bool const &bin_drmax) {
  vector<TString> vc_drmax;
  vc_drmax.push_back("hig_cand_drmax[0]<=1.1");
  vc_drmax.push_back("hig_cand_drmax[0]>1.1");
  if (!bin_drmax) vc_drmax = {"hig_cand_drmax[0]>0"};
  vector<string> res_abcd;
  vector<TString> vc_met;

  string c_2b, c_3b, c_4b, hig, sbd;
  if (sample=="qcd") {
    c_2b = "nbm==0";
    c_3b = "nbm==1";
    hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
    res_abcd = {c_2b+"&&"+sbd, c_2b+"&&"+hig, c_3b+"&&"+sbd, c_3b+"&&"+hig};

    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300 && met<=400");
    vc_met.push_back("met>400");
    if (!bin_met) vc_met = {"met>150"};
  } else if (sample=="zll") {
    c_2b = "nbm==0";
    c_3b = "nbm==1";
    hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
    res_abcd = {c_2b+"&&"+sbd, c_2b+"&&"+hig, c_3b+"&&"+sbd, c_3b+"&&"+hig};

    vc_met.push_back("ll_pt[0]>0   && ll_pt[0]<=75");
    vc_met.push_back("ll_pt[0]>75  && ll_pt[0]<=150");
    vc_met.push_back("ll_pt[0]>150 && ll_pt[0]<=200");
    vc_met.push_back("ll_pt[0]>200 && ll_pt[0]<=300");
    vc_met.push_back("ll_pt[0]>300");
    if (!bin_met) vc_met = {"ll_pt[0]>0"};
  } else if (sample=="ttbar") {
    c_2b = "(nbt==2&&nbm==2)";
    c_3b = "(nbt>=2&&nbm==3&&nbl==3)";
    c_4b = "(nbt>=2&&nbm>=3&&nbl>=4)";
    hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
    res_abcd = {c_2b+"&&"+sbd, c_2b+"&&"+hig, c_3b+"&&"+sbd, c_3b+"&&"+hig, c_4b+"&&"+sbd, c_4b+"&&"+hig};

    vc_met.push_back("met>0   && met<=75");
    vc_met.push_back("met>75  && met<=150");
    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300");
    if (!bin_met) vc_met = {"met>0"};
  } else {
    c_2b = "(nbt==2&&nbm==2)";
    c_3b = "(nbt>=2&&nbm==3&&nbl==3)";
    c_4b = "(nbt>=2&&nbm>=3&&nbl>=4)";
    hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
    res_abcd = {c_2b+"&&"+sbd, c_2b+"&&"+hig, c_3b+"&&"+sbd, c_3b+"&&"+hig, c_4b+"&&"+sbd, c_4b+"&&"+hig};

    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300 && met<=400");
    vc_met.push_back("met>400");
    if (!bin_met) vc_met = {"met>150"};
  }


  //        Define and fill tables
  //------------------------------------
  unsigned nbins(0);
  for (auto &imet: vc_met) {
    for (auto &idrmax: vc_drmax) {
      // Label
      tablerows.push_back(TableRow("$"+CodeToLatex((imet+"&&"+idrmax).Data())+"$"));
      for (auto &iabcd: res_abcd) {
        // Cut for row
        NamedFunc res_reg_ = imet+"&&"+idrmax+"&&"+iabcd;
        // Add row to table
        tablerows.push_back(TableRow("$"+CodeToLatex(iabcd)+"$", res_reg_, 0,0, weight));
        nbins++;
      }
    }
  }

}

void constructCutflowTable(vector<TableRow> & tablerows, 
  NamedFunc weight) {

  tablerows = vector<TableRow>{
  TableRow("No selection", 
    "met>0",0,0, weight),
  TableRow("$p_\\text{t}^\\text{miss}>150 \\text{ GeV}$", 
    "met>150",0,0,weight),
  TableRow("$N_\\text{vl}=0$, $4\\leq N_\\text{j}\\leq 5$", 
    "met>150 && nvlep==0 && njet>=4 && njet<=5",0,0,weight),
  TableRow("$N_\\text{b}\\geq 2$", 
    "met>150 && nvlep==0 && njet>=4 && njet<=5 && nbt>=2",0,0,weight),
  TableRow("$N_\\text{tk}=0$", 
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0",0,0,weight),
  TableRow("Fake MET Cuts",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2",0,0,weight),
  TableRow("$\\Delta m<40 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40",0,0,weight),
  TableRow("$\\Delta R_\\text{max}<2.2$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2",0,0,weight),
  TableRow("$100<\\langle m\\rangle <140 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140",0,0,weight),
  TableRow("$N_\\text{b}\\geq 3$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3",0,0,weight),
  TableRow("$N_\\text{b}=4$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,weight),
  TableRow("$p_\\text{T}^\\text{miss}>200 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>200 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,weight),
  TableRow("$p_\\text{T}^\\text{miss}>300 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>300 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,weight),
  TableRow("$p_\\text{T}^\\text{miss}>400 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>400 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,weight),
  //TableRow("$\\Delta R_\\text{max}>1.1$",        
  //  "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>450 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]>=1.1 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,weight),
  };

}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);
  time_t begtime, endtime;
  time(&begtime);

  NamedFunc search_resolved_cuts = 
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved_cuts = 
                         "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc zll_resolved_cuts =
                         "nlep==2&&njet>=4&&njet<=5&&met<50&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  NamedFunc qcd_resolved_cuts =
                         "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  //NamedFunc study_cut = Higfuncs::lead_signal_lepton_pt>20 && Higfuncs::lead_signal_lepton_pt<25;
  NamedFunc study_cut = "1";
  // table_mode = regions, cutflow
  string table_mode = "cutflow";
  // Only for regions table_mode
  bool bin_met = false;
  bool bin_drmax = false;
  // Should only be used when sample is search
  bool use_met150_skim = true;

  // dataType: mc/data/signal
  // fileTag: (tt, single_t, zjets, wjets, qcd, other, all), (all), (200_1, 600_1, 950_1, ...)
  // sample: search/ttbar/zll/qcd
  // year_string: 2016/2017/2018/run2
  //string production_a = "higgsino_eldorado"; string nanoAodFolder_a = "/net/cms29/cms29r0/pico/NanoAODv5";
  string production_a = "higgsino_humboldt"; string nanoAodFolder_a = "/net/cms25/cms25r5/pico/NanoAODv5";
  string sample_a = "search";
  string year_string_a = "2017";
  TString tablename_a = "FixName:table_resolved_"+production_a+"_"+sample_a+"_"+CopyReplaceAll(year_string_a, ",","_")+"_a";

  string total_luminosity_string_a = getLuminosityString(year_string_a);

  //    Define processes, including intersections
  //--------------------------------------------------
  //NamedFunc base_filters_a = HigUtilities::pass_2016 && "met/met_calo<5"; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters_a = HigUtilities::pass_2016 && "met/met_calo<5"; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters_a = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...
  NamedFunc base_filters_a = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  //NamedFunc weight_a = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years_scaleup;
  //NamedFunc weight_a = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years;
  NamedFunc weight_a = "weight"*Higfuncs::eff_higtrig*w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig_run2*w_years;
  //NamedFunc weight_a = "weight"*Higfuncs::eff_higtrig*w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years;


  string production_b = "higgsino_humboldt"; string nanoAodFolder_b = "/net/cms25/cms25r5/pico/NanoAODv5";
  string sample_b = "search";
  string year_string_b = "2017";
  TString tablename_b = "FixName:table_resolved_"+production_b+"_"+sample_b+"_"+CopyReplaceAll(year_string_b, ",","_")+"_b";

  string total_luminosity_string_b = getLuminosityString(year_string_b);

  //    Define processes, including intersections
  //--------------------------------------------------
  //NamedFunc base_filters_b = HigUtilities::pass_2016 && "met/met_calo<5"; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters_b = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters_b = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...
  NamedFunc base_filters_b = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  //NamedFunc weight_b = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years;
  //NamedFunc weight_b = "w_lumi*w_isr"*Higfuncs::eff_higtrig_run2*w_years;
  NamedFunc weight_b = "weight"*Higfuncs::eff_higtrig*w_years*Functions::w_pileup;
  //NamedFunc weight_b = "weight"*Higfuncs::eff_higtrig_run2*w_years;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Setup plotMaker    //////////////////////////////////////////

  NamedFunc resolved_cuts_a = "1";
  if (sample_a == "ttbar") resolved_cuts_a = ttbar_resolved_cuts;
  else if (sample_a == "zll") resolved_cuts_a = zll_resolved_cuts;
  else if (sample_a == "qcd") resolved_cuts_a = qcd_resolved_cuts;
  else resolved_cuts_a = search_resolved_cuts;

  resolved_cuts_a = resolved_cuts_a && study_cut;
  if (table_mode == "cutflow") resolved_cuts_a = study_cut;

  NamedFunc resolved_cuts_b = "1";
  if (sample_b == "ttbar") resolved_cuts_b = ttbar_resolved_cuts;
  else if (sample_b == "zll") resolved_cuts_b = zll_resolved_cuts;
  else if (sample_b == "qcd") resolved_cuts_b = qcd_resolved_cuts;
  else resolved_cuts_b = search_resolved_cuts;

  resolved_cuts_b = resolved_cuts_b && study_cut;
  if (table_mode == "cutflow") resolved_cuts_b = study_cut;

  vector<shared_ptr<Process> > procs_a;
  addAllMcProcesses("", base_filters_a&&resolved_cuts_a,
    nanoAodFolder_a, production_a, 
    sample_a, year_string_a,
    procs_a, use_met150_skim);
  if (doSignal) {
    addMultipleSignalProcesses("", base_filters_a&&resolved_cuts_a,
      nanoAodFolder_a, production_a, 
      sample_a, year_string_a,
      procs_a, use_met150_skim);
  }

  vector<shared_ptr<Process> > procs_b;
  addAllMcProcesses("", base_filters_b&&resolved_cuts_b,
    nanoAodFolder_b, production_b, 
    sample_b, year_string_b,
    procs_b, use_met150_skim);
  if (doSignal) {
    addMultipleSignalProcesses("", base_filters_b&&resolved_cuts_b,
      nanoAodFolder_b, production_b, 
      sample_b, year_string_b,
      procs_b, use_met150_skim);
  }

  vector<TableRow> tablerows_a;
  vector<TableRow> tablerows_b;
  if (table_mode == "regions") {
    constructRegionTable(tablerows_a, sample_a, weight_a,  bin_met, bin_drmax);
    constructRegionTable(tablerows_b, sample_b, weight_b,  bin_met, bin_drmax);
  } else if (table_mode == "cutflow") {
    constructCutflowTable(tablerows_a, weight_a);
    constructCutflowTable(tablerows_b, weight_b);
  }

  PlotMaker pm;
  bool print_uncertainty = true;
  pm.Push<Table>(tablename_a.Data(), tablerows_a, procs_a, 0, 1, 0, 1, 0, print_uncertainty).LuminosityTag(total_luminosity_string_a).Precision(3);
  pm.Push<Table>(tablename_b.Data(), tablerows_b, procs_b, 0, 1, 0, 1, 0, print_uncertainty).LuminosityTag(total_luminosity_string_b).Precision(3);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Compare tables  //////////////////////////////////////////

  string compareFilename = "tables/table_compare_resolved_"+production_a+"_"+sample_a+"_"+CopyReplaceAll(year_string_a, ",","_")+".tex";

  // Get values from tables from latex file
  string filename_a = "tables/"+CopyReplaceAll(tablename_a.Data(), "FixName:","")+".tex";
  // tableValues[iRow][iColumn]
  vector<vector<float> > tableValues_a;
  vector<vector<float> > tableUncertainties_a;
  getValuesFromTexTable(filename_a, tableValues_a, tableUncertainties_a);
  //cout<<"values"<<endl;
  //for(auto const & row : tableValues_a) {
  //  for(auto const & item: row) cout<<item<<" ";
  //  cout<<endl;
  //}
  //cout<<"uncertainties"<<endl;
  //for(auto const & row : tableUncertainties_a) {
  //  for(auto const & item: row) cout<<item<<" ";
  //  cout<<endl;
  //}
  //tableValues_a[0][0] = 10000;
  //tableUncertainties_a[0][0] = 100;
  //replaceValuesInTexTablename(filename_a, compareFilename, tableValues_a, tableUncertainties_a, 3);
  string filename_b = "tables/"+CopyReplaceAll(tablename_b.Data(), "FixName:","")+".tex";
  vector<vector<float> > tableValues_b;
  vector<vector<float> > tableUncertainties_b;
  getValuesFromTexTable(filename_b, tableValues_b, tableUncertainties_b);

  // Compare table values and calculate uncertainty
  vector<vector<float> > compareTableValues(tableValues_a.size(), vector<float>(tableValues_a[0].size()));
  vector<vector<float> > compareTableUncertainties(tableUncertainties_a.size(), vector<float>(tableUncertainties_a[0].size()));
  for (unsigned iRow = 0; iRow < tableValues_a.size(); ++iRow) {
    for (unsigned iColumn = 0; iColumn < tableValues_a[iRow].size(); ++iColumn) {
      //cout<<"iRow: "<<iRow<<" iColumn: "<<iColumn<<endl;
      float value_a = tableValues_a[iRow][iColumn];
      float value_b = tableValues_b[iRow][iColumn];
      float uncertainty_a = tableUncertainties_a[iRow][iColumn];
      float uncertainty_b = tableUncertainties_b[iRow][iColumn];
      
      float newValue = 0;
      float newUncertainty = 0;
      if (value_a!=0)  {
        newValue = value_b/value_a;
        //newUncertainty = 2;
        newUncertainty = sqrt(pow(1/value_a*uncertainty_b,2)+pow(value_b/value_a/value_a*uncertainty_a,2));
      }
      else {
        newValue = 0;
        newUncertainty = 0;
      }

      compareTableValues[iRow][iColumn] = newValue;
      compareTableUncertainties[iRow][iColumn] = newUncertainty;
    }
  }

  replaceValuesInTexTablename(filename_a, compareFilename, compareTableValues, compareTableUncertainties, 2);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"no_signal", no_argument, 0, 'n'},    
      {"luminosity", required_argument, 0, 'l'},    
      {"tag", required_argument, 0, 't'},
      {"boo", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "nl:t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      doSignal = false;
      break;
    case 'l':
      lumi = atof(optarg);
      break;
    case 0:
      optname = long_options[option_index].name;
      if(0){
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
