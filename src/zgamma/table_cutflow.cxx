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

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  string year_string = "2017";
  int lepton_int = 2; // 0: electron, 1: muon, 2: el+mu
  float lumi = 1;
  string tag = "";
}

std::string removeSpaces(std::string inString){
  inString.erase(std::remove(inString.begin(),inString.end(),' '),inString.end());
  return inString;
}
int stringToVectorString(std::string const& inString, std::vector<std::string>& outputVector, std::string const & delimiter)
{
    if( outputVector.size() != 0) {
    std::cout<<"[Error] StringFunctions::divideString() => outputVector size is not 0."<<std::endl;
    return 1;
  }

  if(delimiter.size()==0){
    std::cout<<"[Error] StringFunctions::divideString() => delimiter size is 0."<<std::endl;
    return 1;
  }

  size_t start = 0;
  size_t end = 0;
  while((end=inString.find(delimiter,start)) != string::npos) {
    if(start!=end) {
      //std::cout<<"found: "<<removeSpaces(inString.substr(start, end-start))<<std::endl;
      outputVector.push_back(removeSpaces(inString.substr(start, end-start)));
    }
    start = end + delimiter.length();
  }
  if(start!=inString.length()) {
    //std::cout<<"found: "<<removeSpaces(inString.substr(start))<<std::endl;
    outputVector.push_back(removeSpaces(inString.substr(start)));
  }
  return 0;
}

NamedFunc add_cut(NamedFunc & current_cut, NamedFunc additional_cut) {
  current_cut = current_cut && additional_cut;
  return current_cut;
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
        stringToVectorString(valueString, splitString, "$\\pm$");
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

void parseYears(string years_string, set<string> & years)
{
  if (years_string == "run2") years = {"2016", "2016APV", "2017", "2018"};
  else if (years_string!="")
  {
    size_t found = years_string.find(",");
    size_t start = 0;
    string _tmp = "";
    while (found != string::npos) {
      _tmp = years_string.substr(start, found-start);
      years.insert(_tmp.c_str());
      start = found+1;
      found = years_string.find(",",start);
    } 
    _tmp = years_string.substr(start, found-start);
    years.insert(_tmp.c_str());
  }
}

string getLuminosityString(string const & year_string) {
  set<string> years;
  parseYears(year_string, years);
  float total_luminosity = 0;
  for (auto const & year : years) {
    // https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDatasetsUL2016
    if (year == "2016APV") total_luminosity += 19.5;
    if (year == "2016") total_luminosity += 16.8;
    //if (year == 2016) total_luminosity += 35.9;
    if (year == "2017") total_luminosity += 41.5;
    if (year == "2018") total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();
  return total_luminosity_string;
}

// Luminosities
const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleTypeString().Contains("-")) return 1.;
  if (b.SampleTypeString()=="2016") return 19.5;
  else if (b.SampleTypeString()=="2016APV") return 16.8;
  else if (b.SampleTypeString()=="2017") return 41.48;
  else if (b.SampleTypeString()=="2018") return 59.83;
  else return 0.;
});

const NamedFunc w_lumi("w_lumi", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleTypeString().Contains("-")) return 1.;
  else return b.w_lumi();
});

const NamedFunc signal_lead_muon_pt("signal_lead_muon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  for (unsigned iPart = 0; iPart<b.mu_pt()->size(); iPart++) {
    if (b.mu_sig()->at(iPart)) return b.mu_pt()->at(iPart);
  }
  return -999;
});
const NamedFunc signal_lead_electron_pt("signal_lead_electron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  for (unsigned iPart = 0; iPart<b.el_pt()->size(); iPart++) {
    if (b.el_sig()->at(iPart)) return b.el_pt()->at(iPart);
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
const NamedFunc stitch("stitch",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type() == 6200) return b.stitch_dy(); // DY
  else if (b.type() == 17100) return !b.stitch_dy(); // SM
  return 1;
});

void constructCutflowTable(vector<TableRow> & tablerows, NamedFunc weight, int electron_or_muon) {
    NamedFunc current_cut = "1";
    NamedFunc el_trigger = "ll_lepid[0]==11 && (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL)";
    NamedFunc mu_trigger = "ll_lepid[0]==13 && (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)";
    if (electron_or_muon == 0) { // el
      tablerows = vector<TableRow>{
        TableRow("No selection", add_cut(current_cut,"1"),0,0, weight),
        //TableRow("overlap veto", add_cut(current_cut, stitch),0,0, weight),
        TableRow("$N_e \\geq 2 \\text{ }Z(ee)$", add_cut(current_cut, "nel>=2&&ll_lepid[0]==11"),0,0, weight),
        TableRow("ee trigger", add_cut(current_cut, el_trigger),0,0, weight),
        TableRow("$p_{T}^{\\text{lead }e} \\geq 25\\text{ GeV}$", add_cut(current_cut, signal_lead_electron_pt > 25),0,0, weight),
        TableRow("$p_{T}^{\\text{sublead }e} \\geq 15\\text{ GeV}$", add_cut(current_cut, signal_sublead_electron_pt > 15),0,1, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$m_{ee} \\geq 50 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 50"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ee\\gamma} \\geq 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>=(15./110)"),0,0, weight),
        TableRow("$m_{ee\\gamma}+m_{ee} \\geq 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0]+ll_m[0]) >= 185"),0,0, weight),
        TableRow("$100 < m_{ee\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
      };
    } else if (electron_or_muon == 1) { // mu
      tablerows = vector<TableRow>{
        TableRow("No selection", add_cut(current_cut,"1"),0,0, weight),
        //TableRow("overlap veto", add_cut(current_cut, stitch),0,0, weight),
        TableRow("$N_{\\mu} \\geq 2 \\text{ }Z(\\mu\\mu)$", add_cut(current_cut, "nmu>=2&&ll_lepid[0]==13"),0,0, weight),
        TableRow("$\\mu\\mu$ trigger", add_cut(current_cut, mu_trigger),0,0, weight),
        TableRow("$p_{T}^{\\text{lead }\\mu} \\geq 20\\text{ GeV}$", add_cut(current_cut, signal_lead_muon_pt > 20),0,0, weight),
        TableRow("$p_{T}^{\\text{sublead }\\mu} \\geq 10\\text{ GeV}$", add_cut(current_cut, signal_sublead_muon_pt > 10),0,1, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$m_{\\mu\\mu} \\geq 50 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 50"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{\\mu\\mu\\gamma} \\geq 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>=(15./110)"),0,0, weight),
        TableRow("$m_{\\mu\\mu\\gamma}+m_{\\mu\\mu} \\geq 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0]+ll_m[0]) >= 185"),0,0, weight),
        TableRow("$100 < m_{\\mu\\mu\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
      };
    } else { // mu + el
      tablerows = vector<TableRow>{
        TableRow("No selection", add_cut(current_cut,"1"),0,0, weight),
        //TableRow("overlap veto", add_cut(current_cut, stitch),0,0, weight),
        TableRow("$N_e \\geq 2 \\text{ }Z(ee) || N_{\\mu} \\geq 2 \\text{ }Z(\\mu\\mu)$", add_cut(current_cut, "(nel>=2&&ll_lepid[0]==11) || (nmu>=2&&ll_lepid[0]==13)"),0,0, weight),
        TableRow("ee trigger (Z(ee)) $||$ $\\mu\\mu$ trigger (Z($\\mu\\mu$))", add_cut(current_cut, el_trigger || mu_trigger),0,0, weight),
        TableRow("$p_{T}^{\\text{lead }e} \\geq 25\\text{ GeV} (Z(ee))|| p_{T}^{\\text{lead }\\mu} \\geq 20\\text{ GeV} (Z(\\mu\\mu))$", add_cut(current_cut, ("ll_lepid[0]==13"&&signal_lead_muon_pt > 20) || ("ll_lepid[0]==11"&&signal_lead_electron_pt > 25)),0,0, weight),
        TableRow("$p_{T}^{\\text{sublead }e} \\geq 15\\text{ GeV} (Z(ee))|| p_{T}^{\\text{sublead }\\mu} \\geq 10\\text{ GeV} (Z(\\mu\\mu))$", add_cut(current_cut, ("ll_lepid[0]==13"&&signal_sublead_muon_pt > 10) || ("ll_lepid[0]==11"&&signal_sublead_electron_pt > 15)),0,1, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$m_{ll} \\geq 50 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 50"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ll\\gamma} \\geq 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>=(15./110)"),0,0, weight),
        TableRow("$m_{ll\\gamma}+m_{ll} \\geq 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0]+ll_m[0]) >= 185"),0,0, weight),
        TableRow("$100 < m_{ll\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
      };
    }
}


int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);
  time_t begtime, endtime;
  time(&begtime);

  //NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  //NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  //NamedFunc baseline_cuts = (el_trigs||mu_trigs) && "(nmu>=2||nel>=2)&&nphoton>=1&&ll_m[0] > 50" && (signal_lead_muon_pt > 25 || signal_lead_electron_pt > 25) && "(photon_pt[0]/llphoton_m[0])>=(15./110) && (llphoton_m[0]+ll_m[0]) >= 185 && llphoton_m[0] > 100 && llphoton_m[0] < 180"; 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  bool doSignal = true;
  bool doMc = false;
  bool doData = false;
  string table_mode = "cutflow";

  //string path_to_production_a = "/homes/jbkim/analysis/nano2pico.htozgamma_deathvalley/unit_test_htozgamma_nanoaodv9";
  //string path_to_production_a = "/data/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v2/";
  string path_to_production_a = "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_deathvalley_v3/";
  //string path_to_production_a = "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_deathvalley_v3/";
  //string tag_a = "htozgamma_deathvalley_v2"; 
  string tag_a = "htozgamma_deathvalley_v3"; 
  if (lepton_int==0) tag_a += "_electron";
  else if (lepton_int==1) tag_a += "_muon";
  else  tag_a += "_lepton";
  //string skim_a = "merged_*_llg";
  string skim_a = "unskimmed";
  string year_string_a = year_string;
  set<string> years_a; parseYears(year_string_a, years_a);
  TString tablename_a = "FixName:table_"+tag_a+"_"+CopyReplaceAll(year_string_a, ",","_");
  string total_luminosity_string_a = getLuminosityString(year_string_a);

  //    Define processes, including intersections
  //--------------------------------------------------
  NamedFunc base_filters_a = "1";
  NamedFunc weight_a = "1";
  //NamedFunc weight_a = w_lumi*w_years;
  NamedFunc resolved_cuts_a = "1";
  // 0: electron, 1: muon, 2: el+mu
  //int electron_or_muon_a = 2;
  int electron_or_muon_a = lepton_int;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Setup plotMaker    //////////////////////////////////////////

  Palette colors("txt/colors.txt", "default");
  set<string> pathNames; 
  vector<shared_ptr<Process> > procs_a;
  if (doMc) {
    pathNames = attach_folder(path_to_production_a, years_a, "mc/"+skim_a, {"*DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8*.root"});
    //procs_a.push_back(Process::MakeShared<Baby_pico>("DYJets", Process::Type::background, colors("dy"), pathNames, "stitch_dy"));
    procs_a.push_back(Process::MakeShared<Baby_pico>("DYJets", Process::Type::background, colors("dy"), pathNames, "1"));
    //pathNames = attach_folder(path_to_production_a, years_a, "mc/"+skim_a, {"*ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX*.root"});
    ////procs_a.push_back(Process::MakeShared<Baby_pico>("ZG", Process::Type::background, colors("single_t"), pathNames, "!stitch_dy"));
    //procs_a.push_back(Process::MakeShared<Baby_pico>("ZG", Process::Type::background, colors("single_t"), pathNames, "1"));
  }
  if (doSignal) {
    pathNames = attach_folder(path_to_production_a, years_a, "mc/"+skim_a, {"*GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8*.root"});
    procs_a.push_back(Process::MakeShared<Baby_pico>("ggH(HToZG)", Process::Type::signal, colors("t1tttt"), pathNames, "1"));
  }
  if (doData) {
    string t_skim_a = (skim_a == "unskimmed") ? "raw_pico" : skim_a;
    pathNames = attach_folder(path_to_production_a, years_a, "data/"+t_skim_a, {"*DoubleEG*.root", "*EGamma*.root", "*DoubleMuon*.root"});
    procs_a.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, colors("data"), pathNames, "1"));
  }

  vector<TableRow> tablerows_a;
  constructCutflowTable(tablerows_a, weight_a, electron_or_muon_a);

  PlotMaker pm;
  bool print_uncertainty = false;
  pm.Push<Table>(tablename_a.Data(), tablerows_a, procs_a, 0, 1, 0, 1, 0, print_uncertainty).LuminosityTag(total_luminosity_string_a).Precision(3);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"year", required_argument, 0, 'y'},    
      {"lepton", required_argument, 0, 'l'},    
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "y:l:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'y':
      year_string = optarg;
      break;
    case 'l':
      lepton_int = atoi(optarg);
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
