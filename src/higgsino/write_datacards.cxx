#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <unistd.h> // getopt in Macs
#include <getopt.h>
#include <dirent.h>

#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TError.h" // Controls error level reporting
#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/table.hpp"
#include "core/palette.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/write_datacards.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
namespace 
{
  string outFolder = getenv("PWD");
  string mass_points_string = ""; // run over all the points, otherwise specify, e.g. "175_1,200_1"
  string years_string = "2016";
  float luminosity = 137.;
  // float luminosity = 35.9;
  string dimensionFilePath = "";
  bool unblind = false;
  string tag = "resolved";
  bool do_met_average = true;
  string higgsino_model = "N1N2";
}

const NamedFunc min_jet_dphi("min_jet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  for (size_t ijet = 0; ijet < (*b.jet_phi()).size(); ++ijet) {
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.jet_phi())[ijet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});

int main(int argc, char *argv[])
{
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);
  HigWriteDataCards::GetOptions(argc, argv);

  string baseFolder = ""; string hostName = execute("echo $HOSTNAME");
  if(Contains(hostName, "cms") || Contains(hostName, "compute-")) baseFolder = "/net/cms29";
  //baseFolder = "";

  cout<<"INFO:: Systematics are ON. Make sure to run on appropriate babies, i.e. unskimmed or skim_sys_abcd!!"<<endl;
  gSystem->mkdir(outFolder.c_str(), kTRUE);
  // massPoints = { {"1000","1"} }
  set<int> years;
  HigUtilities::parseYears(years_string, years);

  map<string, string> samplePaths;
  samplePaths["mc_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/merged_higmc_preselect/";
  samplePaths["signal_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  // samplePaths["data_2016"] = baseFolder + "/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";

  vector<pair<string, string> > massPoints;
  if (mass_points_string == "") HigUtilities::findMassPoints(samplePaths["signal_2016"], massPoints);
  else HigUtilities::parseMassPoints(mass_points_string, massPoints);

  NamedFunc filters = HigUtilities::pass_2016;

  // sampleProcesses[mc, data, signal]
  map<string, vector<shared_ptr<Process> > > sampleProcesses;
  HigUtilities::setMcProcesses(years, samplePaths, filters && "stitch", sampleProcesses);
  // Higfuncs::trig_hig will only work for 2016
  //HigUtilities::setDataProcesses(years, samplePaths, filters&&Higfuncs::trig_hig>0., sampleProcesses);
  HigUtilities::setDataProcesses(years, samplePaths, filters, sampleProcesses);
  HigUtilities::setSignalProcesses(massPoints, years, samplePaths, filters, sampleProcesses);
  
  NamedFunc weight = "w_lumi*w_isr";
  if (higgsino_model=="N1N2") weight *= HigUtilities::w_CNToN1N2;
  //weight *=  (min_jet_dphi>0.5);
  string baseline = "!low_dphi_met && nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200";
  if (tag=="resolved") baseline += "&& njet>=4 && njet<=5 && nbt>=2 && "+higtrim;
  else if (tag=="boosted") {
    baseline += " && ht>600 && nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300";
    baseline += " && fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 && fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250";
  }

  // Set ABCDbins
  // "bkg" and "sig" should be kept for setDataCardBackground()
  map<string, string> xBins;
  if (tag=="resolved") {
    xBins["bkg"] = "nbm==2";
    xBins["sig0"] = "nbm==3&&nbl==3";
    xBins["sig1"] = "nbm>=3&&nbl>=4";
  } else if (tag=="boosted") {
    string j1bb = "((met<=300 && fjet_mva_hbb_btv[0]>0.6)||(met>300 && fjet_mva_hbb_btv[0]>0.3))";
    string j2bb = "((met<=300 && fjet_mva_hbb_btv[1]>0.6)||(met>300 && fjet_mva_hbb_btv[1]>0.3))";
    xBins["bkg"] = "(!"+j1bb+"&& !"+j2bb+")";
    xBins["sig0"] = "((!"+j1bb+"&&"+j2bb+")||("+j1bb+"&& !"+j2bb+"))";
    xBins["sig1"] = j1bb+"&&"+j2bb;
  }
  
  map<string, string> yBins; // Shares sideband
  if (tag=="resolved") {
    yBins["sig"] = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    yBins["bkg"] = "!("+yBins["sig"]+")";
  } else if (tag=="boosted") {
    yBins["sig"] = "fjet_msoftdrop[0]>95 && fjet_msoftdrop[0]<=145 && fjet_msoftdrop[1]>95 && fjet_msoftdrop[1]<=145";
    yBins["bkg"] = "!("+yBins["sig"]+")";
  }
  
  map<string, vector<pair<string, string> > > dimensionBins;
  if (dimensionFilePath==""){
    if (tag=="resolved") {

      dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      dimensionBins["met"].push_back({"met2", "met>300 && met<=400"});
      dimensionBins["met"].push_back({"met3", "met>400"});
      dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});
    } else {
      dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      dimensionBins["met"].push_back({"met2", "met>300 && met<=500"});
      dimensionBins["met"].push_back({"met3", "met>500 && met<=700"});
      dimensionBins["met"].push_back({"met4", "met>700"});
    }
  } else {
    HigWriteDataCards::readDimensionFile(dimensionFilePath, dimensionBins);
    for (auto & dimension : dimensionBins){
      cout<<"dimension bins"<<endl;
      cout<<dimension.first<<endl;
      for (auto & entry : dimension.second){
        cout<<" "<<entry.first<<" "<<entry.second<<endl;
      }
    }
  }


  // sampleBins = { {label, cut} }
  vector<pair<string, string> > sampleBins;
  HigUtilities::setABCDBins(xBins, yBins, dimensionBins, sampleBins);
  //sampleBins = {{"test","1"}};

  // Set systematics somehow..
  // Consider weight for nonHH somehow..

  // cuts[mc,data,signal] = RowInformation(labels, tableRows, yields)
  map<string, HigUtilities::RowInformation > cutTable;
  //HigUtilities::addBinCuts(sampleBins, baseline, weight, "data", cutTable["data"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signal", cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signalGenMet", HigUtilities::nom2genmet, cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "mc", cutTable["mc"]);

  PlotMaker pm;
  // Luminosity used for labeling for table
  // Luminosity used for scaling for hist1d
  bool verbose = false;
  HigUtilities::makePlots(cutTable, sampleProcesses, luminosity, pm, verbose);

  // mYields[process_tag_sampleBinLabel] = GammaParams, TableRow
  map<string, pair<GammaParams, TableRow> > mYields;
  //HigUtilities::fillDataYields(pm, cutTable["data"], mYields);
  // Luminosity used for scaling
  HigUtilities::fillMcYields(pm, luminosity, cutTable["mc"], mYields, true);
  // Luminosity used for scaling
  HigUtilities::fillSignalYieldsProcesses(pm, luminosity, sampleProcesses["signal"], cutTable["signal"], mYields);
  HigUtilities::fillAverageGenMetYields(sampleProcesses["signal"], sampleBins, "signal", "signalGenMet", "signalAverageGenMet", mYields);

  cout<<"cut name: "<<cutTable["signal"].tableRows[5].cut_.Name()<<endl;

  for (auto & process : sampleProcesses["signal"]){
    cout<<"Process name: "<<process->name_<<endl;
    string model;
    int mGluino=0, mLSP=0;
    HigUtilities::getInfoFromProcessName(process->name_, model, mGluino, mLSP);
    TString outPath = outFolder+"/datacard-"+HigUtilities::setProcessNameLong(model, mGluino, mLSP)+"_"+years_string+"_"+tag+".txt";
    cout<<"open "<<outPath<<endl;
    ofstream cardFile(outPath);
    HigWriteDataCards::writeDataCardHeader(sampleBins,cardFile);

    vector<vector<string> > tableValues;
    if (unblind) HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "data", tableValues);
    else HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "mc", tableValues);
    if (do_met_average) {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
    } else {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signal", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signal", mYields, sampleBins, tableValues);      
    }
    HigWriteDataCards::writeTableValues(tableValues,cardFile);
    tableValues.clear();
    HigWriteDataCards::setDataCardBackground(mYields, sampleBins, "mc", tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile, true);
    tableValues.clear();
    //writeDataCardClosure()
    //writeDataCardUncertainties()
    //writeDataCardParam()
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

namespace HigWriteDataCards{
  void readDimensionFile(std::string const & dimensionFilePath, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins){
    map<string, int> dimensionIndex;
    ifstream dimensionFile(dimensionFilePath.c_str());
    if (dimensionFile.is_open()){
      string dimension, dimensionCut;
      while (dimensionFile >> dimension >> dimensionCut){
        dimensionBins[dimension].push_back(std::make_pair(dimension+to_string(dimensionIndex[dimension]), dimensionCut));
        dimensionIndex[dimension]++;
      }
      dimensionFile.close();
    }
  }

  void writeDataCardHeader(vector<pair<string, string> > sampleBins, ofstream & cardFile)
  {
    cardFile<<"imax "<<sampleBins.size()<<"  number of channels\n";
    cardFile<<"jmax 1  number of backgrounds\n";
    cardFile<<"kmax *  number of nuisance parameters\n";
    cardFile<<"shapes * * FAKE\n\n";
  }

  void setDataCardObserved(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, string const & dataTag, vector<vector<string> > & tableValues){
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    row[0] = "bin";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+3] = sampleBins[iBin].first;
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "Observation";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+3] = to_string(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield());
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }

  void setDataCardSignalBackground(string const & processName, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues){
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    row[0] = "bin";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = sampleBins[iBin].first;
      row[iBin*2+3] = sampleBins[iBin].first;
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "process";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = "sig";
      row[iBin*2+3] = "bkg";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "process";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = "0";
      row[iBin*2+3] = "1";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "rate";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      string label = processName + "_"+signalAverageGenMetTag+"_" +sampleBins[iBin].first;
      row[iBin*2+2] = to_string(mYields.at(label).first.Yield());
      row[iBin*2+3] = "1";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }
  
  void setDataCardSignalStatistics(string const & processName, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues)
  {
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    string label;
  
    for (unsigned rBin=0;rBin<sampleBins.size();++rBin)
    {
      setRow(row,"-");
      row[0] = "stat_"+sampleBins[rBin].first;
      row[1] = "lnN";
      label = processName + "_"+signalAverageGenMetTag+"_" +sampleBins[rBin].first;
      string uncertainty;
      if (mYields.at(label).first.Yield() == 0) uncertainty = "2.00";
      else uncertainty = to_string(1+mYields.at(label).first.Uncertainty()/mYields.at(label).first.Yield());
      row[2+2*rBin] = uncertainty;
      tableValues.push_back(row);
    }
  
    setRow(row,"");
    tableValues.push_back(row);
  }
  
  // returns count of non-overlapping occurrences of 'sub' in 'str'
  int countSubstring(const std::string& str, const std::string& sub)
  {
    if (sub.length() == 0) return 0;
    int count = 0;
    for (size_t offset = str.find(sub); offset != std::string::npos; offset = str.find(sub, offset + sub.length()))
    {
        ++count;
    }
    return count;
  }

  void setDataCardBackground(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, string const & mcTag, vector<vector<string> > & tableValues)
  {
    // title + type + tag + tagType + (yield,equation), equation arguments
    vector<string> row(6);
    for (unsigned rBin=0;rBin<sampleBins.size();++rBin)
    {
      row[0] = "rp_"+sampleBins[rBin].first;
      row[1] = "rateParam";
      row[2] = sampleBins[rBin].first;
      row[3] = "bkg";
      if (countSubstring(sampleBins[rBin].first,"sig") != 2) 
      {
        row[4] = to_string(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield());
        // Infinity defined in https://root.cern.ch/doc/master/RooNumber_8cxx_source.html
        row[5] = "[0,1.0e30]";
      }
      else 
      {
        row[4] = "(@0*@1/@2)";
        vector<string> xydimension;
        HigUtilities::stringToVectorString(sampleBins[rBin].first, xydimension, "_");
        vector<string> aNameVector = xydimension;
        aNameVector[0] = "xbkg";
        aNameVector[1] = "ybkg";
        string aName;
        HigUtilities::vectorStringToString(aNameVector, aName, "_");
        vector<string> bNameVector = xydimension;
        bNameVector[1] = "ybkg";
        string bName;
        HigUtilities::vectorStringToString(bNameVector, bName, "_");
        vector<string> cNameVector = xydimension;
        cNameVector[0] = "xbkg";
        string cName;
        HigUtilities::vectorStringToString(cNameVector, cName, "_");
        //cout<<aName<<" "<<bName<<" "<<cName<<endl;
        row[5] = "rp_"+cName+",rp_"+bName+",rp_"+aName;
      }
      tableValues.push_back(row);
      setRow(row,"");
    }
  
  }
  
  void writeTableValues(vector<vector<string> > & tableValues, ofstream & cardFile, bool alignLeft)
  {
    // Set space values
    vector<unsigned> xSpace(tableValues.back().size());
    for (unsigned yIndex = 0; yIndex < tableValues.size(); ++yIndex)
    {
      for (unsigned xIndex = 0; xIndex < tableValues[yIndex].size(); ++xIndex)
      {
        if (xSpace[xIndex] < tableValues[yIndex][xIndex].size()) xSpace[xIndex] = tableValues[yIndex][xIndex].size();
      }
    }
  
    // Write values
    for (unsigned yIndex = 0; yIndex < tableValues.size(); ++yIndex)
    {
      for (unsigned xIndex = 0; xIndex < tableValues[yIndex].size(); ++xIndex)
      {
        if (alignLeft) cardFile << setw(xSpace[xIndex]) << std::left << tableValues[yIndex][xIndex] << " ";
        else cardFile << setw(xSpace[xIndex]) << tableValues[yIndex][xIndex] << " ";
      }
      cardFile<<endl;
    }
  }
  
  void setRow(vector<string> & row, string const & value)
  {
    for(auto & item : row) item = value;
  }

  void GetOptions(int argc, char *argv[])
  {
    string blah;
    while(true){
      static struct option long_options[] = {
        {"output_folder", required_argument, 0, 'o'},
        {"mass_points", required_argument, 0, 'p'},
        {"years", required_argument, 0, 'y'},
        {"luminosity", required_argument, 0, 'l'},
        {"dimension", required_argument, 0, 'd'},
        {"tag", required_argument, 0, 't'},
        {"unblind", no_argument, 0, 'u'},
        {"recomet", no_argument, 0, 0},
        {"higgsino_model", required_argument, 0, 'm'},
        {0, 0, 0, 0}
      };
  
      char opt = -1;
      int option_index;
      opt = getopt_long(argc, argv, "o:p:y:l:d:t:m:nu", long_options, &option_index);
      if( opt == -1) break;
  
      string optname;
      switch(opt){
        case 'o': outFolder = optarg; break;
        case 'p': mass_points_string = optarg; break;
        case 'y': years_string = optarg; break;
        case 'l': luminosity = atof(optarg); break;
        case 't': tag = optarg; break;
        case 'm': higgsino_model = optarg; break;
        case 'd': 
          dimensionFilePath = optarg; 
          if (!FileExists(dimensionFilePath)) 
          {
            cout<<"[Error] No file caled "<<dimensionFilePath<<endl;
            exit(EXIT_FAILURE);
          }
          break;
        case 'u':
          unblind = true;
          break;
        case 0:
          optname = long_options[option_index].name;
          if(optname == "recomet"){
            do_met_average = false;
          }else{
            printf("Bad option! Found option name %s\n", optname.c_str());
          }
          break;
        default: 
          printf("Bad option! getopt_long returned character code 0%o\n", opt); 
          break;
      }
    }
  }

}
