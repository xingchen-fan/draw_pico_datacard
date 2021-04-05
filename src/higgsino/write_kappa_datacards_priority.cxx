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
#include <regex>
#include <utility>

#include "TH1.h"
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
#include "core/hist1d.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/write_kappa_datacards.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
namespace 
{
  string outFolder = getenv("PWD");
  string mass_points_string = ""; // run over all the points, otherwise specify, e.g. "175_1,200_1"
  string mass_point_glob = "";
  string years_string = "2016,2017,2018";
  float luminosity = 1.;
  //float luminosity = 35.9;
  string dimensionFilePath = "";
  bool unblind = false;
  bool unblind_signalregion = false;
  string tag = "resolved";
  bool do_met_average = true;
  string higgsino_model = "N1N2";
  // 1. RSR, RCR, BSR, BCR
  // 2. BSR, BCR, RSR, RCR
  // 3. RSR, BSR, RCR, BCR
  // 4. BSR, RSR, BCR, RCR
  int priority = 1;
  bool is1D = true;
  int t5hh_range = 0;
}

const NamedFunc min_jet_dphi("min_jet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  for (size_t ijet = 0; ijet < (*b.jet_phi()).size(); ++ijet) {
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.jet_phi())[ijet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});

string regSearch(string inString, string const & regExPattern) {
    std::smatch match;
    std::regex pattern(regExPattern);
    std::regex_search (inString, match, pattern);
    return match.str();
}

string getPartialDatasetName(string filename) {
  string baseFilename = filename.substr(filename.find_last_of("/\\")+1);
  string fullDatasetName = regSearch(baseFilename, "[A-Z].*|ttHTobb.*");
  //cout<<datasetName.substr(0,datasetName.find_last_of("Tune"))<<endl;;
  string datasetName = fullDatasetName;
  if (datasetName.find("Run") != string::npos) datasetName = "DATA";
  if (datasetName.find("SMS-TChiHH") != string::npos) datasetName = "TChiHH";
  if (datasetName.find("Tune") != string::npos) datasetName = datasetName.substr(0, datasetName.find("Tune")+4);
  if (datasetName.find("_13TeV") != string::npos) datasetName = datasetName.substr(0, datasetName.find("_13TeV"));
  return datasetName;
}

bool findStringIC(const std::string & strHaystack, const std::string & strNeedle)
{
  auto it = std::search(
    strHaystack.begin(), strHaystack.end(),
    strNeedle.begin(),   strNeedle.end(),
    [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); }
  );
  return (it != strHaystack.end() );
}

pair<int, int> getSignalMassValues(string filename) {
  pair<int, int> nlsp_lsp_mass;
  string baseFilename = filename.substr(filename.find_last_of("/\\")+1);
  string datasetName = regSearch(baseFilename, "[A-Z].*|ttHTobb.*");
  if (datasetName.find("SMS-TChiHH") == string::npos) return {-1,-1};
  string NLSPMassString = regSearch(datasetName, "mChi-.*?_").substr(5);
  string LSPMassString = regSearch(datasetName, "mLSP-.*?_").substr(5);
  //cout<<NLSPMassString<<" "<<LSPMassString<<endl;
  return {stoi(NLSPMassString),stoi(LSPMassString)};
}

bool regionCut(const Baby & b, int regionIndex) {
  bool inRegion = false;
  vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = static_cast<vector<map<Long64_t, set< tuple<string, int, int, int, int, int, int> > > > *> (b.EventVetoData());
  // 0: sampleType, 1: year, 2: run, 3: lumiblock, 4: met, 5: nlsp_mass, 6: lsp_mass
  // Check event nuber
  if ((*eventNumberData)[regionIndex].count(b.event())==1) {
    for (auto const & it : (*eventNumberData)[regionIndex][b.event()]) {
      // Check year
      if (get<1>(it) != abs(b.SampleType())) continue;
      // Check run
      if (get<2>(it) != b.run()) continue;
      bool isSignal = ((*b.FileNames().begin()).find("TChiHH") != string::npos) ? true : false;
      // Check met out of 20%
      if (isSignal) {if (get<4>(it) < b.met()*0.8 || get<4>(it) > b.met()*1.2) continue;}
      else if (get<4>(it) != round(b.met())) continue;
      pair<int, int> signal_mass = getSignalMassValues(*b.FileNames().begin());
      // Check nlsp_mass
      if (signal_mass.first != get<5>(it)) continue;
      // Check lsp_mass
      int regionLsp = get<6>(it)==1? 0:get<6>(it);
      if (signal_mass.second != regionLsp) continue;
      // Check sample name
      //if (get<5>(it) != -1) cout<<*b.FileNames().begin()<<" "<<get<0>(it)<<endl;
      if ((*b.FileNames().begin()).find("TChiHH") != string::npos && get<0>(it).find("TChiHH") != string::npos)  {
        if ((*b.FileNames().begin()).find("HToBB_2D") != get<0>(it).find("2D")) continue;
      } else if (!findStringIC(*b.FileNames().begin(), get<0>(it))) continue;

      //if (b.met()<300) cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
      //cout<<"sampleType: "<<get<0>(it)<<" year: "<<get<1>(it)<<" run: "<<get<2>(it)<<" lumiblock: "<<get<3>(it)<<" met: "<<b.event()<<get<4>(it)<<" nlsp_mass: "<<b.event()<<get<5>(it)<<" lsp_mass: "<<b.event()<<get<6>(it)<<endl;

      inRegion = true;
      break;
    }
  }
  return inRegion;
}

const NamedFunc boostSignalRegion("boostSignalRegion",[](const Baby &b) -> NamedFunc::ScalarType{
  return regionCut(b, 0);
});
const NamedFunc boostControlRegion("boostControlRegion",[](const Baby &b) -> NamedFunc::ScalarType{
  return regionCut(b, 1);
});

void addEventNumberData(set<string> eventNumberFilenames, vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData) {
  eventNumberData->push_back(map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > ());
  for (string const & eventNumberFilename : eventNumberFilenames) {
    ifstream eventNumberFile(eventNumberFilename);
    string line;
    while (getline(eventNumberFile, line)) {
      if(line[0]=='#') continue;
      if(line.find("tree") != string::npos) continue;
      vector<string> lineSplit;
      HigUtilities::stringToVectorString(line, lineSplit, ",");
      string const & sampleType = lineSplit[0];
      int year = stoi(lineSplit[1]);
      int run = stoi(lineSplit[2]);
      int lumiblock = stoi(lineSplit[3]);
      Long64_t event = stoll(lineSplit[4]);
      int met = stoi(lineSplit[5]);
      int nlsp_mass = -1; int lsp_mass = -1;
      bool isSignal = false; if (sampleType.find("TChiHH") != string::npos) isSignal = true;
      if (isSignal) {nlsp_mass = stoi(lineSplit[6]); lsp_mass = stoi(lineSplit[7]);}
      // Insert data
      if ((*eventNumberData).back().count(event) == 0) {
        //cout<<"sampleType: "<<sampleType<<" year: "<<year<<" run: "<<run<<" lumiblock: "<<lumiblock<<" met: "<<met<<" nlsp_mass: "<<nlsp_mass<<" lsp_mass: "<<lsp_mass<<endl;
        (*eventNumberData).back()[event] = {{sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass}};
      } else {
        if ((*eventNumberData).back()[event].count({sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass}) != 0) cout<<"Duplicate: "<<sampleType<<" "<<year<<" "<<run<<" "<<lumiblock<<" "<<event<<" "<<met<<" "<<nlsp_mass<<" "<<lsp_mass<<endl;
        (*eventNumberData).back()[event].insert({sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass});
      }
    }
    eventNumberFile.close();
  }

  //// Read out map
  //int regionIndex = 1;
  //for (auto const & regionIt : (*eventNumberData)[regionIndex]) {
  //  Long64_t const & event = regionIt.first;
  //  for (auto const & it : (*eventNumberData)[regionIndex][event]) {
  //      cout<<"sampleType: "<<get<0>(it)<<" year: "<<get<1>(it)<<" run: "<<get<2>(it)<<" lumiblock: "<<get<3>(it)<<" met: "<<get<4>(it)<<" nlsp_mass: "<<get<5>(it)<<" lsp_mass: "<<get<6>(it)<<endl;
  //  }
  //}

}


const NamedFunc weight_ht_sideband_2016("weight_ht_sideband_2016", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.492;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.111;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.044;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.019;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.976;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.982;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.889;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.836;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.859;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.723;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.871;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.802;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.747;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.791;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.852;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.802;
  else if (b.ht()>900 && b.ht()<=950) weight = 1.012;
  else if (b.ht()>950) weight = 0.971;
  return weight;
});
const NamedFunc weight_ht_sideband_2017("weight_ht_sideband_2017", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.5;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.164;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.098;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.051;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.977;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.897;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.775;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.779;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.741;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.697;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.717;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.657;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.661;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.718;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.603;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.8;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.426;
  else if (b.ht()>950) weight = 0.764;
  return weight;
});
const NamedFunc weight_ht_sideband_2018("weight_ht_sideband_2018", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.629;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.251;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.136;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.039;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.956;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.887;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.738;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.704;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.687;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.716;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.686;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.681;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.667;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.61;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.82;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.723;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.47;
  else if (b.ht()>950) weight = 0.624;
  return weight;
});
const NamedFunc weight_ht_sideband("weight_ht_sideband", [](const Baby &b) -> NamedFunc::ScalarType {
  float weight = 1.;
  if (b.SampleType()==2016) weight = weight_ht_sideband_2016.GetScalar(b);
  else if (b.SampleType()==2017) weight = weight_ht_sideband_2017.GetScalar(b);
  else if (b.SampleType()==2018) weight = weight_ht_sideband_2018.GetScalar(b);
  return weight;
});


int main(int argc, char *argv[])
{
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);
  HigWriteDataCards::GetOptions(argc, argv);

  //---------------------------------------------------------------------------
  //                         load appropriate picos
  //---------------------------------------------------------------------------
  string baseFolder = ""; string hostName = execute("echo $HOSTNAME");
  if(Contains(hostName, "cms") || Contains(hostName, "compute-")) baseFolder = "/net/cms29";
  //baseFolder = "";

  cout<<"INFO:: Systematics are ON. Make sure to run on appropriate babies, i.e. unskimmed or skim_higsys!!"<<endl;
  gSystem->mkdir(outFolder.c_str(), kTRUE);

  set<int> years;
  HigUtilities::parseYears(years_string, years);
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  map<string, string> samplePaths;
  ////samplePaths["mc_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/merged_higmc_preselect/";
  ////samplePaths["signal_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  ////samplePaths["data_2016"] = baseFolder + "/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";
  //samplePaths["mc_2016"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/mc/merged_higmc_preselect/";
  //samplePaths["signal_2016"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2017"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/2017/mc/merged_higmc_preselect/";
  //samplePaths["signal_2017"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/2017/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2018"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/mc/merged_higmc_preselect/";
  //samplePaths["signal_2018"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/SMS-TChiHH_2D/merged_higmc_preselect/";

  //samplePaths["mc_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2016/mc/merged_higmc_preselect/";
  //samplePaths["signal_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2017/mc/merged_higmc_preselect/";
  //samplePaths["signal_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2017/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2018/mc/merged_higmc_preselect/";
  //samplePaths["signal_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2018/SMS-TChiHH_2D/merged_higmc_preselect/";

  samplePaths["mc_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2016/mc/merged_higmc_preselect/";
  samplePaths["signal_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v2/2016/SMS-TChiHH_2D_fastSimJmeCorrection/skim_higsys/";
  samplePaths["data_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2016/data/merged_higdata_preselect/";
  samplePaths["mc_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2017/mc/merged_higmc_preselect/";
  samplePaths["signal_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v2/2017/SMS-TChiHH_2D_fastSimJmeCorrection/skim_higsys/";
  samplePaths["data_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2017/data/merged_higdata_preselect/";
  samplePaths["mc_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2018/mc/merged_higmc_preselect/";
  samplePaths["signal_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v2/2018/SMS-TChiHH_2D_fastSimJmeCorrection/skim_higsys/";
  samplePaths["data_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2018/data/merged_higdata_preselect/";
  if (higgsino_model=="T5HH") {
    samplePaths["signal_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2016/SMS-T5qqqqZH_fastSimJmeCorrection/skim_higsys/";
    samplePaths["signal_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2017/SMS-T5qqqqZH_fastSimJmeCorrection/skim_higsys/";
    samplePaths["signal_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2018/SMS-T5qqqqZH_fastSimJmeCorrection/skim_higsys/";
    //samplePaths["signal_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v2/2016/SMS-T5qqqqZH_fastSimJmeCorrection/unskimmed/";
    //samplePaths["signal_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v2/2017/SMS-T5qqqqZH_fastSimJmeCorrection/unskimmed/";
    //samplePaths["signal_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v2/2018/SMS-T5qqqqZH_fastSimJmeCorrection/unskimmed/";
  }

  //// massPoints = { {"1000","1"} }
  vector<pair<string, string> > massPoints;
  //if (higgsino_model == "N1N2") {
  //  string signal_folder = samplePaths["signal_2016"];
  //  set<string> signal_files = Glob(signal_folder+"/*.root");
  //  string mLSP, mChi;
  //  for (string signal_file : signal_files) {
  //    HigUtilities::filenameToMassPoint(signal_file, mChi, mLSP);
  //    if (stoi(mChi)>800) continue; // Ingore points above 800 GeV
  //    massPoints.push_back({mChi, mLSP});
  //  }
  //}
  //else if (mass_points_string == "") HigUtilities::findMassPoints(samplePaths["signal_2016"], massPoints);
  //else HigUtilities::parseMassPoints(mass_points_string, massPoints);

  if (mass_points_string != "") {
    HigUtilities::parseMassPoints(mass_points_string, massPoints);
  } else if (mass_point_glob != "") {
    string signal_folder = samplePaths["signal_2016"];
    set<string> signal_files;
    signal_files = Glob(signal_folder+"/"+mass_point_glob);
    string mLSP, mChi;
    for (string signal_file : signal_files) {
      HigUtilities::filenameToMassPoint(signal_file, mChi, mLSP);
      if (!is1D && stoi(mChi)>800) continue; // Ingore points above 800 GeV for 2D
      massPoints.push_back({mChi, mLSP});
    }
  } else {
    string signal_folder = samplePaths["signal_2016"];
    set<string> signal_files;
    if (higgsino_model=="T5HH") {
      signal_files = Glob(signal_folder+"/*.root");
    } else if (is1D) {
      signal_files = Glob(signal_folder+"/*mLSP-0_*.root");
    } else {
      signal_files = Glob(signal_folder+"/*.root");
    }
    string mGluino, mLSP, mChi;
    for (string signal_file : signal_files) {
      if (higgsino_model=="T5HH") {
        HigUtilities::filenameToMassPoint(signal_file, mGluino, mChi, mLSP);
        if (stoi(mChi) != (stoi(mGluino)-50)) continue; // Ingore points outside the 2D scan
        if (t5hh_range >= 0 && stoi(mGluino)<t5hh_range) continue; //non-negative t5hh is lower bound
        if (t5hh_range < 0 && stoi(mGluino)>=(-1*t5hh_range)) continue; //negative t5hh is upper bound
        massPoints.push_back({mGluino, mLSP});
      }
      else {
        HigUtilities::filenameToMassPoint(signal_file, mChi, mLSP);
        if (!is1D && stoi(mChi)>800) continue; // Ingore points above 800 GeV for 2D
        massPoints.push_back({mChi, mLSP});
      }
    }
  }


  //NamedFunc filters = HigUtilities::pass_2016;
  //NamedFunc filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";
  NamedFunc filters = Higfuncs::final_pass_filters;

  // sampleProcesses[mc, data, signal]
  map<string, vector<shared_ptr<Process> > > sampleProcesses;
  HigUtilities::setMcProcesses(years, samplePaths, filters && "stitch", sampleProcesses);
  // Higfuncs::trig_hig will only work for 2016
  //HigUtilities::setDataProcesses(years, samplePaths, filters&&Higfuncs::trig_hig>0., sampleProcesses);
  NamedFunc met_triggers = Higfuncs::met_trigger;
  if(unblind) HigUtilities::setDataProcesses(years, samplePaths, filters&&met_triggers, sampleProcesses);
  if (higgsino_model=="T5HH") {
    HigUtilities::setSignalProcessesT5HH(massPoints, years, samplePaths, filters, sampleProcesses);
  } else {
    HigUtilities::setSignalProcesses(massPoints, years, samplePaths, filters, sampleProcesses);
  }


  //named func for studying effect of not having an excess in 3b MET 300-400 (multiply to weight)
  const NamedFunc zero_excess("zero_excess", [](const Baby &b) -> NamedFunc::ScalarType{
    //remove "excess" events in 3b low drmax met 300-400
    if (b.SampleType()>0) return 1.;
    if (b.njet() >= 4 && b.njet() <= 5 && b.nlep()==0 ) {
      if (b.ntk() == 0 && b.met() > 300 && b.met() < 400) {
        if (!b.low_dphi_met() && b.hig_cand_dm()->at(0) < 40 && b.hig_cand_drmax()->at(0)<1.1) {
          if (b.nbt() >=2 && b.nbm() >= 3 && b.hig_cand_am()->at(0)>100) {
            if (b.hig_cand_am()->at(0)<140) {
              return 0.;
            }
          }
        }
      }
    }
    return 1.;
  });
  
  //---------------------------------------------------------------------------
  //                         set weight, selection, bins
  //---------------------------------------------------------------------------
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years*Functions::w_pileup;
  //NamedFunc weight = Higfuncs::final_weight;
  NamedFunc weight = Higfuncs::final_weight;
  NamedFunc weight_genmet = "weight"*Higfuncs::eff_higtrig_mettru*Higfuncs::w_years*Functions::w_pileup;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years*Higfuncs::w_pileup_nosignal;
  //NamedFunc weight = Higfuncs::final_weight*weight_ht_sideband;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*Higfuncs::w_years;

  if (higgsino_model=="N1N2") {
    weight *= HigUtilities::w_CNToN1N2;
    weight_genmet *= HigUtilities::w_CNToN1N2;
  }
  string baseline = "!low_dphi_met && nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200 && met/mht<2 && met/met_calo<2&& weight<1.5";
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
      //// Old
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      //dimensionBins["met"].push_back({"met2", "met>300 && met<=450"});
      //dimensionBins["met"].push_back({"met3", "met>450"});

      // Original
      dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      dimensionBins["met"].push_back({"met2", "met>300 && met<=400"});
      dimensionBins["met"].push_back({"met3", "met>400"});
      dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// Slightly lower met binning
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=250"});
      //dimensionBins["met"].push_back({"met2", "met>250 && met<=325"});
      //dimensionBins["met"].push_back({"met3", "met>325"});
      //dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      //dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// 3 met bin
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      //dimensionBins["met"].push_back({"met2", "met>300"});
      //dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      //dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// mix binning
      //dimensionBins["mix"].push_back({"met0_drmax1", "met>150 && met<=200 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met1_drmax1", "met>200 && met<=300 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met2_drmax1", "met>300 && met<=400 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met3_drmax1", "met>400 && hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met0_drmax0", "met>150 && met<=200 &&hig_cand_drmax[0]<=1.1"});
      //dimensionBins["mix"].push_back({"met1_drmax0", "met>200 && met<=300 &&hig_cand_drmax[0]<=1.1"});
      //dimensionBins["mix"].push_back({"met2_drmax0", "met>300 &&hig_cand_drmax[0]<=1.1"});

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
  // 1. RSR, RCR, BSR, BCR
  // 2. BSR, BCR, RSR, RCR
  // 3. RSR, BSR, RCR, BCR
  // 4. BSR, RSR, BCR, RCR
  HigUtilities::setABCDBinsPriority(xBins, yBins, dimensionBins, sampleBins, priority);
  //sampleBins = {{"test","1"}};

  //---------------------------------------------------------------------------
  //                         set systematics
  //---------------------------------------------------------------------------
  // Set systematics somehow..
  // controlSystematics[xsig0_ybkg_met0_drmax0]["ttbar"] = systematic value
  map<string, map<string, float> > controlSystematics;
  HigWriteDataCards::setControlSystematics(controlSystematics);
  // Consider weight for nonHH somehow..
  
  NamedFunc weight_notrig = Higfuncs::final_weight_notrgeff;
  if (higgsino_model=="N1N2") weight_notrig *= HigUtilities::w_CNToN1N2;
  
  //vector holding systematic variations as a name and weights. May need to be extended to a class
  vector<pair<string, vector<NamedFunc>>> systematics_vector;
  systematics_vector.push_back(make_pair(string("LumiSyst"),
      vector<NamedFunc>({})));
  systematics_vector.push_back(make_pair(string("TrigSyst"),
      vector<NamedFunc>({weight_notrig*Higfuncs::eff_higtrig_run2_syst_up, 
      weight_notrig*Higfuncs::eff_higtrig_run2_syst_down})));
  systematics_vector.push_back(make_pair(string("SignalJetID"),
      vector<NamedFunc>({})));
  systematics_vector.push_back(make_pair(string("SignalBCTag"),
      vector<NamedFunc>({weight*"sys_bchig[0]/w_bhig", weight*"sys_bchig[1]/w_bhig"})));
  systematics_vector.push_back(make_pair(string("SignalBCTagFastSIM"),
      vector<NamedFunc>({weight*"sys_fs_bchig[0]/w_bhig", weight*"sys_fs_bchig[1]/w_bhig"})));
  systematics_vector.push_back(make_pair(string("SignalUDSGTag"),
      vector<NamedFunc>({weight*"sys_udsghig[0]/w_bhig", weight*"sys_udsghig[1]/w_bhig"})));
  systematics_vector.push_back(make_pair(string("SignalUDSGTagFastSIM"),
      vector<NamedFunc>({weight*"sys_fs_udsghig[0]/w_bhig", weight*"sys_fs_udsghig[1]/w_bhig"})));
  systematics_vector.push_back(make_pair(string("SignalPrefire"),
      vector<NamedFunc>({weight*"sys_prefire[0]/w_prefire", weight*"sys_prefire[1]/w_prefire"})));
  systematics_vector.push_back(make_pair(string("ISRSyst"),
      vector<NamedFunc>({weight*"sys_isr[0]/w_isr", weight*"sys_isr[1]/w_isr"})));
  systematics_vector.push_back(make_pair(string("SignalMETFastSIM"),
      vector<NamedFunc>({})));
  systematics_vector.push_back(make_pair(string("SignalPU"),
      vector<NamedFunc>({weight/Functions::w_pileup*Functions::w_syst_pileup_up,
      weight/Functions::w_pileup*Functions::w_syst_pileup_down})));
  systematics_vector.push_back(make_pair(string("SignalScale"),
      vector<NamedFunc>({weight*"sys_murf[0]", weight*"sys_murf[1]", weight*"sys_murf[2]", 
      weight*"sys_murf[3]",weight*"sys_murf[5]",weight*"sys_murf[6]",weight*"sys_murf[7]",
      weight*"sys_murf[8]"})));
  systematics_vector.push_back(make_pair(string("SignalJEC"),
      vector<NamedFunc>({weight,weight})));
  systematics_vector.push_back(make_pair(string("SignalJER"),
      vector<NamedFunc>({weight,weight})));

  vector<pair<string, vector<NamedFunc>>> systematics_vector_genmet;
  systematics_vector_genmet.push_back(make_pair(string("LumiSyst"),
      vector<NamedFunc>({})));
  systematics_vector_genmet.push_back(make_pair(string("TrigSyst"),
      vector<NamedFunc>({weight_notrig*Higfuncs::eff_higtrig_run2_syst_up_mettru, 
      weight_notrig*Higfuncs::eff_higtrig_run2_syst_down_mettru})));
  systematics_vector_genmet.push_back(make_pair(string("SignalJetID"),
      vector<NamedFunc>({})));
  systematics_vector_genmet.push_back(make_pair(string("SignalBCTag"),
      vector<NamedFunc>({weight_genmet*"sys_bchig[0]/w_bhig", weight_genmet*"sys_bchig[1]/w_bhig"})));
  systematics_vector_genmet.push_back(make_pair(string("SignalBCTagFastSIM"),
      vector<NamedFunc>({weight_genmet*"sys_fs_bchig[0]/w_bhig", weight_genmet*"sys_fs_bchig[1]/w_bhig"})));
  systematics_vector_genmet.push_back(make_pair(string("SignalUDSGTag"),
      vector<NamedFunc>({weight_genmet*"sys_udsghig[0]/w_bhig", weight_genmet*"sys_udsghig[1]/w_bhig"})));
  systematics_vector_genmet.push_back(make_pair(string("SignalUDSGTagFastSIM"),
      vector<NamedFunc>({weight_genmet*"sys_fs_udsghig[0]/w_bhig", weight_genmet*"sys_fs_udsghig[1]/w_bhig"})));
  systematics_vector_genmet.push_back(make_pair(string("SignalPrefire"),
      vector<NamedFunc>({weight_genmet*"sys_prefire[0]/w_prefire", weight_genmet*"sys_prefire[1]/w_prefire"})));
  systematics_vector_genmet.push_back(make_pair(string("ISRSyst"),
      vector<NamedFunc>({weight_genmet*"sys_isr[0]/w_isr", weight_genmet*"sys_isr[1]/w_isr"})));
  systematics_vector_genmet.push_back(make_pair(string("SignalMETFastSIM"),
      vector<NamedFunc>({})));
  systematics_vector_genmet.push_back(make_pair(string("SignalPU"),
      vector<NamedFunc>({weight_genmet/Functions::w_pileup*Functions::w_syst_pileup_up, 
      weight_genmet/Functions::w_pileup*Functions::w_syst_pileup_down})));
  systematics_vector_genmet.push_back(make_pair(string("SignalScale"),
      vector<NamedFunc>({weight_genmet*"sys_murf[0]", weight_genmet*"sys_murf[1]", weight_genmet*"sys_murf[2]", 
      weight_genmet*"sys_murf[3]",weight_genmet*"sys_murf[5]",weight_genmet*"sys_murf[6]",weight_genmet*"sys_murf[7]",
      weight_genmet*"sys_murf[8]"})));
  systematics_vector_genmet.push_back(make_pair(string("SignalJEC"),
      vector<NamedFunc>({weight_genmet,weight_genmet})));
  systematics_vector_genmet.push_back(make_pair(string("SignalJER"),
      vector<NamedFunc>({weight_genmet,weight_genmet})));

  //---------------------------------------------------------------------------
  //                           book tables and plots
  //---------------------------------------------------------------------------
  // cuts[mc,data,signal] = RowInformation(labels, tableRows, yields)
  map<string, HigUtilities::RowInformation > cutTable;
  if(unblind) HigUtilities::addBinCuts(sampleBins, baseline, weight, "data", cutTable["data"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signal", cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight_genmet, "signalGenMet", HigUtilities::nom2genmet, cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "mc", cutTable["mc"]);
  for (unsigned sys_idx = 0; sys_idx < systematics_vector.size(); sys_idx++) {
    pair<string, vector<NamedFunc>> sys = systematics_vector[sys_idx];
    for (unsigned wgt_idx = 0; wgt_idx < sys.second.size(); wgt_idx++) {
      //for JEC and JER systematics, switch out analysis variables
      if (sys.first == "SignalJEC") {
        string sys_name = sys.first+to_string(wgt_idx);
        HigUtilities::addBinCuts(HigUtilities::nom2sys_bins(sampleBins,to_string(wgt_idx+2)), 
                                 HigUtilities::nom2sys_string(baseline,to_string(wgt_idx+2)), weight, 
                                 "signal_"+sys_name, cutTable["signal"]);
        HigUtilities::addBinCuts(HigUtilities::nom2sys_bins(sampleBins, to_string(wgt_idx+2)), 
                                 HigUtilities::nom2sys_string(baseline,to_string(wgt_idx+2)), weight_genmet, 
                                 "signalGenMet_"+sys_name, HigUtilities::nom2genmet,
                                 cutTable["signal"]);
      }
      else if (sys.first == "SignalJER") {
        string sys_name = sys.first+to_string(wgt_idx);
        HigUtilities::addBinCuts(HigUtilities::nom2sys_bins(sampleBins,to_string(wgt_idx)), 
                                 HigUtilities::nom2sys_string(baseline,to_string(wgt_idx)), weight, 
                                 "signal_"+sys_name, cutTable["signal"]);
        HigUtilities::addBinCuts(HigUtilities::nom2sys_bins(sampleBins, to_string(wgt_idx)), 
                                 HigUtilities::nom2sys_string(baseline,to_string(wgt_idx)), weight_genmet, 
                                 "signalGenMet_"+sys_name, HigUtilities::nom2genmet,
                                 cutTable["signal"]);
      }
      else {
        string sys_name = sys.first+to_string(wgt_idx);
        HigUtilities::addBinCuts(sampleBins, baseline, sys.second[wgt_idx], 
                                 "signal_"+sys_name, cutTable["signal"]);
        HigUtilities::addBinCuts(sampleBins, baseline, systematics_vector_genmet[sys_idx].second[wgt_idx], 
                                 "signalGenMet_"+sys_name, HigUtilities::nom2genmet,
                                 cutTable["signal"]);
      }
    }
  }

  PlotMaker pm;
  set<string> boostedSignalRegion = {
                                     "BoostedEvents/boostedEvts_MC2016bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig1D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig1D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig1D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig2D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig2D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig2D_SR_noVeto.txt",
                                     };
  set<string> boostedControlRegion = {
                                     "BoostedEvents/boostedEvts_MC2016bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig1D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig1D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig1D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig2D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig2D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig2D_CR_noVeto.txt",
                                     };
  //set<string> resolvedSignalRegion = {
  //                                    "eventlist/processed_resolved_list_SR_SCAN_MC.txt",
  //                                    "eventlist/processed_resolved_list_SR_SCAN_TChiHH1D.txt",
  //                                    //"eventlist/processed_resolved_list_SR_SCAN_TChiHH2D.txt",
  //                                   };
  //set<string> resolvedControlRegion = {
  //                                    "eventlist/processed_resolved_list_CR_SCAN_MC.txt",
  //                                    "eventlist/processed_resolved_list_CR_SCAN_TChiHH1D.txt",
  //                                    //"eventlist/processed_resolved_list_CR_SCAN_TChiHH2D.txt",
  //                                    };
  //eventNumberData[0][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Boosted SR
  //eventNumberData[1][event] = {[sampleType, year, run, lumiblock, met, NSLP, LSP]} // Boosted CR
  //eventNumberData[2][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Resolved SR
  //eventNumberData[3][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Resolved CR
  vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = new vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > >();
  addEventNumberData(boostedSignalRegion, eventNumberData);
  addEventNumberData(boostedControlRegion, eventNumberData);
  //addEventNumberData(resolvedSignalRegion, eventNumberData);
  //addEventNumberData(resolvedControlRegion, eventNumberData);
  pm.SetEventVetoData(static_cast<void * >(eventNumberData));
  // Luminosity used for labeling for table
  // Luminosity used for scaling for hist1d
  bool verbose = false;
  HigUtilities::makePlots(cutTable, sampleProcesses, luminosity, pm, verbose);

  // fill mYields
  // mYields[process_tag_sampleBinLabel] = GammaParams, TableRow
  map<string, pair<GammaParams, TableRow> > mYields;
  if(unblind) HigUtilities::fillDataYields(pm, cutTable["data"], mYields, true);
  // Luminosity used for scaling
  HigUtilities::fillMcYields(pm, luminosity, cutTable["mc"], mYields, true);
  // Luminosity used for scaling
  HigUtilities::fillSignalYieldsProcesses(pm, luminosity, sampleProcesses["signal"], cutTable["signal"], mYields);
  HigUtilities::fillAverageGenMetYields(sampleProcesses["signal"], sampleBins, "signal", "signalGenMet", "signalAverageGenMet", mYields);
  for (pair<string, vector<NamedFunc>> sys : systematics_vector) {
    for (unsigned wgt_idx = 0; wgt_idx < sys.second.size(); wgt_idx++) {
      string sys_name = sys.first+to_string(wgt_idx);
      HigUtilities::fillAverageGenMetYields(sampleProcesses["signal"], sampleBins, 
          "signal_"+sys_name, "signalGenMet_"+sys_name, 
          "signalAverageGenMet_"+sys_name, mYields);
    }
  }

  // Calculate kappas from mc
  // kappas[xsigX_ysigY_metA_drmaxB]
  map<string, double> kappas;
  // kappa_uncertainties[metX_drmaxY]
  map<string, double> kappa_uncertainties;
  HigWriteDataCards::calculateKappas(mYields, dimensionBins, kappas, kappa_uncertainties);
  //// For testing
  //for (auto & item : kappas) item.second = 2;
  //for (auto & item : kappa_uncertainties) item.second = 0.001;

  //cout<<"cut name: "<<cutTable["signal"].tableRows[5].cut_.Name()<<endl;

  for (auto & process : sampleProcesses["signal"]){
    cout<<"Process name: "<<process->name_<<endl;
    string model;
    int mGluino=0, mLSP=0;
    HigUtilities::getInfoFromProcessName(process->name_, model, mGluino, mLSP);
    TString outPath = outFolder+"/datacard-"+HigUtilities::setProcessNameLong(model, mGluino, mLSP)+"_"+years_string+"_priority"+to_string(priority)+"_"+tag+".txt";
    cout<<"open "<<outPath<<endl;
    ofstream cardFile(outPath);
    HigWriteDataCards::writeDataCardHeader(sampleBins,cardFile);

    vector<vector<string> > tableValues;
    if (unblind && unblind_signalregion) HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "data", tableValues);
    else if (unblind && !unblind_signalregion) HigWriteDataCards::setDataCardObservedBlind(mYields, sampleBins, tableValues);
    else HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "mc", tableValues);
    if (do_met_average) {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalSystematics(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues, systematics_vector);
    } else {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signal", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signal", mYields, sampleBins, tableValues);      
      HigWriteDataCards::setDataCardSignalSystematics(process->name_, "signal", mYields, sampleBins, tableValues, systematics_vector);
    }
    HigWriteDataCards::setDataCardControlSystematics(controlSystematics, sampleBins, tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile);
    tableValues.clear();

    if (unblind) HigWriteDataCards::setDataCardBackground(mYields, sampleBins, "data", tableValues);
    else HigWriteDataCards::setDataCardBackground(mYields, sampleBins, "mc", tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile, true);
    tableValues.clear();
    HigWriteDataCards::setDataCardKappa(kappas, kappa_uncertainties, dimensionBins, tableValues);
    HigWriteDataCards::writeTableValues(tableValues, cardFile);
    tableValues.clear();
    //writeDataCardClosure()
    //writeDataCardUncertainties()
    //writeDataCardParam()
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

namespace HigWriteDataCards{
  void calculateKappas(std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins, std::map<std::string, double> & kappas, std::map<std::string, double> & kappa_uncertainties) {
    map<string, vector<pair<string, string> > > t_dimensionBins = dimensionBins;
    // Combine dimensions to one list of dimensions.
    HigUtilities::combineDimensionBins(t_dimensionBins);

    auto dimension = t_dimensionBins.begin();
    for (unsigned iBin =0; iBin < dimension->second.size(); ++iBin) {
      for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
        GammaParams xsig_ysig = mYields.at("mc_xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first).first;
        GammaParams xsig_ybkg = mYields.at("mc_xsig"+to_string(iXbin)+"_ybkg_"+(dimension->second)[iBin].first).first;
        GammaParams xbkg_ysig = mYields.at("mc_xbkg_ysig_"+(dimension->second)[iBin].first).first;
        GammaParams xbkg_ybkg = mYields.at("mc_xbkg_ybkg_"+(dimension->second)[iBin].first).first;
        vector<vector<float> > entries = {{static_cast<float>(xbkg_ybkg.NEffective())}, {static_cast<float>(xsig_ybkg.NEffective())}, {static_cast<float>(xbkg_ysig.NEffective())}, {static_cast<float>(xsig_ysig.NEffective())}};
        vector<vector<float> > entry_weights = {{static_cast<float>(xbkg_ybkg.Weight())}, {static_cast<float>(xsig_ybkg.Weight())}, {static_cast<float>(xbkg_ysig.Weight())}, {static_cast<float>(xsig_ysig.Weight())}};
        vector<float> powers = {1, -1, -1, 1};
        float kappa_up = 0, kappa_down = 0;
        double kappa = calcKappa(entries, entry_weights, powers, kappa_down, kappa_up);
        double kappa_uncertainty = kappa_up>kappa_down ? kappa_up:kappa_down;
        //cout<<xsig_ysig.Yield()<<" "<<xsig_ybkg.Yield()<<" "<<xbkg_ysig.Yield()<<" "<<xbkg_ybkg.Yield()<<endl;
        //cout<<"  kappa: "<<kappa<<" stat. uncertainty: "<<kappa_uncertainty<<endl;
        string kappaName = "xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first;
        kappas[kappaName] = kappa;
        kappa_uncertainties[kappaName] = kappa_uncertainty;
      }
    }
  }

  void setControlSystematics(map<string, map<string, float> > & controlSystematics) {
    controlSystematics["xsig0_ybkg_met0_drmax0"]["ttbar"] = 1.12; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["ttbar"] = 1.17;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["ttbar"] = 1.02;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["ttbar"] = 1.09;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["ttbar"] = 1.10;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["ttbar"] = 1.16;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["ttbar"] = 1.02;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["ttbar"] = 1.09;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["ttbar"] = 1.06;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["ttbar"] = 1.12;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["ttbar"] = 1.01;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["ttbar"] = 1.08;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["ttbar"] = 1.06;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["ttbar"] = 1.08;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["ttbar"] = 1.01;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["ttbar"] = 1.08;

    controlSystematics["xsig0_ybkg_met0_drmax0"]["vjets"] = 1.01; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["vjets"] = 1.01;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["vjets"] = 1.04;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["vjets"] = 1.05;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["vjets"] = 1.10;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["vjets"] = 1.06;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["vjets"] = 1.06;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["vjets"] = 1.10;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["vjets"] = 1.02;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["vjets"] = 1.02;

    controlSystematics["xsig0_ybkg_met0_drmax0"]["qcd"] = 1.00; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["qcd"] = 1.02;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["qcd"] = 1.01;
  }

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

  void make_npv_plots(string tag, map<string, HigUtilities::HistInformation> & hist_info) {
    //add npv plot with no cuts past higloose+filters
    PlotOpt plot_opt("txt/plot_styles.txt", "CMSPaper");
    plot_opt = plot_opt.Title(PlotOptTypes::TitleType::info)
      .Bottom(PlotOptTypes::BottomType::off)
      .Overflow(PlotOptTypes::OverflowType::none)
      .YAxis(PlotOptTypes::YAxisType::linear)
      .Stack(PlotOptTypes::StackType::lumi_shapes).LegendColumns(3);
    //if (tag=="mc") plot_opt = plot_opt.Stack(StackType::signal_overlay);
    HigUtilities::HistInformation this_hist_info;
    this_hist_info.axis_ = new Axis(100, -0.5, 99.5, "npv", "N_{pv}", {});
    this_hist_info.cut_ = new NamedFunc("1");
    //this_hist_info.weight_ = new NamedFunc("weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years*Higfuncs::w_pileup_nosignal);
    this_hist_info.weight_ = new NamedFunc(Higfuncs::final_weight);
    this_hist_info.plot_opt_ = new PlotOpt(plot_opt);
    this_hist_info.figure_index = -1;
    hist_info[tag] = this_hist_info;
  }

  void fill_npv_hists(PlotMaker& pm, map<string, HigUtilities::HistInformation> & hist_info, map<string, vector<shared_ptr<Process>>> & sample_processes, map<string, TH1D> &h_mc_npv, TH1D **h_data_npv) {
    Hist1D * sig_h1d = static_cast<Hist1D*>(pm.Figures()[hist_info["signal"].figure_index].get());
    for (auto & process : sample_processes["signal"]) {
      Hist1D::SingleHist1D* sig_sh1d = static_cast<Hist1D::SingleHist1D*>(sig_h1d->GetComponent(process.get()));
      h_mc_npv[process->name_] = sig_sh1d->scaled_hist_;
      //h_mc_npv[process->name_].Scale(1.0/h_mc_npv[process->name_].Integral());
    }

    TH1D* h_ret_npv = new TH1D("h_ret_npv","(Pseudo)data N_{pv}", 100, -0.5, 99.5);

    string signal_name = "mc";
    if (unblind) {
      signal_name = "data";
    }
    Hist1D * data_h1d = static_cast<Hist1D*>(pm.Figures()[hist_info[signal_name].figure_index].get());
    for (auto & process : sample_processes[signal_name]) {
      Hist1D::SingleHist1D* data_sh1d = static_cast<Hist1D::SingleHist1D*>(data_h1d->GetComponent(process.get()));
      h_ret_npv->Add(&(data_sh1d->scaled_hist_));
    }
    //h_ret_npv->Scale(1.0/h_ret_npv->Integral());
    (*h_data_npv) = h_ret_npv;
  }

  void delete_hist_info(map<string, HigUtilities::HistInformation> & hist_info) {
    for (auto single_hist_info : hist_info) {
      string const & type = single_hist_info.first;
      delete hist_info[type].axis_;
      delete hist_info[type].cut_;
      delete hist_info[type].weight_;
      delete hist_info[type].plot_opt_;
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
      //row[iBin*2+3] = to_string(int(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield())); // Integerize. Round down
      //row[iBin*2+3] = to_string(int(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield())+1); // Integerize. Round up
      //row[iBin*2+3] = RoundNumber(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield(), 0, 1).Data(); // Integerize. Round
      row[iBin*2+3] = to_string(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield());
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }

  void setDataCardObservedBlind(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues){
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
      //row[iBin*2+3] = to_string(int(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield())); // Integerize. Round down
      //row[iBin*2+3] = to_string(int(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield())+1); // Integerize. Round up
      //row[iBin*2+3] = RoundNumber(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield(), 0, 1).Data(); // Integerize. Round
      row[iBin*2+3] = to_string(mYields.at("data_"+sampleBins[iBin].first).first.Yield());
      if (iBin%6==4 || iBin%6==5) row[iBin*2+3] = to_string(mYields.at("mc_"+sampleBins[iBin].first).first.Yield());
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

  void setDataCardSignalSystematics(string const & processName, string const & signalAverageGenMetTag, 
      map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, 
      vector<vector<string> > & tableValues, vector<pair<string, vector<NamedFunc>>> systematics_vector)
  {

    vector<string> row(2+2*sampleBins.size());

    for (pair<string, vector<NamedFunc>> sys : systematics_vector) {
      setRow(row,"-");
      row[0] = sys.first;
      row[1] = "lnN";
      for (unsigned int bin_idx = 0; bin_idx < (sampleBins.size()); bin_idx++) {
        if (sys.first == "LumiSyst") {
          row[2+2*bin_idx] = "1.0243";
        }
        else if (sys.first == "SignalJetID") {
          row[2+2*bin_idx] = "1.01";
        }
        else if (sys.first == "SignalMETFastSIM") {
          //reco MET vs gen MET systematic
          string reclabel = processName + "_signal_" + sampleBins[bin_idx].first;
          string genlabel = processName + "_signalGenMet_" + sampleBins[bin_idx].first;
          string avrlabel = processName + "_signalAverageGenMet_" + sampleBins[bin_idx].first;
          float recmet_value = mYields.at(reclabel).first.Yield();
          float avrmet_value = mYields.at(avrlabel).first.Yield();
          float genmet_value = mYields.at(genlabel).first.Yield();
          float syst_value = 2.0;
          if (avrmet_value != 0)
            syst_value = (recmet_value-genmet_value)/2.0/avrmet_value;
          //if (fabs(syst_value) > 0.5 && syst_value != 2.0) 
          //  std::cout << "DEBUG: large fs met systematic. RECO value: "
          //       << recmet_value << ", GEN value: " << genmet_value << std::endl;
          if (syst_value > 1.0) syst_value = 1.0;
          if (syst_value <= -1.0) syst_value = -0.99;
          //convention for sign
          row[2+2*bin_idx] = to_string(syst_value+1.0);
        }
        else if (sys.first == "SignalScale") {
          //renormalization and factorization scale systematic
          string label = processName + "_" + signalAverageGenMetTag + "_" + sampleBins[bin_idx].first;
          float nominal_value = mYields.at(label).first.Yield();
          if (nominal_value == 0) {
            //no signal, poor stats
            row[2+2*bin_idx] = "2.00/2.00";
          }
          else {
            float max_variation_up = 0.0;
            float max_variation_down = 0.0;
            for (unsigned wgt_idx = 0; wgt_idx < sys.second.size(); wgt_idx++) {
              label = processName + "_" + signalAverageGenMetTag + "_" + sys.first + to_string(wgt_idx) + "_" + sampleBins[bin_idx].first;
              float variation_value = mYields.at(label).first.Yield();
              float variation = variation_value/nominal_value-1.0;
              if (variation > max_variation_up) max_variation_up = variation;
              if (variation < max_variation_down) max_variation_down = variation;
            }
            if (max_variation_up > 1.0) max_variation_up = 1.0;
            if (max_variation_down <= -1.0) max_variation_down = -0.99;
            row[2+2*bin_idx] = to_string(max_variation_down+1.0)+"/"+to_string(max_variation_up+1.0);
          }
        }
        else {
          //normal up/down systematic
          string label = processName + "_" + signalAverageGenMetTag + "_" + sampleBins[bin_idx].first;
          float nominal_value = mYields.at(label).first.Yield();
          if (nominal_value == 0 || std::isnan(nominal_value) || std::isinf(nominal_value)) {
            //no signal, poor stats - make all variations up??
            row[2+2*bin_idx] = "2.00/2.00";
          }
          else {
            label = processName + "_" + signalAverageGenMetTag + "_" + sys.first + "0_" + sampleBins[bin_idx].first;
            float value_up = mYields.at(label).first.Yield();
            label = processName + "_" + signalAverageGenMetTag + "_" + sys.first + "1_" + sampleBins[bin_idx].first;
            float value_down = mYields.at(label).first.Yield();
            float syst_up = (value_up-nominal_value)/nominal_value;
            float syst_down = (value_down-nominal_value)/nominal_value;
            if (syst_up > 1.0 || std::isnan(syst_up) || std::isinf(syst_up)) syst_up = 1.0;
            if (syst_up <= -1.0) syst_up = -0.99;
            if (syst_down > 1.0 || std::isnan(syst_down) || std::isinf(syst_down)) syst_down = 1.0;
            if (syst_down <= -1.0) syst_down = -0.99;
            row[2+2*bin_idx] = to_string(syst_down+1.0)+"/"+to_string(syst_up+1.0);
          }
        } //sys name
      } //bin loop
      tableValues.push_back(row);
    } //systematic loop

    setRow(row,"");
    tableValues.push_back(row);
  }

  void setDataCardControlSystematics(map<string, map<string, float> > controlSystematics, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues) 
  {
    // controlSystematics[xsig0_ybkg_met0_drmax0]["ttbar"] = 1.03

    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    
    // controlSystematicNames[0] = ttbar
    vector<string> controlSystematicNames;
    for (auto it = controlSystematics.begin()->second.begin(); it != controlSystematics.begin()->second.end(); ++it) controlSystematicNames.push_back(it->first);

    // Set control sample systematics
    for (unsigned iControlSystematic = 0; iControlSystematic < controlSystematicNames.size(); ++iControlSystematic) {
      string systematicName = controlSystematicNames[iControlSystematic]; // ex: ttbar
      setRow(row,"-");
      row[0] = "cr_"+systematicName;
      row[1] = "lnN";
      for (unsigned rBin=0;rBin<sampleBins.size();++rBin) {
        string regionName = sampleBins[rBin].first;
        if (controlSystematics.find(regionName) != controlSystematics.end()) {
          row[2+2*rBin+1] = to_string(controlSystematics[regionName][systematicName]);
        }
      }
      tableValues.push_back(row);
    }

    setRow(row, "");
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
        float mcYield = mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield();
        //row[4] = to_string(int(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield())); //Integerize. Round down
        //row[4] = to_string(int(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield())+1); //Integerize. Round up
        //row[4] = RoundNumber(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield(), 0, 1).Data(); //Integerize. Round
        row[4] = to_string(mcYield);
        // Infinity defined in https://root.cern.ch/doc/master/RooNumber_8cxx_source.html
        row[5] = "[0,"+to_string(mcYield*10)+"]";
      }
      else 
      {
        //row[4] = "(@0*@1/@2)";
        //vector<string> xydimension;
        //HigUtilities::stringToVectorString(sampleBins[rBin].first, xydimension, "_");
        //vector<string> aNameVector = xydimension;
        //aNameVector[0] = "xbkg";
        //aNameVector[1] = "ybkg";
        //string aName;
        //HigUtilities::vectorStringToString(aNameVector, aName, "_");
        //vector<string> bNameVector = xydimension;
        //bNameVector[1] = "ybkg";
        //string bName;
        //HigUtilities::vectorStringToString(bNameVector, bName, "_");
        //vector<string> cNameVector = xydimension;
        //cNameVector[0] = "xbkg";
        //string cName;
        //HigUtilities::vectorStringToString(cNameVector, cName, "_");
        ////cout<<aName<<" "<<bName<<" "<<cName<<endl;
        //row[5] = "rp_"+cName+",rp_"+bName+",rp_"+aName;

        row[4] = "(@0*@1/@2)*@3";
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
        string kappaName;
        vector<string> kappaNameVector = xydimension;
        HigUtilities::vectorStringToString(kappaNameVector, kappaName, "_");
        //cout<<aName<<" "<<bName<<" "<<cName<<endl;
        row[5] = "rp_"+cName+",rp_"+bName+",rp_"+aName+",kappa_"+kappaName;
      }
      tableValues.push_back(row);
      setRow(row,"");
    }
  
  }

  void setDataCardKappa(map<string, double > & kappas, map<string, double > & kappa_uncertainties, map<string, vector<pair<string, string> > > & dimensionBins, vector<vector<string> > & tableValues)
  {
    map<string, vector<pair<string, string> > > t_dimensionBins = dimensionBins;
    // Combine dimensions to one list of dimensions.
    HigUtilities::combineDimensionBins(t_dimensionBins);

    auto dimension = t_dimensionBins.begin();
    vector<string> row(4);
    for (unsigned iBin =0; iBin < dimension->second.size(); ++iBin) {
      for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
        //string kappaName = "xsig"+to_string(iXbin)+"_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax);
        string kappaName = "xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first;
        row[0] = "kappa_"+kappaName;
        row[1] = "param";
        row[2] = to_string(kappas[kappaName]);
        row[3] = to_string(kappa_uncertainties[kappaName]);
        tableValues.push_back(row);
        setRow(row,"");
      }
    }
    cout<<"End"<<endl;

    //for (unsigned iMet = 0; iMet < dimensionBins["met"].size(); iMet++) {
    //  for (unsigned iDrmax = 0 ; iDrmax < dimensionBins["drmax"].size(); iDrmax++) {
    //    for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
    //      string kappaName = "xsig"+to_string(iXbin)+"_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax);
    //      row[0] = "kappa_"+kappaName;
    //      row[1] = "param";
    //      row[2] = to_string(kappas[kappaName]);
    //      row[3] = to_string(kappa_uncertainties[kappaName]);
    //      tableValues.push_back(row);
    //      setRow(row,"");
    //    }
    //  }
    //}
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
    while(true){
      static struct option long_options[] = {
        {"output_folder", required_argument, 0, 'o'},
        {"mass_points", required_argument, 0, 'p'},
        {"years", required_argument, 0, 'y'},
        {"luminosity", required_argument, 0, 'l'},
        {"dimension", required_argument, 0, 'd'},
        {"tag", required_argument, 0, 't'},
        {"unblind", no_argument, 0, 'u'},
        {"unblind_signalregion", no_argument, 0, 0},
        {"recomet", no_argument, 0, 0},
        {"higgsino_model", required_argument, 0, 'm'},
        {"priority", required_argument, 0, 'r'},
        {"t5hhrange", required_argument, 0, 's'},
        {0, 0, 0, 0}
      };
  
      char opt = -1;
      int option_index;
      opt = getopt_long(argc, argv, "o:p:g:y:l:d:t:m:r:s:nu12", long_options, &option_index);
      if( opt == -1) break;

      string optname;
      switch(opt){
        case 'o': outFolder = optarg; break;
        case 'p': mass_points_string = optarg; break;
        case 'g': mass_point_glob = optarg; break;
        case 'y': years_string = optarg; break;
        case 'l': luminosity = atof(optarg); break;
        case 't': tag = optarg; break;
        case 'm': higgsino_model = optarg; break;
        case 'r': priority = atoi(optarg); break;
        case 's': t5hh_range = atoi(optarg); break;
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
        case '1':
          is1D = true;
          break;
        case '2':
          is1D = false;
          break;
        case 0:
          optname = long_options[option_index].name;
          if(optname == "recomet"){
            do_met_average = false;
          }else if(optname == "unblind_signalregion"){
            unblind_signalregion = true;
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
