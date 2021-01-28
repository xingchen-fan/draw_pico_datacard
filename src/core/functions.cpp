#include "core/functions.hpp"
#include <regex>

#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace Functions{

  const NamedFunc hem_veto("hem_veto",[](const Baby &b) -> NamedFunc::ScalarType{
    if(abs(b.SampleType()) == 2018) {
      if (b.SampleType()<0) {
        if (b.run() >= 319077) { 
          if(b.nel() > 0) {
            for(size_t i = 0; i < b.lep_pt()->size(); i++) {
              if (abs(b.lep_pdgid()->at(i))==13) continue;
              if(b.lep_eta()->at(i) < -1.5 && (b.lep_phi()->at(i) > -1.6 && b.lep_phi()->at(i) < -0.8)) 
                return static_cast<float>(0);
            }
          }
          for(size_t i = 0; i < b.jet_pt()->size(); i++) {
            if(Functions::IsGoodJet(b,i) && b.jet_eta()->at(i) < -1.5 && (b.jet_phi()->at(i) > -1.6 && b.jet_phi()->at(i) < -0.8)) 
              return static_cast<float>(0);
          }
        }
      } else {
        if ((b.event()%1961) < 1296) { 
          if(b.nel() > 0) {
            for(size_t i = 0; i < b.lep_pt()->size(); i++) {
              if (abs(b.lep_pdgid()->at(i))==13) continue;
              if(b.lep_eta()->at(i) < -1.5 && (b.lep_phi()->at(i) > -1.6 && b.lep_phi()->at(i) < -0.8)) 
                return static_cast<float>(0);
            }
          }
          for(size_t i = 0; i < b.jet_pt()->size(); i++) {
            if(Functions::IsGoodJet(b,i) && b.jet_eta()->at(i) < -1.5 && (b.jet_phi()->at(i) > -1.6 && b.jet_phi()->at(i) < -0.8)) 
              return static_cast<float>(0);
          }
        }
      }
    }
    return static_cast<float>(1);
  });
  
  const NamedFunc lowDphiFix("lowDphiFix",[](const Baby &b) -> NamedFunc::ScalarType{
    bool lowDphi(false);
    for (unsigned i(0); i<b.jet_pt()->size(); i++){
      if (i>3) break;
      if (i<=1 && b.jet_met_dphi()->at(i)<=0.5) lowDphi = true;
      else if (i>=2 && b.jet_met_dphi()->at(i)<=0.3) lowDphi = true;
    }
    return lowDphi;
  });

  const NamedFunc boostedRegionIdx("boostedRegionIdx",[](const Baby &b) -> NamedFunc::ScalarType{
    float bb_wp = 0.7;
    if (b.met()<=300) bb_wp = 0.86;
    // 0 = D, 1 = B1, 2 = B2, 3 = C, 4 = A1, 5 = A2
    bool j1mass = b.fjet_msoftdrop()->at(0)>95 && b.fjet_msoftdrop()->at(0)<=145;
    bool j2mass = b.fjet_msoftdrop()->at(1)>95 && b.fjet_msoftdrop()->at(1)<=145;
    bool j1bb = b.fjet_deep_md_hbb_btv()->at(0)>bb_wp;
    bool j2bb = b.fjet_deep_md_hbb_btv()->at(1)>bb_wp;

    if (j1mass && j2mass) {
      if (j1bb && j2bb) return 5; 
      else if (j1bb || j2bb) return 4; 
      else return 3;
    } else { 
      if (j1bb && j2bb) return 2;
      else if (j1bb || j2bb) return 1;
      else return 0;
    }
  });

  const NamedFunc amBoostedRegionIdx("amBoostedRegionIdx",[](const Baby &b) -> NamedFunc::ScalarType{
    float bb_wp = 0.7;
    if (b.met()<=300) bb_wp = 0.86;
    // 0 = D, 1 = B1, 2 = B2, 3 = C, 4 = A1, 5 = A2
    float am = (b.fjet_msoftdrop()->at(0)+b.fjet_msoftdrop()->at(1))/2;
    bool j1bb = b.fjet_deep_md_hbb_btv()->at(0)>bb_wp;
    bool j2bb = b.fjet_deep_md_hbb_btv()->at(1)>bb_wp;

    if (am>95 && am<=145) {
      if (j1bb && j2bb) return 5; 
      else if (j1bb || j2bb) return 4; 
      else return 3;
    } else { 
      if (j1bb && j2bb) return 2;
      else if (j1bb || j2bb) return 1;
      else return 0;
    }
  });

  const NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
    int tmp_ntrub(0);
    for (unsigned i(0); i<b.jet_pt()->size(); i++){
      if (!b.jet_h1d()->at(i) && !b.jet_h2d()->at(i)) continue;
      if (b.jet_hflavor()->at(i)==5) tmp_ntrub++;
    }
    return tmp_ntrub;
  });

  const NamedFunc w_pileup("w_pileup",[](const Baby &b) -> NamedFunc::ScalarType{
    float weight = 1.;

    if(b.type()>0 && b.type()<1000) weight = 1; // data
    else { // mc
      vector<double> weights;
      // weights are from get_puweights.py
      if (b.SampleType()==2016) {
        weights = vector<double>({3.661e-01, 8.939e-01, 1.198e+00, 9.627e-01, 1.121e+00, 1.165e+00, 7.956e-01, 4.958e-01, 7.422e-01, 8.789e-01, 9.642e-01, 1.072e+00, 1.125e+00, 1.176e+00, 1.202e+00, 1.208e+00, 1.200e+00, 1.183e+00, 1.144e+00, 1.097e+00, 1.066e+00, 1.051e+00, 1.052e+00, 1.051e+00, 1.050e+00, 1.058e+00, 1.072e+00, 1.083e+00, 1.096e+00, 1.108e+00, 1.095e+00, 1.083e+00, 1.041e+00, 9.858e-01, 9.108e-01, 8.209e-01, 7.168e-01, 6.100e-01, 5.031e-01, 4.048e-01, 3.092e-01, 2.279e-01, 1.637e-01, 1.132e-01, 7.730e-02, 5.092e-02, 3.189e-02, 2.009e-02, 1.226e-02, 7.426e-03, 4.380e-03, 2.608e-03, 1.566e-03, 9.714e-04, 7.292e-04, 6.727e-04, 7.305e-04, 9.488e-04, 1.355e-03, 1.894e-03, 3.082e-03, 4.097e-03, 4.874e-03, 5.256e-03, 5.785e-03, 5.515e-03, 5.000e-03, 4.410e-03, 4.012e-03, 3.548e-03, 3.108e-03, 2.702e-03, 2.337e-03, 2.025e-03, 1.723e-03});
      } else if (b.SampleType()==2017) {
        weights =  vector<double>({1.835e-01, 3.933e+00, 3.471e+00, 2.492e+00, 1.625e+00, 1.515e+00, 1.288e+00, 1.278e+00, 6.158e-01, 1.452e+00, 1.498e+00, 1.487e+00, 1.331e+00, 1.164e+00, 1.079e+00, 1.054e+00, 1.081e+00, 1.129e+00, 1.166e+00, 1.189e+00, 1.213e+00, 1.238e+00, 1.260e+00, 1.271e+00, 1.273e+00, 1.271e+00, 1.271e+00, 1.267e+00, 1.274e+00, 1.252e+00, 1.221e+00, 1.170e+00, 1.109e+00, 1.038e+00, 9.694e-01, 9.119e-01, 8.668e-01, 8.345e-01, 7.884e-01, 7.508e-01, 7.592e-01, 7.932e-01, 8.583e-01, 9.588e-01, 1.094e+00, 1.256e+00, 1.420e+00, 1.496e+00, 1.531e+00, 1.462e+00, 1.337e+00, 1.154e+00, 9.506e-01, 7.497e-01, 5.698e-01, 4.107e-01, 2.900e-01, 1.987e-01, 1.376e-01, 9.664e-02, 6.923e-02, 5.086e-02, 3.844e-02, 2.996e-02, 2.409e-02, 1.712e-02, 1.248e-02, 1.077e-02, 9.597e-03, 8.812e-03, 8.310e-03, 8.018e-03, 7.885e-03, 7.875e-03, 6.337e-03, 5.334e-03, 5.444e-03, 5.584e-03, 5.744e-03, 5.916e-03, 6.090e-03, 6.257e-03, 6.408e-03, 5.095e-03, 4.228e-03, 4.259e-03, 4.266e-03, 4.247e-03, 4.203e-03, 4.137e-03, 4.052e-03, 3.950e-03, 2.943e-03, 2.296e-03, 2.200e-03, 2.101e-03, 2.002e-03, 1.902e-03, 1.804e-03});
      } else if (b.SampleType()==2018) {
        weights = vector<double>({9.687e+03, 1.327e+01, 4.377e+01, 1.869e+01, 1.251e+01, 9.042e+00, 6.579e+00, 4.873e+00, 3.629e+00, 2.763e+00, 2.224e+00, 1.900e+00, 1.702e+00, 1.579e+00, 1.504e+00, 1.465e+00, 1.450e+00, 1.451e+00, 1.459e+00, 1.464e+00, 1.458e+00, 1.437e+00, 1.401e+00, 1.353e+00, 1.301e+00, 1.250e+00, 1.205e+00, 1.167e+00, 1.137e+00, 1.115e+00, 1.100e+00, 1.089e+00, 1.082e+00, 1.077e+00, 1.073e+00, 1.069e+00, 1.063e+00, 1.052e+00, 1.036e+00, 1.012e+00, 9.794e-01, 9.373e-01, 8.862e-01, 8.272e-01, 7.618e-01, 6.923e-01, 6.211e-01, 5.505e-01, 4.827e-01, 4.194e-01, 3.618e-01, 3.105e-01, 2.659e-01, 2.276e-01, 1.952e-01, 1.679e-01, 1.450e-01, 1.259e-01, 1.097e-01, 9.590e-02, 8.401e-02, 7.361e-02, 6.442e-02, 5.623e-02, 4.889e-02, 4.231e-02, 3.643e-02, 3.122e-02, 2.662e-02, 2.262e-02, 1.915e-02, 1.618e-02, 1.363e-02, 1.147e-02, 9.633e-03, 8.075e-03, 6.750e-03, 5.623e-03, 4.661e-03, 3.838e-03, 3.134e-03, 2.531e-03, 2.016e-03, 1.578e-03, 1.208e-03, 8.991e-04, 6.465e-04, 4.458e-04, 2.929e-04, 1.824e-04, 1.094e-04, 6.258e-05, 3.428e-05, 1.807e-05, 9.221e-06, 4.578e-06, 2.221e-06, 1.058e-06, 4.964e-07, 2.308e-07});
      }
      unsigned npu_tru_mean = static_cast<unsigned> (b.npu_tru_mean());
      if (npu_tru_mean>=weights.size()) npu_tru_mean = weights.size()-1;

      weight = weights.at(npu_tru_mean);
    }

    return weight;
  });

  const NamedFunc SampleType("SampleType",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.SampleType();
  });


  bool IsGoodJet(const Baby &b, size_t ijet){
    return ijet<b.jet_pt()->size()
      && b.jet_pt()->at(ijet) > 30.
      && fabs(b.jet_eta()->at(ijet))<2.4
      && !b.jet_islep()->at(ijet);
  }

  NamedFunc MismeasurementCorrection(const string &file_path,
                                     const string &mismeas_scenario,
                                     Variation variation){
    ConfigParser cp;
    cp.Load(file_path, mismeas_scenario);

    string name = "w_mismeas";
    NamedFunc reweight_cut = cp.GetOpt<std::string>("reweight_cut");
    double wgt = 1.;
    switch(variation){
    case Variation::central:
      wgt = cp.GetOpt<double>("w_central");
      break;
    case Variation::up:
      wgt = cp.GetOpt<double>("w_up");
      break;
    case Variation::down:
      wgt = cp.GetOpt<double>("w_down");
      break;
    default:
      wgt = 1.;
      ERROR("Invalid variation: "+to_string(static_cast<int>(variation)));
      break;
    }

    if(reweight_cut.IsScalar()){
      return NamedFunc(name, [reweight_cut, wgt](const Baby &b){
          return reweight_cut.GetScalar(b) ? wgt : 1.;
        });
    }else{
      return NamedFunc(name, [reweight_cut, wgt](const Baby &b){
          NamedFunc::VectorType cuts = reweight_cut.GetVector(b);
          NamedFunc::VectorType wgts(cuts.size());
          for(size_t i = 0; i < cuts.size(); ++i){
            wgts.at(i) = cuts.at(i) ? wgt : 1.;
          }
          return wgts;
        });
    }
  }

  NamedFunc MismeasurementWeight(const std::string &file_path,
                                 const std::string &mismeas_scenario){
    ConfigParser cp;
    cp.Load(file_path, mismeas_scenario);
    NamedFunc cut = cp.GetOpt<std::string>("mismeas_cut");
    NamedFunc wgt = cp.GetOpt<std::string>("mismeas_wgt");
    string name = "w_sys_mm";
    if(cut.IsScalar()){
      if(wgt.IsScalar()){
        return NamedFunc(name, [cut, wgt](const Baby &b){
            return cut.GetScalar(b) ? wgt.GetScalar(b) : 1.;
          });
      }else{
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType wgts = wgt.GetVector(b);
            return cut.GetScalar(b) ? wgts : NamedFunc::VectorType(wgts.size(), 1.);
          });
      }
    }else{
      if(wgt.IsScalar()){
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType cuts = cut.GetVector(b);
            NamedFunc::VectorType wgts(cuts.size());
            NamedFunc::ScalarType scalar_wgt = wgt.GetScalar(b);
            for(size_t i = 0; i < wgts.size(); ++i){
              wgts.at(i) = cuts.at(i) ? scalar_wgt : 1.;
            }
            return wgts;
          });
      }else{
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType cuts = cut.GetVector(b);
            NamedFunc::VectorType wgts = wgt.GetVector(b);
            size_t min_size = min(cuts.size(), wgts.size());
            NamedFunc::VectorType out(min_size);
            for(size_t i = 0; i < min_size; ++i){
              out.at(i) = cuts.at(i) ? wgts.at(i) : 1.;
            }
            return out;
          });
      }
    }
  }

  string regSearch(string inString, string const & regExPattern) {
      std::smatch match;
      std::regex pattern(regExPattern);
      std::regex_search (inString, match, pattern);
      return match.str();
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
  bool findStringIC(const std::string & strHaystack, const std::string & strNeedle)
  {
    auto it = std::search(
      strHaystack.begin(), strHaystack.end(),
      strNeedle.begin(),   strNeedle.end(),
      [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); }
    );
    return (it != strHaystack.end() );
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
          // Need to reject if 1D is compared with 2D
          int babyDimension = (*b.FileNames().begin()).find("HToBB_2D") == string::npos ? 1 : 2;
          int eventListDimension = get<0>(it).find("2D") == string::npos ? 1 : 2;
          if (babyDimension != eventListDimension) continue;
          //if ((*b.FileNames().begin()).find("HToBB_2D") != get<0>(it).find("2D")) continue;
        } else if (!findStringIC(*b.FileNames().begin(), get<0>(it))) continue;
        //if (isSignal) cout<<"event: "<<b.event()<<" => sampleType: "<<get<0>(it)<<" year: "<<get<1>(it)<<" run: "<<get<2>(it)<<" lumiblock: "<<get<3>(it)<<" met: "<<get<4>(it)<<" nlsp_mass: "<<get<5>(it)<<" lsp_mass: "<<get<6>(it)<<endl;
  
        //if (b.met()<300) cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
        //cout<<"event: "<<b.event()<<" => sampleType: "<<get<0>(it)<<" year: "<<get<1>(it)<<" run: "<<get<2>(it)<<" lumiblock: "<<get<3>(it)<<" met: "<<get<4>(it)<<" nlsp_mass: "<<get<5>(it)<<" lsp_mass: "<<get<6>(it)<<endl;
  
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

}
