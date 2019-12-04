#include "core/functions.hpp"

#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace Functions{

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
    float bb_wp = 0.3;
    if (b.met()<=300) bb_wp = 0.6;
    // 0 = D, 1 = B1, 2 = B2, 3 = C, 4 = A1, 5 = A2
    bool j1mass = b.fjet_msoftdrop()->at(0)>95 && b.fjet_msoftdrop()->at(0)<=145;
    bool j2mass = b.fjet_msoftdrop()->at(1)>95 && b.fjet_msoftdrop()->at(1)<=145;
    bool j1bb = b.fjet_mva_hbb_btv()->at(0)>bb_wp;
    bool j2bb = b.fjet_mva_hbb_btv()->at(1)>bb_wp;

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

  const NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
    int tmp_ntrub(0);
    for (unsigned i(0); i<b.jet_pt()->size(); i++){
      if (!b.jet_h1d()->at(i) && !b.jet_h2d()->at(i)) continue;
      if (b.jet_hflavor()->at(i)==5) tmp_ntrub++;
    }
    return tmp_ntrub;
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
}
