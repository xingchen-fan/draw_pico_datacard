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

  bool IsGoodJet(const Baby &b, size_t ijet){
    return ijet<b.jet_pt()->size()
      && b.jet_pt()->at(ijet) > 30.
      && fabs(b.jet_eta()->at(ijet))<2.4
      && !b.jet_islep()->at(ijet);
  }

//// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
const NamedFunc eff_higtrig("eff_higtrig", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup, errdown; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup=0;errdown=0;
  errup+=errdown;
  if(b.type()>0 && b.type()<1000) eff = 1;

  //// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
  //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31])";
  else if(b.nvlep()==0){
    if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
      if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.242; errup = 0.030; errdown = 0.028;}
      else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.320; errup = 0.007; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.359; errup = 0.005; errdown = 0.005;}
      else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {eff = 0.379; errup = 0.002; errdown = 0.002;}
      else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {eff = 0.306; errup = 0.002; errdown = 0.002;}
      else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.217; errup = 0.034; errdown = 0.031;}
      else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.344; errup = 0.007; errdown = 0.007;}
      else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.394; errup = 0.005; errdown = 0.005;}
      else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {eff = 0.421; errup = 0.003; errdown = 0.003;}
      else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {eff = 0.331; errup = 0.002; errdown = 0.002;}
      else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.211; errup = 0.036; errdown = 0.033;}
      else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.387; errup = 0.008; errdown = 0.008;}
      else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.434; errup = 0.006; errdown = 0.006;}
      else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {eff = 0.464; errup = 0.003; errdown = 0.003;}
      else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {eff = 0.363; errup = 0.002; errdown = 0.002;}
      else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.263; errup = 0.039; errdown = 0.035;}
      else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.406; errup = 0.009; errdown = 0.009;}
      else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.469; errup = 0.006; errdown = 0.006;}
      else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {eff = 0.503; errup = 0.003; errdown = 0.003;}
      else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {eff = 0.397; errup = 0.002; errdown = 0.002;}
      else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.328; errup = 0.045; errdown = 0.042;}
      else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.434; errup = 0.010; errdown = 0.010;}
      else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.507; errup = 0.007; errdown = 0.007;}
      else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {eff = 0.545; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {eff = 0.422; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.256; errup = 0.047; errdown = 0.042;}
      else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.463; errup = 0.011; errdown = 0.011;}
      else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.520; errup = 0.007; errdown = 0.007;}
      else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {eff = 0.582; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {eff = 0.455; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.234; errup = 0.048; errdown = 0.043;}
      else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.497; errup = 0.011; errdown = 0.011;}
      else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.558; errup = 0.008; errdown = 0.008;}
      else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {eff = 0.609; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {eff = 0.478; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.330; errup = 0.055; errdown = 0.051;}
      else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.479; errup = 0.012; errdown = 0.012;}
      else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.576; errup = 0.008; errdown = 0.008;}
      else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {eff = 0.646; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {eff = 0.504; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.406; errup = 0.053; errdown = 0.051;}
      else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.537; errup = 0.013; errdown = 0.013;}
      else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.592; errup = 0.009; errdown = 0.009;}
      else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {eff = 0.682; errup = 0.005; errdown = 0.005;}
      else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {eff = 0.528; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.530; errup = 0.068; errdown = 0.069;}
      else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.557; errup = 0.014; errdown = 0.014;}
      else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.609; errup = 0.009; errdown = 0.009;}
      else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {eff = 0.705; errup = 0.005; errdown = 0.005;}
      else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {eff = 0.557; errup = 0.004; errdown = 0.004;}
      else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.407; errup = 0.045; errdown = 0.043;}
      else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.549; errup = 0.011; errdown = 0.011;}
      else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.626; errup = 0.007; errdown = 0.007;}
      else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {eff = 0.729; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {eff = 0.584; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.462; errup = 0.045; errdown = 0.044;}
      else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.576; errup = 0.012; errdown = 0.012;}
      else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.660; errup = 0.008; errdown = 0.008;}
      else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {eff = 0.771; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {eff = 0.626; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.505; errup = 0.052; errdown = 0.052;}
      else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.575; errup = 0.013; errdown = 0.014;}
      else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.667; errup = 0.009; errdown = 0.009;}
      else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {eff = 0.796; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {eff = 0.657; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.462; errup = 0.058; errdown = 0.057;}
      else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.595; errup = 0.014; errdown = 0.015;}
      else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.705; errup = 0.009; errdown = 0.009;}
      else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {eff = 0.823; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {eff = 0.686; errup = 0.004; errdown = 0.004;}
      else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.500; errup = 0.057; errdown = 0.057;}
      else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.614; errup = 0.016; errdown = 0.016;}
      else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.722; errup = 0.010; errdown = 0.011;}
      else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {eff = 0.840; errup = 0.005; errdown = 0.005;}
      else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {eff = 0.711; errup = 0.004; errdown = 0.004;}
      else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.535; errup = 0.042; errdown = 0.043;}
      else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.620; errup = 0.011; errdown = 0.011;}
      else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.736; errup = 0.008; errdown = 0.008;}
      else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {eff = 0.859; errup = 0.003; errdown = 0.003;}
      else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {eff = 0.745; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.569; errup = 0.049; errdown = 0.051;}
      else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.617; errup = 0.014; errdown = 0.014;}
      else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.737; errup = 0.010; errdown = 0.010;}
      else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {eff = 0.882; errup = 0.004; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {eff = 0.777; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 300 && met<= 350) {eff = 0.558; errup = 0.040; errdown = 0.041;}
      else if(ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.638; errup = 0.013; errdown = 0.013;}
      else if(ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.770; errup = 0.010; errdown = 0.010;}
      else if(ht> 800 && ht<=1000 && met> 300 && met<= 350) {eff = 0.902; errup = 0.003; errdown = 0.004;}
      else if(ht>1000 && ht<=9999 && met> 300 && met<= 350) {eff = 0.804; errup = 0.003; errdown = 0.003;}
      else if(ht>   0 && ht<= 200 && met> 350 && met<= 400) {eff = 0.548; errup = 0.056; errdown = 0.057;}
      else if(ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.642; errup = 0.019; errdown = 0.019;}
      else if(ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.777; errup = 0.014; errdown = 0.015;}
      else if(ht> 800 && ht<=1000 && met> 350 && met<= 400) {eff = 0.925; errup = 0.005; errdown = 0.005;}
      else if(ht>1000 && ht<=9999 && met> 350 && met<= 400) {eff = 0.838; errup = 0.004; errdown = 0.004;}
      else if(ht>   0 && ht<= 200 && met> 400 && met<= 450) {eff = 0.557; errup = 0.070; errdown = 0.072;}
      else if(ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.710; errup = 0.024; errdown = 0.025;}
      else if(ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.864; errup = 0.018; errdown = 0.020;}
      else if(ht> 800 && ht<=1000 && met> 400 && met<= 450) {eff = 0.945; errup = 0.006; errdown = 0.006;}
      else if(ht>1000 && ht<=9999 && met> 400 && met<= 450) {eff = 0.862; errup = 0.005; errdown = 0.005;}
      else if(ht>   0 && ht<= 200 && met> 450 && met<= 500) {eff = 0.735; errup = 0.081; errdown = 0.097;}
      else if(ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.695; errup = 0.036; errdown = 0.039;}
      else if(ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.913; errup = 0.022; errdown = 0.028;}
      else if(ht> 800 && ht<=1000 && met> 450 && met<= 500) {eff = 0.970; errup = 0.006; errdown = 0.007;}
      else if(ht>1000 && ht<=9999 && met> 450 && met<= 500) {eff = 0.886; errup = 0.006; errdown = 0.006;}
      else if(ht>   0 && ht<= 200 && met> 500 && met<=9999) {eff = 0.550; errup = 0.088; errdown = 0.091;}
      else if(ht> 200 && ht<= 600 && met> 500 && met<=9999) {eff = 0.731; errup = 0.029; errdown = 0.031;}
      else if(ht> 600 && ht<= 800 && met> 500 && met<=9999) {eff = 0.956; errup = 0.014; errdown = 0.019;}
      else if(ht> 800 && ht<=1000 && met> 500 && met<=9999) {eff = 0.981; errup = 0.004; errdown = 0.005;}
      else if(ht>1000 && ht<=9999 && met> 500 && met<=9999) {eff = 0.907; errup = 0.004; errdown = 0.004;}

    } else { // TRUE MET
      if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.532; errup = 0.013; errdown = 0.013;}
      else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.612; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.589; errup = 0.023; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {eff = 0.515; errup = 0.042; errdown = 0.042;}
      else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {eff = 0.588; errup = 0.052; errdown = 0.054;}
      else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.591; errup = 0.014; errdown = 0.014;}
      else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.684; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.678; errup = 0.022; errdown = 0.022;}
      else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {eff = 0.537; errup = 0.042; errdown = 0.042;}
      else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {eff = 0.511; errup = 0.057; errdown = 0.057;}
      else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.619; errup = 0.016; errdown = 0.016;}
      else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.727; errup = 0.006; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.699; errup = 0.022; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {eff = 0.690; errup = 0.039; errdown = 0.041;}
      else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {eff = 0.568; errup = 0.048; errdown = 0.049;}
      else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.678; errup = 0.017; errdown = 0.018;}
      else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.769; errup = 0.006; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.732; errup = 0.022; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {eff = 0.609; errup = 0.042; errdown = 0.044;}
      else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {eff = 0.685; errup = 0.058; errdown = 0.064;}
      else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.670; errup = 0.019; errdown = 0.020;}
      else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.811; errup = 0.005; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.779; errup = 0.022; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {eff = 0.736; errup = 0.041; errdown = 0.045;}
      else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {eff = 0.663; errup = 0.056; errdown = 0.061;}
      else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.730; errup = 0.020; errdown = 0.021;}
      else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.838; errup = 0.005; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.820; errup = 0.022; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {eff = 0.819; errup = 0.037; errdown = 0.043;}
      else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {eff = 0.736; errup = 0.055; errdown = 0.062;}
      else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.745; errup = 0.023; errdown = 0.024;}
      else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.874; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.848; errup = 0.021; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {eff = 0.869; errup = 0.038; errdown = 0.048;}
      else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {eff = 0.759; errup = 0.048; errdown = 0.055;}
      else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.777; errup = 0.024; errdown = 0.026;}
      else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.903; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.850; errup = 0.021; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {eff = 0.839; errup = 0.041; errdown = 0.049;}
      else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {eff = 0.847; errup = 0.044; errdown = 0.055;}
      else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.792; errup = 0.026; errdown = 0.028;}
      else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.907; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.884; errup = 0.020; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {eff = 0.870; errup = 0.036; errdown = 0.045;}
      else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {eff = 0.781; errup = 0.051; errdown = 0.059;}
      else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.757; errup = 0.033; errdown = 0.036;}
      else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.924; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.921; errup = 0.016; errdown = 0.020;}
      else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {eff = 0.936; errup = 0.027; errdown = 0.041;}
      else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {eff = 0.803; errup = 0.049; errdown = 0.059;}
      else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.841; errup = 0.022; errdown = 0.025;}
      else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.949; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.927; errup = 0.013; errdown = 0.015;}
      else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {eff = 0.894; errup = 0.023; errdown = 0.027;}
      else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {eff = 0.839; errup = 0.036; errdown = 0.042;}
      else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.850; errup = 0.028; errdown = 0.032;}
      else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.966; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.952; errup = 0.011; errdown = 0.013;}
      else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {eff = 0.919; errup = 0.024; errdown = 0.031;}
      else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {eff = 0.959; errup = 0.018; errdown = 0.027;}
      else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.896; errup = 0.029; errdown = 0.037;}
      else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.973; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.979; errup = 0.008; errdown = 0.011;}
      else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {eff = 0.956; errup = 0.017; errdown = 0.025;}
      else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {eff = 0.971; errup = 0.016; errdown = 0.027;}
      else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.844; errup = 0.042; errdown = 0.053;}
      else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.983; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.976; errup = 0.009; errdown = 0.013;}
      else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {eff = 0.983; errup = 0.011; errdown = 0.022;}
      else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {eff = 0.942; errup = 0.028; errdown = 0.043;}
      else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.880; errup = 0.047; errdown = 0.065;}
      else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.985; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.992; errup = 0.005; errdown = 0.010;}
      else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {eff = 0.989; errup = 0.010; errdown = 0.026;}
      else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {eff = 0.931; errup = 0.030; errdown = 0.044;}
      else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.915; errup = 0.033; errdown = 0.047;}
      else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.989; errup = 0.002; errdown = 0.002;}
      else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.992; errup = 0.004; errdown = 0.006;}
      else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {eff = 0.984; errup = 0.009; errdown = 0.016;}
      else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {eff = 0.965; errup = 0.015; errdown = 0.023;}
      else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.862; errup = 0.065; errdown = 0.096;}
      else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.992; errup = 0.002; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.989; errup = 0.005; errdown = 0.008;}
      else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {eff = 0.963; errup = 0.016; errdown = 0.024;}
      else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {eff = 0.991; errup = 0.007; errdown = 0.020;}
      else if(ht>   0 && ht<= 200 && met> 300 && met<=9999) {eff = 0.744; errup = 0.074; errdown = 0.089;}
      else if(ht> 200 && ht<= 600 && met> 300 && met<=9999) {eff = 0.994; errup = 0.001; errdown = 0.002;}
      else if(ht> 600 && ht<= 800 && met> 300 && met<=9999) {eff = 0.996; errup = 0.002; errdown = 0.003;}
      else if(ht> 800 && ht<=1000 && met> 300 && met<=9999) {eff = 1.000; errup = 0.000; errdown = 0.003;}
      else if(ht>1000 && ht<=9999 && met> 300 && met<=9999) {eff = 0.987; errup = 0.005; errdown = 0.008;}
    }

    //// MET || Ele27 || Ele105 || Ele115
    //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]||trig[22]||trig[40]||trig[24]||trig[41])"
  // } else if(b.nels()==1 && b.nmus()==0){
  //   vector<float> leps_pt; 
  //   if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
  //   else leps_pt.push_back(0);
  //   if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met>=0 && met<= 110) {eff=0.160; errup=0.019; errdown = 0.017;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>=0 && met<= 110) {eff=0.400; errup=0.024; errdown=0.024;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>=0 && met<= 110) {eff=0.728; errup=0.006; errdown=0.006;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>=0 && met<= 110) {eff=0.880; errup=0.017; errdown=0.019;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>=0 && met<= 110) {eff=0.950; errup=0.003; errdown=0.003;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>110 && met<= 120) {eff=0.244; errup=0.024; errdown=0.023;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>110 && met<= 120) {eff=0.420; errup=0.027; errdown=0.027;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>110 && met<= 120) {eff=0.761; errup=0.007; errdown=0.007;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>110 && met<= 120) {eff=0.918; errup=0.015; errdown=0.017;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>110 && met<= 120) {eff=0.958; errup=0.003; errdown=0.003;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>120 && met<= 130) {eff=0.331; errup=0.030; errdown=0.029;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>120 && met<= 130) {eff=0.500; errup=0.031; errdown=0.031;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>120 && met<= 130) {eff=0.800; errup=0.007; errdown=0.007;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>120 && met<= 130) {eff=0.928; errup=0.015; errdown=0.018;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>120 && met<= 130) {eff=0.960; errup=0.003; errdown=0.003;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>130 && met<= 140) {eff=0.491; errup=0.031; errdown=0.031;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>130 && met<= 140) {eff=0.608; errup=0.033; errdown=0.034;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>130 && met<= 140) {eff=0.831; errup=0.007; errdown=0.007;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>130 && met<= 140) {eff=0.931; errup=0.016; errdown=0.020;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>130 && met<= 140) {eff=0.967; errup=0.003; errdown=0.003;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>140 && met<= 150) {eff=0.573; errup=0.033; errdown=0.033;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>140 && met<= 150) {eff=0.677; errup=0.035; errdown=0.037;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>140 && met<= 150) {eff=0.856; errup=0.007; errdown=0.007;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>140 && met<= 150) {eff=0.923; errup=0.018; errdown=0.022;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>140 && met<= 150) {eff=0.971; errup=0.003; errdown=0.004;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>150 && met<= 160) {eff=0.643; errup=0.037; errdown=0.039;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 150 && met<= 160) {eff=0.738; errup=0.033; errdown=0.036;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 150 && met<= 160) {eff=0.871; errup=0.007; errdown=0.007;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 150 && met<= 160) {eff=0.935; errup=0.016; errdown=0.021;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 150 && met<= 160) {eff=0.982; errup=0.003; errdown=0.003;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 160 && met<= 170) {eff=0.760; errup=0.034; errdown=0.038;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 160 && met<= 170) {eff=0.792; errup=0.034; errdown=0.038;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 160 && met<= 170) {eff=0.910; errup=0.007; errdown=0.007;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 160 && met<= 170) {eff=0.970; errup=0.013; errdown=0.020;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 160 && met<= 170) {eff=0.981; errup=0.003; errdown=0.004;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 170 && met<= 180) {eff=0.829; errup=0.033; errdown=0.038;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 170 && met<= 180) {eff=0.863; errup=0.031; errdown=0.037;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 170 && met<= 180) {eff=0.937; errup=0.006; errdown=0.006;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 170 && met<= 180) {eff=0.988; errup=0.008; errdown=0.016;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 170 && met<= 180) {eff=0.982; errup=0.003; errdown=0.004;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 180 && met<= 190) {eff=0.761; errup=0.041; errdown=0.046;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 180 && met<= 190) {eff=0.863; errup=0.032; errdown=0.038;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 180 && met<= 190) {eff=0.939; errup=0.006; errdown=0.007;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 180 && met<= 190) {eff=0.969; errup=0.015; errdown=0.024;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 180 && met<= 190) {eff=0.984; errup=0.003; errdown=0.004;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 190 && met<= 200) {eff=0.892; errup=0.033; errdown=0.042;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 190 && met<= 200) {eff=0.902; errup=0.030; errdown=0.039;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 190 && met<= 200) {eff=0.956; errup=0.006; errdown=0.006;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 190 && met<= 200) {eff=0.983; errup=0.011; errdown=0.022;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 190 && met<= 200) {eff=0.993; errup=0.002; errdown=0.003;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 200 && met<= 210) {eff=0.950; errup=0.021; errdown=0.032;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 200 && met<= 210) {eff=0.951; errup=0.023; errdown=0.037;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 200 && met<= 210) {eff=0.973; errup=0.005; errdown=0.005;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 200 && met<= 210) {eff=1.000; errup=0.000; errdown=0.018;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 200 && met<= 210) {eff=0.985; errup=0.003; errdown=0.004;}
  //   else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 210 && met<=9999) {eff=0.974; errup=0.005; errdown=0.006;}
  //   else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 210 && met<=9999) {eff=0.981; errup=0.004; errdown=0.005;}
  //   else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 210 && met<=9999) {eff=0.992; errup=0.001; errdown=0.001;}
  //   else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 210 && met<=9999) {eff=0.997; errup=0.002; errdown=0.003;}
  //   else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 210 && met<=9999) {eff=0.996; errup=0.001; errdown=0.001;}

  //   //// MET || Mu24 || Mu50
  //   //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]||trig[19]||trig[55]||trig[21])"
  // } else if(b.nels()==0 && b.nmus()==1){
  //   vector<float> leps_pt; 
  //   if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
  //   else leps_pt.push_back(0);
  //   if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met>= 0 && met<= 110) {eff=0.271; errup=0.017; errdown = 0.016;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met>= 0 && met<= 110) {eff=0.725; errup=0.017; errdown=0.018;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met>= 0 && met<= 110) {eff=0.814; errup=0.008; errdown=0.009;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met>= 0 && met<= 110) {eff=0.964; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 110 && met<= 120) {eff=0.363; errup=0.020; errdown=0.020;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 110 && met<= 120) {eff=0.755; errup=0.018; errdown=0.019;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 110 && met<= 120) {eff=0.842; errup=0.009; errdown=0.009;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 110 && met<= 120) {eff=0.969; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 120 && met<= 130) {eff=0.452; errup=0.022; errdown=0.022;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 120 && met<= 130) {eff=0.824; errup=0.018; errdown=0.019;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 120 && met<= 130) {eff=0.869; errup=0.009; errdown=0.009;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 120 && met<= 130) {eff=0.971; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 130 && met<= 140) {eff=0.590; errup=0.025; errdown=0.025;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 130 && met<= 140) {eff=0.875; errup=0.017; errdown=0.019;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 130 && met<= 140) {eff=0.904; errup=0.008; errdown=0.009;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 130 && met<= 140) {eff=0.972; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 140 && met<= 150) {eff=0.660; errup=0.026; errdown=0.027;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 140 && met<= 150) {eff=0.891; errup=0.017; errdown=0.019;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 140 && met<= 150) {eff=0.938; errup=0.007; errdown=0.008;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 140 && met<= 150) {eff=0.980; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 150 && met<= 160) {eff=0.778; errup=0.024; errdown=0.026;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 150 && met<= 160) {eff=0.915; errup=0.016; errdown=0.019;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 150 && met<= 160) {eff=0.940; errup=0.008; errdown=0.009;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 150 && met<= 160) {eff=0.984; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 160 && met<= 170) {eff=0.798; errup=0.026; errdown=0.029;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 160 && met<= 170) {eff=0.946; errup=0.015; errdown=0.020;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 160 && met<= 170) {eff=0.967; errup=0.006; errdown=0.008;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 160 && met<= 170) {eff=0.991; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 170 && met<= 180) {eff=0.885; errup=0.022; errdown=0.025;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 170 && met<= 180) {eff=0.937; errup=0.016; errdown=0.021;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 170 && met<= 180) {eff=0.977; errup=0.006; errdown=0.007;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 170 && met<= 180) {eff=0.987; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 180 && met<= 190) {eff=0.927; errup=0.019; errdown=0.024;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 180 && met<= 190) {eff=0.958; errup=0.014; errdown=0.019;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 180 && met<= 190) {eff=0.974; errup=0.006; errdown=0.008;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 180 && met<= 190) {eff=0.992; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 190 && met<= 200) {eff=0.921; errup=0.019; errdown=0.024;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 190 && met<= 200) {eff=0.965; errup=0.014; errdown=0.020;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 190 && met<= 200) {eff=0.991; errup=0.004; errdown=0.006;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 190 && met<= 200) {eff=0.991; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 200 && met<= 210) {eff=0.926; errup=0.022; errdown=0.028;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 200 && met<= 210) {eff=0.994; errup=0.005; errdown=0.015;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 200 && met<= 210) {eff=0.994; errup=0.003; errdown=0.006;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 200 && met<= 210) {eff=0.994; errup=0.002; errdown=0.002;}
  //   else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 210 && met<=9999) {eff=0.981; errup=0.004; errdown=0.004;}
  //   else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 210 && met<=9999) {eff=0.994; errup=0.002; errdown=0.003;}
  //   else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 210 && met<=9999) {eff=0.996; errup=0.001; errdown=0.001;}
  //   else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 210 && met<=9999) {eff=0.997; errup=0.000; errdown=0.000;}

  //   //// Ele27 || Ele105 || Ele115
  //   //// "(trig[22]||trig[40]||trig[24]||trig[41])"
  // } else if(b.nels()==2 && b.nmus()==0){
  //   vector<float> leps_pt; 
  //   if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
  //   else leps_pt.push_back(0);
  //   if(leps_pt[0]>  40 && leps_pt[0]<=  45) {eff = 0.944; errup = 0.015; errdown = 0.019;}
  //   else if(leps_pt[0]>  45 && leps_pt[0]<=  50) {eff = 0.910; errup = 0.015; errdown = 0.017;}
  //   else if(leps_pt[0]>  50 && leps_pt[0]<=  55) {eff = 0.927; errup = 0.013; errdown = 0.015;}
  //   else if(leps_pt[0]>  55 && leps_pt[0]<=  60) {eff = 0.912; errup = 0.013; errdown = 0.015;}
  //   else if(leps_pt[0]>  60 && leps_pt[0]<=  65) {eff = 0.941; errup = 0.011; errdown = 0.013;}
  //   else if(leps_pt[0]>  65 && leps_pt[0]<=  70) {eff = 0.901; errup = 0.014; errdown = 0.016;}
  //   else if(leps_pt[0]>  70 && leps_pt[0]<=  75) {eff = 0.921; errup = 0.013; errdown = 0.016;}
  //   else if(leps_pt[0]>  75 && leps_pt[0]<=  80) {eff = 0.947; errup = 0.011; errdown = 0.014;}
  //   else if(leps_pt[0]>  80 && leps_pt[0]<=  85) {eff = 0.954; errup = 0.011; errdown = 0.013;}
  //   else if(leps_pt[0]>  85 && leps_pt[0]<=  90) {eff = 0.939; errup = 0.012; errdown = 0.014;}
  //   else if(leps_pt[0]>  90 && leps_pt[0]<=  95) {eff = 0.940; errup = 0.012; errdown = 0.015;}
  //   else if(leps_pt[0]>  95 && leps_pt[0]<= 100) {eff = 0.932; errup = 0.014; errdown = 0.017;}
  //   else if(leps_pt[0]> 100 && leps_pt[0]<= 105) {eff = 0.934; errup = 0.014; errdown = 0.017;}
  //   else if(leps_pt[0]> 105 && leps_pt[0]<= 110) {eff = 0.965; errup = 0.010; errdown = 0.014;}
  //   else if(leps_pt[0]> 110 && leps_pt[0]<=9999) {eff = 0.994; errup = 0.001; errdown = 0.001;}

  //   //// Mu24 || Mu50
  //   //// "(trig[19]||trig[55]||trig[21])"
  // } else if(b.nels()==0 && b.nmus()==2){
  //   vector<float> leps_pt; 
  //   if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
  //   else leps_pt.push_back(0);
  //   if(leps_pt[0]>  40 && leps_pt[0]<=  45) {eff = 0.959; errup = 0.010; errdown = 0.012;}
  //   if(leps_pt[0]>  45 && leps_pt[0]<=  50) {eff = 0.970; errup = 0.006; errdown = 0.007;}
  //   if(leps_pt[0]>  50 && leps_pt[0]<=9999) {eff = 0.982; errup = 0.000; errdown = 0.000;}
    }
    return eff;
  });

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
