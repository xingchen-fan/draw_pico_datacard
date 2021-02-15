#include "higgsino/hig_functions.hpp"
#include "core/functions.hpp"
#include "higgsino/apply_trigeffs_v0.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"

#include "TVector2.h"
#include "TMath.h"
#include "Math/Vector4D.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace Higfuncs{

const NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
  int tmp_ntrub(0);
  for (unsigned i(0); i<b.jet_pt()->size(); i++){
    if (!b.jet_h1d()->at(i) && !b.jet_h2d()->at(i)) continue;
    if (b.jet_hflavor()->at(i)==5) tmp_ntrub++;
  }
  return tmp_ntrub;
});

const NamedFunc nbfake("nbfake",[](const Baby &b) -> NamedFunc::ScalarType{
  int nb_tru(0);
  for (unsigned i(0); i<b.jet_pt()->size(); i++){
    if (!b.jet_h1d()->at(i) && !b.jet_h2d()->at(i)) continue;
    if (b.jet_hflavor()->at(i)==5) nb_tru++;
  }
  int nb_reco(0);
  if (b.nbt()>=2) {
    if (b.nbm()==2) nb_reco = 2;
    if (b.nbm()==3&&b.nbl()==3) nb_reco = 3;
    if (b.nbm()>=3&&b.nbl()>=4) nb_reco = 4;
  }
  return nb_reco-nb_tru;
});

const NamedFunc stitch_htmet("stitch_htmet",[](const Baby &b) -> NamedFunc::ScalarType{
  bool stitch_(b.stitch());
  if (b.ht_isr_me()>600 && b.type()>=1000 && b.type()<1300) 
    stitch_ = false; 
  return stitch_;
});

const NamedFunc hig_pt1("hig_pt1",[](const Baby &b) -> NamedFunc::ScalarType{
    float higpt(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)!=25) continue;
      if (b.mc_pt()->at(i)>higpt) higpt = b.mc_pt()->at(i);
    }
    return higpt;
});

const NamedFunc hig_pt2("hig_pt2",[](const Baby &b) -> NamedFunc::ScalarType{
    float higpt1(0), higpt2(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)!=25) continue;
      if (b.mc_pt()->at(i)>higpt1) {
	higpt2 = higpt1;
	higpt1 = b.mc_pt()->at(i);
      } else if (b.mc_pt()->at(i)>higpt2) {
	higpt2 = b.mc_pt()->at(i);
      }
    }
    return higpt2;
});

const NamedFunc hig_bcat("hig_bcat",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbt()==2 && b.nbm()==2) return 2;
  else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) return 3;
  else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) return 4;
  else return 0;
});

const NamedFunc higd_bcat("higd_bcat",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbt()==2 && b.nbm()==2) return 2;
  else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) return 3;
  else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) return 4;
  else return 0;
});

const NamedFunc higd_bcat_extended("higd_bcat_extended",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbm()==0) return 0;
  else if (b.nbm()==1) return 1;
  else if (b.nbt()==2 && b.nbm()==2) return 2;
  else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) return 3;
  else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) return 4;
  else return 99;
});

const NamedFunc higd_bcat_mmmm("higd_bcat_mmmm",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbm()>=2) return min(4,b.nbm());
  else return 0;
});
const NamedFunc higd_bcat_ttll("higd_bcat_ttll",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbt()==2 && b.nbl()==2) return 2;
  else if (b.nbt()>=2 && b.nbl()==3) return 3;
  else if (b.nbt()>=2 && b.nbl()>=4) return 4;
  else return 0;
});
const NamedFunc higd_bcat_tmml("higd_bcat_tmml",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbt()>=1 && b.nbm()==2) return 2;
  else if (b.nbt()>=1 && b.nbm()==3 && b.nbl()==3) return 3;
  else if (b.nbt()>=1 && b.nbm()>=3 && b.nbl()>=4) return 4;
  else return 0;
});

// apply weights found from the nb and MET data/MC comparisons
const NamedFunc wgt_comp("wgt_comp",[](const Baby &b) -> NamedFunc::ScalarType{
  float wgt = 1;
  if ( (b.type()>=1000 && b.type()<2000) ||  // ttbar
    // (b.type()>=3000 && b.type()<4000) ||     // single top
    (b.type()>=4000 && b.type()<6000) ||     // ttw and ttz
    b.type()==9000  ||                       // ttH
    b.type()==10000  ||                      // ttgamma
    b.type()==11000) {                       // tttt
    //apply normalization factor from MET distribution
    wgt = 1.052;
    if (b.met()<=50)                      wgt*=1.005;
    else if (b.met()>50 && b.met()<=100)  wgt*=1.018;
    else if (b.met()>100 && b.met()<=150) wgt*=0.976;
    else if (b.met()>150 && b.met()<=200) wgt*=0.914;
    else if (b.met()>200 && b.met()<=300) wgt*=0.854;
    else if (b.met()>300)                 wgt*=0.865;
    // nb correction from nb distribution data/mc, inclusive in MET
    if (b.nbt()==2 && b.nbm()==2)                     wgt*=0.992;
    else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) wgt*=1.040;
    else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) wgt*=1.165;
  } 
  else if((b.type()>=8000 && b.type()<9000) || // zjets
  (b.type()>=2000 && b.type()<3000) ||      // wjets
  (b.type()>=6000 && b.type()<7000)) {      // dyjets
    //apply normalization factor from MET distribution
    wgt = 1.416;
    if (b.met()>150 && b.met()<=200)      wgt*=1.121;
    else if (b.met()>200 && b.met()<=300) wgt*=0.951;
    else if (b.met()>300)                 wgt*=0.722;
    // nb correction from nb distribution data/mc, inclusive in MET
    if (b.nbt()==2 && b.nbm()==2)                     wgt*=0.985;
    else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) wgt*=1.176;
    else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) wgt*=1.097;
  } 
  else if ( (b.type()>=7000 && b.type()<8000)) { // qcd
  //apply normalization factor from MET distribution
    wgt = 1.700;
    if (b.met()>150 && b.met()<=200)      wgt*=0.927;
    else if (b.met()>200 && b.met()<=300) wgt*=1.199;
    else if (b.met()>300)                 wgt*=1.301;
    // nb correction from nb distribution data/mc, inclusive in MET
    if (b.nbt()==2 && b.nbm()==2)                    wgt*=0.982;
    else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) wgt*=1.142;
    else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) wgt*=1.069;
  }
  return wgt;
});

//// subtract ttbar based on MC prediction reweighted to data in 1l CR
//// since ttbar has to be combined in the same process def with data, 
//// also apply stitch, json and trigger here
//const NamedFunc wgt_subtr_ttx("wgt_subtr_ttx",[](const Baby &b) -> NamedFunc::ScalarType{
//  if ( (b.type()>=1000 && b.type()<2000) ||  // ttbar
//    // (b.type()>=3000 && b.type()<4000) ||     // single top
//    (b.type()>=4000 && b.type()<6000) ||     // ttw and ttz
//    b.type()==9000  ||                       // ttH
//    b.type()==10000  ||                      // ttgamma
//    b.type()==11000) {
//    float wgt = -1; 
//    // apply lumi 
//    wgt*= 35.9;
//    // apply ttx weights derived from 1l CR
//    wgt = 1.052;
//    if (b.met()<=50)                      wgt*=1.005;
//    else if (b.met()>50 && b.met()<=100)  wgt*=1.018;
//    else if (b.met()>100 && b.met()<=150) wgt*=0.976;
//    else if (b.met()>150 && b.met()<=200) wgt*=0.914;
//    else if (b.met()>200 && b.met()<=300) wgt*=0.854;
//    else if (b.met()>300)                 wgt*=0.865;
//    // nb correction from nb distribution data/mc, inclusive in MET
//    if (b.nbt()==2 && b.nbm()==2)                     wgt*=0.992;
//    else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) wgt*=1.040;
//    else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) wgt*=1.165;
//    // apply stitch
//    if (b.stitch()) return wgt;
//    else return 0;
//  } else if (b.type()>0 && b.type()<1000){ // apply trigger and json for data
//    return trig_hig_decision(b);
//  }
//  // for all other backgrounds, chill (they are not in the "data" process so no need to apply lumi)
//  return 1;
//});

// calculate effect of systematics calculated for each background 
// in the data control regions on the total bkg. kappa
const NamedFunc wgt_syst_ttx("wgt_syst_ttx",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=1000 && b.type()<2000) ||  // ttbar
    // (b.type()>=3000 && b.type()<4000) ||     // single top
    (b.type()>=4000 && b.type()<6000) ||     // ttw and ttz
    b.type()==9000  ||                       // ttH
    b.type()==10000  ||                      // ttgamma
    b.type()==11000) {                       // tttt
    if ((*b.hig_cand_am())[0]<=100 || ((*b.hig_cand_am())[0]>140 && (*b.hig_cand_am())[0]<=200)) {
      if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3)  {// 3b
        if ((*b.hig_cand_drmax())[0]<=1.1) return 0.13; // low dr
        else return 0.02; // high dr
      } else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) {// 4b
        if ((*b.hig_cand_drmax())[0]<=1.1) return 0.16; // low dr 
        else return 0.07; // high dr
      }
    }
  }
  return 0;
});

const NamedFunc wgt_syst_vjets("wgt_syst_vjets",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=8000 && b.type()<9000) || // zjets
    (b.type()>=2000 && b.type()<3000) ||    // wjets
    (b.type()>=6000 && b.type()<7000)) {   // dyjets
    if ((*b.hig_cand_am())[0]<=100 || ((*b.hig_cand_am())[0]>140 && (*b.hig_cand_am())[0]<=200))
      if (b.nbt()>=2 && b.nbm()>=3) { // high b
        if ((*b.hig_cand_drmax())[0]<=1.1) return 0.16; // low dr
        else return 0.05; // high dr
      }
  }
  return 0;
});

const NamedFunc wgt_syst_qcd("wgt_syst_qcd",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=7000 && b.type()<8000)) { // qcd
    if ((*b.hig_cand_am())[0]<=100 || ((*b.hig_cand_am())[0]>140 && (*b.hig_cand_am())[0]<=200))
      if (b.nbt()>=2 && b.nbm()>=3) { // high b
        if ((*b.hig_cand_drmax())[0]<=1.1) return 0.11; // low dr
        else return 0.09; // high dr
      }
  }
  return 0;
});

//return weight from composition systematic
//note this returns the weight scaling, not the variation
const NamedFunc wgt_syst_comp_unnormed("wgt_syst_comp_unnormed",[](const Baby &b) -> NamedFunc::ScalarType{
  double weight = 1.0;
  int sample = b.type()/1000;
  switch (sample) {
    case 1:
      if (Higfuncs::jetid_nb.GetScalar(b) > 1.5 && Higfuncs::jetid_nb.GetScalar(b) < 2.5) 
        weight = weight*1.0588;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 2.5 && Higfuncs::jetid_nb.GetScalar(b) < 3.5) 
        weight = weight*1.12644;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 3.5 && Higfuncs::jetid_nb.GetScalar(b) < 4.5) 
        weight = weight*1.09696;
      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 0 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 1.1) 
        weight = weight*1.04028;
      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 1.1 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 2.2) 
        weight = weight*1.06802;
      if (b.met() > 150 && b.met() < 200) 
        weight = weight*0.769807;
      else if (b.met() > 200 && b.met() < 300) 
        weight = weight*0.757657;
      else if (b.met() > 300 && b.met() < 400) 
        weight = weight*0.677293;
      else if (b.met() > 400 && b.met() < 9999) 
        weight = weight*0.668484;
      break;
    case 8:
      if (Higfuncs::jetid_nb.GetScalar(b) > 1.5 && Higfuncs::jetid_nb.GetScalar(b) < 2.5) 
        weight = weight*0.845564;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 2.5 && Higfuncs::jetid_nb.GetScalar(b) < 3.5) 
        weight = weight*1.77212;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 3.5 && Higfuncs::jetid_nb.GetScalar(b) < 4.5) 
        weight = weight*1.19877;
      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 0 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 1.1) 
        weight = weight*0.89623;
      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 1.1 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 2.2) 
        weight = weight*0.928458;
      if (b.met() > 150 && b.met() < 200) 
        weight = weight*0.955291;
      else if (b.met() > 200 && b.met() < 300) 
        weight = weight*0.863349;
      else if (b.met() > 300 && b.met() < 400) 
        weight = weight*0.827839;
      else if (b.met() > 400 && b.met() < 9999) 
        weight = weight*0.659814;
      break;
    case 7:
      if (Higfuncs::jetid_nb.GetScalar(b) > 1.5 && Higfuncs::jetid_nb.GetScalar(b) < 2.5) 
        weight = weight*1.28086;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 2.5 && Higfuncs::jetid_nb.GetScalar(b) < 3.5) 
        weight = weight*1.22963;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 3.5 && Higfuncs::jetid_nb.GetScalar(b) < 4.5) 
        weight = weight*1.00255;
      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 0 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 1.1) 
        weight = weight*1.18912;
      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 1.1 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 2.2) 
        weight = weight*1.41141;
      if (b.met() > 150 && b.met() < 200) 
        weight = weight*1.34946;
      else if (b.met() > 200 && b.met() < 300) 
        weight = weight*1.25167;
      else if (b.met() > 300 && b.met() < 400) 
        weight = weight*1.07007;
      else if (b.met() > 400 && b.met() < 9999) 
        weight = weight*0.844552;
    default:
      break;
  }
  return weight;
});

//return weight from composition systematic
//note this returns the weight scaling, not the variation
const NamedFunc wgt_syst_comp("wgt_syst_comp",[](const Baby &b) -> NamedFunc::ScalarType{
  double weight = 1.0;
  int sample = b.type()/1000;
  switch (sample) {
    case 1:
      if (Higfuncs::jetid_nb.GetScalar(b) > 1.5 && Higfuncs::jetid_nb.GetScalar(b) < 2.5) 
        weight = weight*0.994324;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 2.5 && Higfuncs::jetid_nb.GetScalar(b) < 3.5) 
        weight = weight*1.05784;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 3.5 && Higfuncs::jetid_nb.GetScalar(b) < 4.5) 
        weight = weight*1.03016;
      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 0 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 1.1) 
        weight = weight*0.97693;
      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 1.1 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 2.2) 
        weight = weight*1.00298;
      if (b.met() > 150 && b.met() < 200) 
        weight = weight*0.769807;
      else if (b.met() > 200 && b.met() < 300) 
        weight = weight*0.757657;
      else if (b.met() > 300 && b.met() < 400) 
        weight = weight*0.677293;
      else if (b.met() > 400 && b.met() < 9999) 
        weight = weight*0.668484;
      break;
    case 8:
      if (Higfuncs::jetid_nb.GetScalar(b) > 1.5 && Higfuncs::jetid_nb.GetScalar(b) < 2.5) 
        weight = weight*0.928212;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 2.5 && Higfuncs::jetid_nb.GetScalar(b) < 3.5) 
        weight = weight*1.94533;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 3.5 && Higfuncs::jetid_nb.GetScalar(b) < 4.5) 
        weight = weight*1.31595;
      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 0 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 1.1) 
        weight = weight*0.974129;
      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 1.1 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 2.2) 
        weight = weight*1.00916;
      if (b.met() > 150 && b.met() < 200) 
        weight = weight*0.955291;
      else if (b.met() > 200 && b.met() < 300) 
        weight = weight*0.863349;
      else if (b.met() > 300 && b.met() < 400) 
        weight = weight*0.827839;
      else if (b.met() > 400 && b.met() < 9999) 
        weight = weight*0.659814;
      break;
    case 7:
      if (Higfuncs::jetid_nb.GetScalar(b) > 1.5 && Higfuncs::jetid_nb.GetScalar(b) < 2.5) 
        weight = weight*1.01916;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 2.5 && Higfuncs::jetid_nb.GetScalar(b) < 3.5) 
        weight = weight*0.978398;
      else if (Higfuncs::jetid_nb.GetScalar(b) > 3.5 && Higfuncs::jetid_nb.GetScalar(b) < 4.5) 
        weight = weight*0.797707;
      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 0 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 1.1) 
        weight = weight*0.906144;
      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > 1.1 && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < 2.2) 
        weight = weight*1.07554;
      if (b.met() > 150 && b.met() < 200) 
        weight = weight*1.34946;
      else if (b.met() > 200 && b.met() < 300) 
        weight = weight*1.25167;
      else if (b.met() > 300 && b.met() < 400) 
        weight = weight*1.07007;
      else if (b.met() > 400 && b.met() < 9999) 
        weight = weight*0.844552;
    default:
      break;
  }
  return weight;
});

const NamedFunc wgt_syst_lumi_up("wgt_syst_lumi_up",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType() < 0) return 1.0; //data
    else if (b.SampleType()==2017) return 1.023;
    return 1.025;
});

const NamedFunc wgt_syst_lumi_down("wgt_syst_lumi_down",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType() < 0) return 1.0; //data
    else if (b.SampleType()==2017) return 0.977;
    return 0.975;
});

//// Definition of analysis trigger
// NamedFunc::ScalarType trig_hig_decision(const Baby &b){
//    bool mettrig = b.trig()->at(13)||b.trig()->at(33)||b.trig()->at(14)||b.trig()->at(15)
//      ||b.trig()->at(30)||b.trig()->at(31);
//    bool eltrig = b.trig()->at(22)||b.trig()->at(40)||b.trig()->at(24)||b.trig()->at(41);
//    bool mutrig = b.trig()->at(19)||b.trig()->at(55)||b.trig()->at(21);

//    if(b.nel()==1 && b.nmu()==0){
//      if(mettrig || eltrig) return 1;
//      else return -1;
//    } else if(b.nel()==0 && b.nmu()==1){
//      if(mettrig || mutrig) return 1;
//      else return -1;
//    } else if(b.nel()==2 && b.nmu()==0){
//      if(eltrig) return 1;
//      else return -1;
//    } else if(b.nel()==0 && b.nmu()==2){
//      if(mutrig) return 1;
//      else return -1;
//    } else if(b.nel()==1 && b.nmu()==1){
//      if(eltrig || mutrig) return 1;
//      else return -1;
//    } else if(b.nvlep()==0){
//      if(mettrig) return 1;
//      else return -1;
//    }

//    return -1;
// }

// const NamedFunc trig_hig("trig_hig", [](const Baby &b) -> NamedFunc::ScalarType{
//  return trig_hig_decision(b);
//  });
  
// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
const NamedFunc err_higtrig("err_higtrig", [](const Baby &b) -> NamedFunc::ScalarType{
    float errup, errdown; // Stat uncertainty. Not used, but for reference
    float uncert = 0., met = b.met(), ht = b.ht();
    errup=0;errdown=0;
    errup+=errdown;

    if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {uncert = 0.072; errup = 0.013; errdown = 0.013;}
    else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {uncert = 0.071; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {uncert = 0.075; errup = 0.023; errdown = 0.024;}
    else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {uncert = 0.083; errup = 0.042; errdown = 0.042;}
    else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {uncert = 0.089; errup = 0.052; errdown = 0.054;}
    else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {uncert = 0.072; errup = 0.014; errdown = 0.014;}
    else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {uncert = 0.071; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {uncert = 0.074; errup = 0.022; errdown = 0.022;}
    else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {uncert = 0.083; errup = 0.042; errdown = 0.042;}
    else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {uncert = 0.091; errup = 0.057; errdown = 0.057;}
    else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {uncert = 0.047; errup = 0.016; errdown = 0.016;}
    else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {uncert = 0.044; errup = 0.006; errdown = 0.006;}
    else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {uncert = 0.050; errup = 0.022; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {uncert = 0.060; errup = 0.039; errdown = 0.041;}
    else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {uncert = 0.066; errup = 0.048; errdown = 0.049;}
    else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {uncert = 0.047; errup = 0.017; errdown = 0.018;}
    else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {uncert = 0.044; errup = 0.006; errdown = 0.006;}
    else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {uncert = 0.050; errup = 0.022; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {uncert = 0.062; errup = 0.042; errdown = 0.044;}
    else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {uncert = 0.077; errup = 0.058; errdown = 0.064;}
    else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {uncert = 0.048; errup = 0.019; errdown = 0.020;}
    else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {uncert = 0.044; errup = 0.005; errdown = 0.006;}
    else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {uncert = 0.050; errup = 0.022; errdown = 0.024;}
    else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {uncert = 0.063; errup = 0.041; errdown = 0.045;}
    else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {uncert = 0.075; errup = 0.056; errdown = 0.061;}
    else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {uncert = 0.049; errup = 0.020; errdown = 0.021;}
    else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {uncert = 0.044; errup = 0.005; errdown = 0.006;}
    else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {uncert = 0.050; errup = 0.022; errdown = 0.024;}
    else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {uncert = 0.062; errup = 0.037; errdown = 0.043;}
    else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {uncert = 0.076; errup = 0.055; errdown = 0.062;}
    else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {uncert = 0.049; errup = 0.023; errdown = 0.024;}
    else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {uncert = 0.048; errup = 0.021; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {uncert = 0.064; errup = 0.038; errdown = 0.048;}
    else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {uncert = 0.069; errup = 0.048; errdown = 0.055;}
    else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {uncert = 0.049; errup = 0.024; errdown = 0.026;}
    else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {uncert = 0.048; errup = 0.021; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {uncert = 0.065; errup = 0.041; errdown = 0.049;}
    else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {uncert = 0.069; errup = 0.044; errdown = 0.055;}
    else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {uncert = 0.051; errup = 0.026; errdown = 0.028;}
    else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {uncert = 0.048; errup = 0.020; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {uncert = 0.062; errup = 0.036; errdown = 0.045;}
    else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {uncert = 0.073; errup = 0.051; errdown = 0.059;}
    else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {uncert = 0.055; errup = 0.033; errdown = 0.036;}
    else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {uncert = 0.046; errup = 0.016; errdown = 0.020;}
    else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {uncert = 0.059; errup = 0.027; errdown = 0.041;}
    else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {uncert = 0.072; errup = 0.049; errdown = 0.059;}
    else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {uncert = 0.035; errup = 0.022; errdown = 0.025;}
    else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {uncert = 0.028; errup = 0.013; errdown = 0.015;}
    else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {uncert = 0.036; errup = 0.023; errdown = 0.027;}
    else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {uncert = 0.049; errup = 0.036; errdown = 0.042;}
    else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {uncert = 0.040; errup = 0.028; errdown = 0.032;}
    else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {uncert = 0.027; errup = 0.011; errdown = 0.013;}
    else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {uncert = 0.039; errup = 0.024; errdown = 0.031;}
    else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {uncert = 0.036; errup = 0.018; errdown = 0.027;}
    else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {uncert = 0.044; errup = 0.029; errdown = 0.037;}
    else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {uncert = 0.026; errup = 0.008; errdown = 0.011;}
    else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {uncert = 0.035; errup = 0.017; errdown = 0.025;}
    else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {uncert = 0.036; errup = 0.016; errdown = 0.027;}
    else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {uncert = 0.055; errup = 0.042; errdown = 0.053;}
    else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {uncert = 0.015; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {uncert = 0.020; errup = 0.009; errdown = 0.013;}
    else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {uncert = 0.027; errup = 0.011; errdown = 0.022;}
    else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {uncert = 0.046; errup = 0.028; errdown = 0.043;}
    else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {uncert = 0.067; errup = 0.047; errdown = 0.065;}
    else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {uncert = 0.015; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {uncert = 0.018; errup = 0.005; errdown = 0.010;}
    else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {uncert = 0.030; errup = 0.010; errdown = 0.026;}
    else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {uncert = 0.047; errup = 0.030; errdown = 0.044;}
    else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {uncert = 0.049; errup = 0.033; errdown = 0.047;}
    else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {uncert = 0.012; errup = 0.002; errdown = 0.002;}
    else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {uncert = 0.013; errup = 0.004; errdown = 0.006;}
    else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {uncert = 0.020; errup = 0.009; errdown = 0.016;}
    else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {uncert = 0.026; errup = 0.015; errdown = 0.023;}
    else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {uncert = 0.096; errup = 0.065; errdown = 0.096;}
    else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {uncert = 0.012; errup = 0.002; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {uncert = 0.015; errup = 0.005; errdown = 0.008;}
    else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {uncert = 0.027; errup = 0.016; errdown = 0.024;}
    else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {uncert = 0.023; errup = 0.007; errdown = 0.020;}
    else if(ht>   0 && ht<= 200 && met> 300 && met<=9999) {uncert = 0.089; errup = 0.074; errdown = 0.089;}
    else if(ht> 200 && ht<= 600 && met> 300 && met<=9999) {uncert = 0.006; errup = 0.001; errdown = 0.002;}
    else if(ht> 600 && ht<= 800 && met> 300 && met<=9999) {uncert = 0.007; errup = 0.002; errdown = 0.003;}
    else if(ht> 800 && ht<=1000 && met> 300 && met<=9999) {uncert = 0.007; errup = 0.000; errdown = 0.003;}
    else if(ht>1000 && ht<=9999 && met> 300 && met<=9999) {uncert = 0.010; errup = 0.005; errdown = 0.008;}

    return uncert;
  });
  
const NamedFunc eff_higtrig_run2("eff_higtrig_run2", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup, errdown; // Not used, but for reference
  float eff = 1.;
  errup=0;errdown=0;
  errup+=errdown;
  if(b.type()>0 && b.type()<1000) eff = 1; // data

  else if(b.nvlep()==0){ // search MC sample and qcd MC control sample
    if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
      if (b.SampleType()==2016) eff = get_0l_fakemet_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_0l_fakemet_trigeff2017.GetVector(b)[0];
      else if (b.SampleType()==2018) eff = get_0l_fakemet_trigeff2018.GetVector(b)[0];
    } else { // TRUE MET
      if (b.SampleType()==2016) eff = get_0l_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_0l_trigeff2017.GetVector(b)[0];
      else if (b.SampleType()==2018) eff = get_0l_trigeff2018.GetVector(b)[0];
    }
  } else if (b.nel()==1 && b.nmu()==0) { // 1 lepton MC control sample
    if (b.SampleType()==2016) eff = get_1el_trigeff2016.GetVector(b)[0];
    else if (b.SampleType()==2017) eff = get_1el_trigeff2017.GetVector(b)[0];
    else if (b.SampleType()==2018) eff = get_1el_trigeff2018.GetVector(b)[0];
  } else if (b.nel()==0 && b.nmu()==1) { // 1 lepton MC control sample
    if (b.SampleType()==2016) eff = get_1mu_trigeff2016.GetVector(b)[0];
    else if (b.SampleType()==2017) eff = get_1mu_trigeff2017.GetVector(b)[0];
    else if (b.SampleType()==2018) eff = get_1mu_trigeff2018.GetVector(b)[0];
  } else if (b.nel()==2 && b.nmu()==0) { // 2 lepton MC control sample
    if (b.SampleType()==2016) eff = get_2el_trigeff2016.GetVector(b)[0];
    else if (b.SampleType()==2017) eff = get_2el_trigeff2017.GetVector(b)[0];
    else if (b.SampleType()==2018) eff = get_2el_trigeff2018.GetVector(b)[0];
  } else if (b.nel()==0 && b.nmu()==2) { // 2 lepton MC control sample
    if (b.SampleType()==2016) eff = get_2mu_trigeff2016.GetVector(b)[0];
    else if (b.SampleType()==2017) eff = get_2mu_trigeff2017.GetVector(b)[0];
    else if (b.SampleType()==2018) eff = get_2mu_trigeff2018.GetVector(b)[0];
  }
  return eff;
});

const NamedFunc eff_higtrig_run2_syst_up("eff_higtrig_run2_syst_up", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup, errdown; // Not used, but for reference
  float eff = 1.;
  std::vector<double> trig_vec;
  errup=0;errdown=0;
  errup+=errdown;
  if(b.type()>0 && b.type()<1000) eff = 1; // data

  else if(b.nvlep()==0){ // search MC sample and qcd MC control sample
    if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
      if (b.SampleType()==2016) trig_vec = get_0l_fakemet_trigeff2016.GetVector(b);
      else if (b.SampleType()==2017) trig_vec = get_0l_fakemet_trigeff2017.GetVector(b);
      else if (b.SampleType()==2018) trig_vec = get_0l_fakemet_trigeff2018.GetVector(b);
    } else { // TRUE MET
      if (b.SampleType()==2016) trig_vec = get_0l_trigeff2016.GetVector(b);
      else if (b.SampleType()==2017) trig_vec = get_0l_trigeff2017.GetVector(b);
      else if (b.SampleType()==2018) trig_vec = get_0l_trigeff2018.GetVector(b);
    }
  } else if (b.nel()==1 && b.nmu()==0) { // 1 lepton MC control sample
    if (b.SampleType()==2016) trig_vec = get_1el_trigeff2016.GetVector(b);
    else if (b.SampleType()==2017) trig_vec = get_1el_trigeff2017.GetVector(b);
    else if (b.SampleType()==2018) trig_vec = get_1el_trigeff2018.GetVector(b);
  } else if (b.nel()==0 && b.nmu()==1) { // 1 lepton MC control sample
    if (b.SampleType()==2016) trig_vec = get_1mu_trigeff2016.GetVector(b);
    else if (b.SampleType()==2017) trig_vec = get_1mu_trigeff2017.GetVector(b);
    else if (b.SampleType()==2018) trig_vec = get_1mu_trigeff2018.GetVector(b);
  } else if (b.nel()==2 && b.nmu()==0) { // 2 lepton MC control sample
    if (b.SampleType()==2016) trig_vec = get_2el_trigeff2016.GetVector(b);
    else if (b.SampleType()==2017) trig_vec = get_2el_trigeff2017.GetVector(b);
    else if (b.SampleType()==2018) trig_vec = get_2el_trigeff2018.GetVector(b);
  } else if (b.nel()==0 && b.nmu()==2) { // 2 lepton MC control sample
    if (b.SampleType()==2016) trig_vec = get_2mu_trigeff2016.GetVector(b);
    else if (b.SampleType()==2017) trig_vec = get_2mu_trigeff2017.GetVector(b);
    else if (b.SampleType()==2018) trig_vec = get_2mu_trigeff2018.GetVector(b);
  }
  eff = trig_vec[0]+trig_vec[1];
  return eff;
});

const NamedFunc eff_higtrig_run2_syst_down("eff_higtrig_run2_syst_down", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup, errdown; // Not used, but for reference
  float eff = 1.;
  std::vector<double> trig_vec;
  errup=0;errdown=0;
  errup+=errdown;
  if(b.type()>0 && b.type()<1000) eff = 1; // data

  else if(b.nvlep()==0){ // search MC sample and qcd MC control sample
    if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
      if (b.SampleType()==2016) trig_vec = get_0l_fakemet_trigeff2016.GetVector(b);
      else if (b.SampleType()==2017) trig_vec = get_0l_fakemet_trigeff2017.GetVector(b);
      else if (b.SampleType()==2018) trig_vec = get_0l_fakemet_trigeff2018.GetVector(b);
    } else { // TRUE MET
      if (b.SampleType()==2016) trig_vec = get_0l_trigeff2016.GetVector(b);
      else if (b.SampleType()==2017) trig_vec = get_0l_trigeff2017.GetVector(b);
      else if (b.SampleType()==2018) trig_vec = get_0l_trigeff2018.GetVector(b);
    }
  } else if (b.nel()==1 && b.nmu()==0) { // 1 lepton MC control sample
    if (b.SampleType()==2016) trig_vec = get_1el_trigeff2016.GetVector(b);
    else if (b.SampleType()==2017) trig_vec = get_1el_trigeff2017.GetVector(b);
    else if (b.SampleType()==2018) trig_vec = get_1el_trigeff2018.GetVector(b);
  } else if (b.nel()==0 && b.nmu()==1) { // 1 lepton MC control sample
    if (b.SampleType()==2016) trig_vec = get_1mu_trigeff2016.GetVector(b);
    else if (b.SampleType()==2017) trig_vec = get_1mu_trigeff2017.GetVector(b);
    else if (b.SampleType()==2018) trig_vec = get_1mu_trigeff2018.GetVector(b);
  } else if (b.nel()==2 && b.nmu()==0) { // 2 lepton MC control sample
    if (b.SampleType()==2016) trig_vec = get_2el_trigeff2016.GetVector(b);
    else if (b.SampleType()==2017) trig_vec = get_2el_trigeff2017.GetVector(b);
    else if (b.SampleType()==2018) trig_vec = get_2el_trigeff2018.GetVector(b);
  } else if (b.nel()==0 && b.nmu()==2) { // 2 lepton MC control sample
    if (b.SampleType()==2016) trig_vec = get_2mu_trigeff2016.GetVector(b);
    else if (b.SampleType()==2017) trig_vec = get_2mu_trigeff2017.GetVector(b);
    else if (b.SampleType()==2018) trig_vec = get_2mu_trigeff2018.GetVector(b);
  }
  eff = trig_vec[0]+trig_vec[2];
  return eff;
});

const NamedFunc eff_higtrig_run2_v0("eff_higtrig_run2_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup, errdown; // Not used, but for reference
  float eff = 1.;
  errup=0;errdown=0;
  errup+=errdown;
  if(b.type()>0 && b.type()<1000) eff = 1; // data

  else if(b.nvlep()==0){ // search MC sample and qcd MC control sample
    if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
      if (b.SampleType()==2016) eff = get_0l_fakemet_trigeff2016_v0.GetScalar(b);
      else if (b.SampleType()==2017) eff = get_0l_fakemet_trigeff2017_v0.GetScalar(b);
      else if (b.SampleType()==2018) eff = get_0l_fakemet_trigeff2018_v0.GetScalar(b);
    } else { // TRUE MET
      if (b.SampleType()==2016) eff = get_0l_trigeff2016_v0.GetScalar(b);
      else if (b.SampleType()==2017) eff = get_0l_trigeff2017_v0.GetScalar(b);
      else if (b.SampleType()==2018) eff = get_0l_trigeff2018_v0.GetScalar(b);
    }
  } else if (b.nel()==1 && b.nmu()==0) { // 1 lepton MC control sample
    if (b.SampleType()==2016) eff = get_1el_trigeff2016_v0.GetScalar(b);
    else if (b.SampleType()==2017) eff = get_1el_trigeff2017_v0.GetScalar(b);
    else if (b.SampleType()==2018) eff = get_1el_trigeff2018_v0.GetScalar(b);
  } else if (b.nel()==0 && b.nmu()==1) { // 1 lepton MC control sample
    if (b.SampleType()==2016) eff = get_1mu_trigeff2016_v0.GetScalar(b);
    else if (b.SampleType()==2017) eff = get_1mu_trigeff2017_v0.GetScalar(b);
    else if (b.SampleType()==2018) eff = get_1mu_trigeff2018_v0.GetScalar(b);
  } else if (b.nel()==2 && b.nmu()==0) { // 2 lepton MC control sample
    if (b.SampleType()==2016) eff = get_2el_trigeff2016_v0.GetScalar(b);
    else if (b.SampleType()==2017) eff = get_2el_trigeff2017_v0.GetScalar(b);
    else if (b.SampleType()==2018) eff = get_2el_trigeff2018_v0.GetScalar(b);
  } else if (b.nel()==0 && b.nmu()==2) { // 2 lepton MC control sample
    if (b.SampleType()==2016) eff = get_2mu_trigeff2016_v0.GetScalar(b);
    else if (b.SampleType()==2017) eff = get_2mu_trigeff2017_v0.GetScalar(b);
    else if (b.SampleType()==2018) eff = get_2mu_trigeff2018_v0.GetScalar(b);
  }
  return eff;
});
  
////// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
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
   } else if(b.nel()==1 && b.nmu()==0){
     vector<float> lep_pt; 
     if (b.lep_pt()->size()>0) lep_pt.push_back(b.lep_pt()->at(0));
     else lep_pt.push_back(0);
     if(lep_pt[0]>  20 && lep_pt[0]<=  25 && met>=0 && met<= 110) {eff=0.160; errup=0.019; errdown = 0.017;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met>=0 && met<= 110) {eff=0.400; errup=0.024; errdown=0.024;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met>=0 && met<= 110) {eff=0.728; errup=0.006; errdown=0.006;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met>=0 && met<= 110) {eff=0.880; errup=0.017; errdown=0.019;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met>=0 && met<= 110) {eff=0.950; errup=0.003; errdown=0.003;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met>110 && met<= 120) {eff=0.244; errup=0.024; errdown=0.023;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met>110 && met<= 120) {eff=0.420; errup=0.027; errdown=0.027;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met>110 && met<= 120) {eff=0.761; errup=0.007; errdown=0.007;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met>110 && met<= 120) {eff=0.918; errup=0.015; errdown=0.017;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met>110 && met<= 120) {eff=0.958; errup=0.003; errdown=0.003;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met>120 && met<= 130) {eff=0.331; errup=0.030; errdown=0.029;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met>120 && met<= 130) {eff=0.500; errup=0.031; errdown=0.031;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met>120 && met<= 130) {eff=0.800; errup=0.007; errdown=0.007;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met>120 && met<= 130) {eff=0.928; errup=0.015; errdown=0.018;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met>120 && met<= 130) {eff=0.960; errup=0.003; errdown=0.003;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met>130 && met<= 140) {eff=0.491; errup=0.031; errdown=0.031;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met>130 && met<= 140) {eff=0.608; errup=0.033; errdown=0.034;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met>130 && met<= 140) {eff=0.831; errup=0.007; errdown=0.007;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met>130 && met<= 140) {eff=0.931; errup=0.016; errdown=0.020;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met>130 && met<= 140) {eff=0.967; errup=0.003; errdown=0.003;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met>140 && met<= 150) {eff=0.573; errup=0.033; errdown=0.033;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met>140 && met<= 150) {eff=0.677; errup=0.035; errdown=0.037;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met>140 && met<= 150) {eff=0.856; errup=0.007; errdown=0.007;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met>140 && met<= 150) {eff=0.923; errup=0.018; errdown=0.022;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met>140 && met<= 150) {eff=0.971; errup=0.003; errdown=0.004;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met>150 && met<= 160) {eff=0.643; errup=0.037; errdown=0.039;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met> 150 && met<= 160) {eff=0.738; errup=0.033; errdown=0.036;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met> 150 && met<= 160) {eff=0.871; errup=0.007; errdown=0.007;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met> 150 && met<= 160) {eff=0.935; errup=0.016; errdown=0.021;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met> 150 && met<= 160) {eff=0.982; errup=0.003; errdown=0.003;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met> 160 && met<= 170) {eff=0.760; errup=0.034; errdown=0.038;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met> 160 && met<= 170) {eff=0.792; errup=0.034; errdown=0.038;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met> 160 && met<= 170) {eff=0.910; errup=0.007; errdown=0.007;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met> 160 && met<= 170) {eff=0.970; errup=0.013; errdown=0.020;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met> 160 && met<= 170) {eff=0.981; errup=0.003; errdown=0.004;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met> 170 && met<= 180) {eff=0.829; errup=0.033; errdown=0.038;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met> 170 && met<= 180) {eff=0.863; errup=0.031; errdown=0.037;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met> 170 && met<= 180) {eff=0.937; errup=0.006; errdown=0.006;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met> 170 && met<= 180) {eff=0.988; errup=0.008; errdown=0.016;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met> 170 && met<= 180) {eff=0.982; errup=0.003; errdown=0.004;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met> 180 && met<= 190) {eff=0.761; errup=0.041; errdown=0.046;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met> 180 && met<= 190) {eff=0.863; errup=0.032; errdown=0.038;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met> 180 && met<= 190) {eff=0.939; errup=0.006; errdown=0.007;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met> 180 && met<= 190) {eff=0.969; errup=0.015; errdown=0.024;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met> 180 && met<= 190) {eff=0.984; errup=0.003; errdown=0.004;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met> 190 && met<= 200) {eff=0.892; errup=0.033; errdown=0.042;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met> 190 && met<= 200) {eff=0.902; errup=0.030; errdown=0.039;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met> 190 && met<= 200) {eff=0.956; errup=0.006; errdown=0.006;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met> 190 && met<= 200) {eff=0.983; errup=0.011; errdown=0.022;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met> 190 && met<= 200) {eff=0.993; errup=0.002; errdown=0.003;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met> 200 && met<= 210) {eff=0.950; errup=0.021; errdown=0.032;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met> 200 && met<= 210) {eff=0.951; errup=0.023; errdown=0.037;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met> 200 && met<= 210) {eff=0.973; errup=0.005; errdown=0.005;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met> 200 && met<= 210) {eff=1.000; errup=0.000; errdown=0.018;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met> 200 && met<= 210) {eff=0.985; errup=0.003; errdown=0.004;}
     else if(lep_pt[0]> 20 && lep_pt[0]<=  25 && met> 210 && met<=9999) {eff=0.974; errup=0.005; errdown=0.006;}
     else if(lep_pt[0]> 25 && lep_pt[0]<=  30 && met> 210 && met<=9999) {eff=0.981; errup=0.004; errdown=0.005;}
     else if(lep_pt[0]> 30 && lep_pt[0]<= 110 && met> 210 && met<=9999) {eff=0.992; errup=0.001; errdown=0.001;}
     else if(lep_pt[0]>110 && lep_pt[0]<= 120 && met> 210 && met<=9999) {eff=0.997; errup=0.002; errdown=0.003;}
     else if(lep_pt[0]>120 && lep_pt[0]<=9999 && met> 210 && met<=9999) {eff=0.996; errup=0.001; errdown=0.001;}
  
     //// MET || Mu24 || Mu50
     //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]||trig[19]||trig[55]||trig[21])"
   } else if(b.nel()==0 && b.nmu()==1){
     vector<float> lep_pt; 
     if (b.lep_pt()->size()>0) lep_pt.push_back(b.lep_pt()->at(0));
     else lep_pt.push_back(0);
     if(lep_pt[0]>  20 && lep_pt[0]<=  25 && met>= 0 && met<= 110) {eff=0.271; errup=0.017; errdown = 0.016;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met>= 0 && met<= 110) {eff=0.725; errup=0.017; errdown=0.018;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met>= 0 && met<= 110) {eff=0.814; errup=0.008; errdown=0.009;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met>= 0 && met<= 110) {eff=0.964; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 110 && met<= 120) {eff=0.363; errup=0.020; errdown=0.020;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 110 && met<= 120) {eff=0.755; errup=0.018; errdown=0.019;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 110 && met<= 120) {eff=0.842; errup=0.009; errdown=0.009;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 110 && met<= 120) {eff=0.969; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 120 && met<= 130) {eff=0.452; errup=0.022; errdown=0.022;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 120 && met<= 130) {eff=0.824; errup=0.018; errdown=0.019;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 120 && met<= 130) {eff=0.869; errup=0.009; errdown=0.009;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 120 && met<= 130) {eff=0.971; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 130 && met<= 140) {eff=0.590; errup=0.025; errdown=0.025;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 130 && met<= 140) {eff=0.875; errup=0.017; errdown=0.019;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 130 && met<= 140) {eff=0.904; errup=0.008; errdown=0.009;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 130 && met<= 140) {eff=0.972; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 140 && met<= 150) {eff=0.660; errup=0.026; errdown=0.027;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 140 && met<= 150) {eff=0.891; errup=0.017; errdown=0.019;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 140 && met<= 150) {eff=0.938; errup=0.007; errdown=0.008;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 140 && met<= 150) {eff=0.980; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 150 && met<= 160) {eff=0.778; errup=0.024; errdown=0.026;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 150 && met<= 160) {eff=0.915; errup=0.016; errdown=0.019;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 150 && met<= 160) {eff=0.940; errup=0.008; errdown=0.009;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 150 && met<= 160) {eff=0.984; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 160 && met<= 170) {eff=0.798; errup=0.026; errdown=0.029;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 160 && met<= 170) {eff=0.946; errup=0.015; errdown=0.020;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 160 && met<= 170) {eff=0.967; errup=0.006; errdown=0.008;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 160 && met<= 170) {eff=0.991; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 170 && met<= 180) {eff=0.885; errup=0.022; errdown=0.025;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 170 && met<= 180) {eff=0.937; errup=0.016; errdown=0.021;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 170 && met<= 180) {eff=0.977; errup=0.006; errdown=0.007;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 170 && met<= 180) {eff=0.987; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 180 && met<= 190) {eff=0.927; errup=0.019; errdown=0.024;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 180 && met<= 190) {eff=0.958; errup=0.014; errdown=0.019;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 180 && met<= 190) {eff=0.974; errup=0.006; errdown=0.008;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 180 && met<= 190) {eff=0.992; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 190 && met<= 200) {eff=0.921; errup=0.019; errdown=0.024;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 190 && met<= 200) {eff=0.965; errup=0.014; errdown=0.020;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 190 && met<= 200) {eff=0.991; errup=0.004; errdown=0.006;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 190 && met<= 200) {eff=0.991; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 200 && met<= 210) {eff=0.926; errup=0.022; errdown=0.028;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 200 && met<= 210) {eff=0.994; errup=0.005; errdown=0.015;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 200 && met<= 210) {eff=0.994; errup=0.003; errdown=0.006;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 200 && met<= 210) {eff=0.994; errup=0.002; errdown=0.002;}
     else if(lep_pt[0]>20 && lep_pt[0]<=  25 && met> 210 && met<=9999) {eff=0.981; errup=0.004; errdown=0.004;}
     else if(lep_pt[0]>25 && lep_pt[0]<=  30 && met> 210 && met<=9999) {eff=0.994; errup=0.002; errdown=0.003;}
     else if(lep_pt[0]>30 && lep_pt[0]<=  50 && met> 210 && met<=9999) {eff=0.996; errup=0.001; errdown=0.001;}
     else if(lep_pt[0]>50 && lep_pt[0]<=9999 && met> 210 && met<=9999) {eff=0.997; errup=0.000; errdown=0.000;}
  
     //// Ele27 || Ele105 || Ele115
     //// "(trig[22]||trig[40]||trig[24]||trig[41])"
   } else if(b.nel()==2 && b.nmu()==0){
     vector<float> lep_pt; 
     if (b.lep_pt()->size()>0) lep_pt.push_back(b.lep_pt()->at(0));
     else lep_pt.push_back(0);
     if(lep_pt[0]>  40 && lep_pt[0]<=  45) {eff = 0.944; errup = 0.015; errdown = 0.019;}
     else if(lep_pt[0]>  45 && lep_pt[0]<=  50) {eff = 0.910; errup = 0.015; errdown = 0.017;}
     else if(lep_pt[0]>  50 && lep_pt[0]<=  55) {eff = 0.927; errup = 0.013; errdown = 0.015;}
     else if(lep_pt[0]>  55 && lep_pt[0]<=  60) {eff = 0.912; errup = 0.013; errdown = 0.015;}
     else if(lep_pt[0]>  60 && lep_pt[0]<=  65) {eff = 0.941; errup = 0.011; errdown = 0.013;}
     else if(lep_pt[0]>  65 && lep_pt[0]<=  70) {eff = 0.901; errup = 0.014; errdown = 0.016;}
     else if(lep_pt[0]>  70 && lep_pt[0]<=  75) {eff = 0.921; errup = 0.013; errdown = 0.016;}
     else if(lep_pt[0]>  75 && lep_pt[0]<=  80) {eff = 0.947; errup = 0.011; errdown = 0.014;}
     else if(lep_pt[0]>  80 && lep_pt[0]<=  85) {eff = 0.954; errup = 0.011; errdown = 0.013;}
     else if(lep_pt[0]>  85 && lep_pt[0]<=  90) {eff = 0.939; errup = 0.012; errdown = 0.014;}
     else if(lep_pt[0]>  90 && lep_pt[0]<=  95) {eff = 0.940; errup = 0.012; errdown = 0.015;}
     else if(lep_pt[0]>  95 && lep_pt[0]<= 100) {eff = 0.932; errup = 0.014; errdown = 0.017;}
     else if(lep_pt[0]> 100 && lep_pt[0]<= 105) {eff = 0.934; errup = 0.014; errdown = 0.017;}
     else if(lep_pt[0]> 105 && lep_pt[0]<= 110) {eff = 0.965; errup = 0.010; errdown = 0.014;}
     else if(lep_pt[0]> 110 && lep_pt[0]<=9999) {eff = 0.994; errup = 0.001; errdown = 0.001;}
  
     //// Mu24 || Mu50
     //// "(trig[19]||trig[55]||trig[21])"
   } else if(b.nel()==0 && b.nmu()==2){
     vector<float> lep_pt; 
     if (b.lep_pt()->size()>0) lep_pt.push_back(b.lep_pt()->at(0));
     else lep_pt.push_back(0);
     if(lep_pt[0]>  40 && lep_pt[0]<=  45) {eff = 0.959; errup = 0.010; errdown = 0.012;}
     if(lep_pt[0]>  45 && lep_pt[0]<=  50) {eff = 0.970; errup = 0.006; errdown = 0.007;}
     if(lep_pt[0]>  50 && lep_pt[0]<=9999) {eff = 0.982; errup = 0.000; errdown = 0.000;}
  }
  return eff;
 });

const NamedFunc wgt_syst_lowptleptrig("wgt_syst_lowptleptrig", [](const Baby &b) -> NamedFunc::ScalarType{
    //for 2017
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 250 && met> 150 && met<= 155) {eff = 0.165414; errup = 0.0386864; errdown = 0.0330951;}
  else if (ht> 250 && ht<= 350 && met> 150 && met<= 155) {eff = 0.307692; errup = 0.0341672; errdown = 0.0323054;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.403614; errup = 0.0415856; errdown = 0.0403373;}
  else if (ht> 450 && ht<= 600 && met> 150 && met<= 155) {eff = 0.393443; errup = 0.0337166; errdown = 0.0327883;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.398734; errup = 0.0426681; errdown = 0.0412877;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.46875; errup = 0.0701234; errdown = 0.0690358;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.647059; errup = 0.129922; errdown = 0.150722;}
  else if (ht> 0 && ht<= 250 && met> 155 && met<= 160) {eff = 0.177966; errup = 0.042371; errdown = 0.0362895;}
  else if (ht> 250 && ht<= 350 && met> 155 && met<= 160) {eff = 0.376682; errup = 0.0352; errdown = 0.0340209;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.403974; errup = 0.0437778; errdown = 0.0424063;}
  else if (ht> 450 && ht<= 600 && met> 155 && met<= 160) {eff = 0.476596; errup = 0.0347357; errdown = 0.0345241;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.48366; errup = 0.0436303; errdown = 0.0434004;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.581818; errup = 0.0731686; errdown = 0.0765096;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.263158; errup = 0.13882; errdown = 0.108565;}
  else if (ht> 0 && ht<= 250 && met> 160 && met<= 165) {eff = 0.22093; errup = 0.0539358; errdown = 0.0466584;}
  else if (ht> 250 && ht<= 350 && met> 160 && met<= 165) {eff = 0.470874; errup = 0.0372593; errdown = 0.0369577;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.588652; errup = 0.0441377; errdown = 0.0454966;}
  else if (ht> 450 && ht<= 600 && met> 160 && met<= 165) {eff = 0.509615; errup = 0.0369219; errdown = 0.0370205;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.586466; errup = 0.0455717; errdown = 0.0469796;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.439024; errup = 0.0901517; errdown = 0.0867656;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.45; errup = 0.135158; errdown = 0.129246;}
  else if (ht> 0 && ht<= 250 && met> 165 && met<= 170) {eff = 0.36; errup = 0.0637251; errdown = 0.0595793;}
  else if (ht> 250 && ht<= 350 && met> 165 && met<= 170) {eff = 0.476471; errup = 0.0412653; errdown = 0.0409684;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.578947; errup = 0.0496353; errdown = 0.0511431;}
  else if (ht> 450 && ht<= 600 && met> 165 && met<= 170) {eff = 0.59009; errup = 0.0347456; errdown = 0.0356101;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.548872; errup = 0.0463198; errdown = 0.0471149;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.638889; errup = 0.0880506; errdown = 0.0969352;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.454545; errup = 0.189662; errdown = 0.179582;}
  else if (ht> 0 && ht<= 250 && met> 170 && met<= 175) {eff = 0.275862; errup = 0.0711112; errdown = 0.0623488;}
  else if (ht> 250 && ht<= 350 && met> 170 && met<= 175) {eff = 0.5; errup = 0.0439995; errdown = 0.0439995;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.59; errup = 0.0529068; errdown = 0.0548769;}
  else if (ht> 450 && ht<= 600 && met> 170 && met<= 175) {eff = 0.60221; errup = 0.0384238; errdown = 0.0396347;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.596491; errup = 0.0491774; errdown = 0.0510213;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.645833; errup = 0.0751004; errdown = 0.0819963;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.818182; errup = 0.116508; errdown = 0.191402;}
  else if (ht> 0 && ht<= 250 && met> 175 && met<= 180) {eff = 0.392857; errup = 0.112356; errdown = 0.103445;}
  else if (ht> 250 && ht<= 350 && met> 175 && met<= 180) {eff = 0.6; errup = 0.0448856; errdown = 0.0464896;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.7; errup = 0.0511019; errdown = 0.0560215;}
  else if (ht> 450 && ht<= 600 && met> 175 && met<= 180) {eff = 0.678363; errup = 0.0374083; errdown = 0.0396568;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.66055; errup = 0.0480947; errdown = 0.0513192;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.733333; errup = 0.0702056; errdown = 0.082142;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.692308; errup = 0.140801; errdown = 0.177171;}
  else if (ht> 0 && ht<= 250 && met> 180 && met<= 185) {eff = 0.421053; errup = 0.0942357; errdown = 0.0894837;}
  else if (ht> 250 && ht<= 350 && met> 180 && met<= 185) {eff = 0.601942; errup = 0.051728; errdown = 0.0538929;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.670588; errup = 0.0542946; errdown = 0.0587372;}
  else if (ht> 450 && ht<= 600 && met> 180 && met<= 185) {eff = 0.751592; errup = 0.0358197; errdown = 0.0393064;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.68; errup = 0.0494059; errdown = 0.0533658;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.59375; errup = 0.0972359; errdown = 0.104003;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.8; errup = 0.128046; errdown = 0.205454;}
  else if (ht> 0 && ht<= 250 && met> 185 && met<= 190) {eff = 0.272727; errup = 0.0992802; errdown = 0.0831565;}
  else if (ht> 250 && ht<= 350 && met> 185 && met<= 190) {eff = 0.674699; errup = 0.0547582; errdown = 0.0594239;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.810811; errup = 0.0472278; errdown = 0.0567879;}
  else if (ht> 450 && ht<= 600 && met> 185 && met<= 190) {eff = 0.69281; errup = 0.0390696; errdown = 0.0417994;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.77907; errup = 0.0466584; errdown = 0.0539358;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.827586; errup = 0.0723692; errdown = 0.0998191;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.785714; errup = 0.1142; errdown = 0.165474;}
  else if (ht> 0 && ht<= 250 && met> 190 && met<= 195) {eff = 0.478261; errup = 0.123779; errdown = 0.121562;}
  else if (ht> 250 && ht<= 350 && met> 190 && met<= 195) {eff = 0.658824; errup = 0.0548887; errdown = 0.0590212;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.806452; errup = 0.0521796; errdown = 0.0635295;}
  else if (ht> 450 && ht<= 600 && met> 190 && met<= 195) {eff = 0.798387; errup = 0.0372733; errdown = 0.0425968;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.771084; errup = 0.0481984; errdown = 0.0555253;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.575758; errup = 0.0967863; errdown = 0.102076;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.5; errup = 0.220457; errdown = 0.220457;}
  else if (ht> 0 && ht<= 250 && met> 195 && met<= 200) {eff = 0.695652; errup = 0.10473; errdown = 0.124968;}
  else if (ht> 250 && ht<= 350 && met> 195 && met<= 200) {eff = 0.724138; errup = 0.050503; errdown = 0.0562302;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.87013; errup = 0.0390847; errdown = 0.0502291;}
  else if (ht> 450 && ht<= 600 && met> 195 && met<= 200) {eff = 0.877358; errup = 0.0324472; errdown = 0.0405703;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.831461; errup = 0.0409024; errdown = 0.0493424;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.882353; errup = 0.0756114; errdown = 0.134601;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.857143; errup = 0.0917089; errdown = 0.158271;}
  else if (ht> 0 && ht<= 250 && met> 200 && met<= 210) {eff = 0.772727; errup = 0.0944237; errdown = 0.12452;}
  else if (ht> 250 && ht<= 350 && met> 200 && met<= 210) {eff = 0.824074; errup = 0.0377742; errdown = 0.0444931;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.895652; errup = 0.028883; errdown = 0.0367676;}
  else if (ht> 450 && ht<= 600 && met> 200 && met<= 210) {eff = 0.876652; errup = 0.0221898; errdown = 0.0258332;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.891566; errup = 0.024485; errdown = 0.0297711;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.878049; errup = 0.0516458; errdown = 0.0740833;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.842105; errup = 0.0846498; errdown = 0.130121;}
  else if (ht> 0 && ht<= 250 && met> 210 && met<= 220) {eff = 0.785714; errup = 0.1142; errdown = 0.165474;}
  else if (ht> 250 && ht<= 350 && met> 210 && met<= 220) {eff = 0.8875; errup = 0.035833; errdown = 0.0471369;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.906977; errup = 0.0315422; errdown = 0.0426732;}
  else if (ht> 450 && ht<= 600 && met> 210 && met<= 220) {eff = 0.956522; errup = 0.0149067; errdown = 0.0207436;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.910569; errup = 0.0259506; errdown = 0.0336394;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 1; errup = 0; errdown = 0.0439098;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.928571; errup = 0.0591648; errdown = 0.145681;}
  else if (ht> 0 && ht<= 250 && met> 220 && met<= 230) {eff = 0.875; errup = 0.103637; errdown = 0.23225;}
  else if (ht> 250 && ht<= 350 && met> 220 && met<= 230) {eff = 0.892308; errup = 0.0388522; errdown = 0.0531815;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.947368; errup = 0.02502; errdown = 0.0396606;}
  else if (ht> 450 && ht<= 600 && met> 220 && met<= 230) {eff = 0.946746; errup = 0.0172067; errdown = 0.0233743;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.932203; errup = 0.0231199; errdown = 0.0317357;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.95; errup = 0.0322293; errdown = 0.0621707;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.944444; errup = 0.046004; errdown = 0.116415;}
  else if (ht> 0 && ht<= 250 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 250 && ht<= 350 && met> 230 && met<= 240) {eff = 0.972973; errup = 0.0223689; errdown = 0.0594217;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.0252456;}
  else if (ht> 450 && ht<= 600 && met> 230 && met<= 240) {eff = 0.941176; errup = 0.0201; errdown = 0.0277309;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.936709; errup = 0.0270672; errdown = 0.0405476;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.952381; errup = 0.0306976; errdown = 0.0593794;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.15411;}
  else if (ht> 0 && ht<= 250 && met> 240 && met<= 250) {eff = 0.5; errup = 0.314699; errdown = 0.314699;}
  else if (ht> 250 && ht<= 350 && met> 240 && met<= 250) {eff = 0.94; errup = 0.0324767; errdown = 0.054936;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.977778; errup = 0.0183906; errdown = 0.0492514;}
  else if (ht> 450 && ht<= 600 && met> 240 && met<= 250) {eff = 0.966102; errup = 0.0161538; errdown = 0.0259913;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.975904; errup = 0.0155489; errdown = 0.0308989;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.0879414;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.168149;}
  else if (ht> 0 && ht<= 250 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 250 && ht<= 350 && met> 250 && met<= 275) {eff = 0.960784; errup = 0.0252892; errdown = 0.0493939;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 275) {eff = 0.980392; errup = 0.0126548; errdown = 0.0252763;}
  else if (ht> 450 && ht<= 600 && met> 250 && met<= 275) {eff = 0.969027; errup = 0.0113532; errdown = 0.0162848;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.978261; errup = 0.0103756; errdown = 0.0168548;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.0317826;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.108691;}
  else if (ht> 0 && ht<= 250 && met> 275 && met<= 300) {eff = 3.08544e-320; errup = 0; errdown = 0;}
  else if (ht> 250 && ht<= 350 && met> 275 && met<= 300) {eff = 0.933333; errup = 0.0552158; errdown = 0.13708;}
  else if (ht> 350 && ht<= 450 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.0323409;}
  else if (ht> 450 && ht<= 600 && met> 275 && met<= 300) {eff = 0.983146; errup = 0.00915875; errdown = 0.0161224;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.0153517;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.0576587;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.142229;}
  else if (ht> 0 && ht<= 250 && met> 300 && met<= 9999) {eff = 2.53456e-321; errup = -1; errdown = -1;}
  else if (ht> 250 && ht<= 350 && met> 300 && met<= 9999) {eff = 0.933333; errup = 0.0552158; errdown = 0.13708;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0312433;}
  else if (ht> 450 && ht<= 600 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0100098;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.993506; errup = 0.00537236; errdown = 0.0147726;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0354547;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0683597;}
  return eff/get_0l_trigeff2017.GetVector(b)[0];
});

const NamedFunc weight_higd("weight_hig_deep",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type()==-999999){ //normalize weights for TChiHZ benchmarks
    if (b.mprod()==225) return b.weight()/b.w_btag()*b.w_bhig()/.9666;
    else if (b.mprod()==400) return b.weight()/b.w_btag()*b.w_bhig()/.9705;
    else if (b.mprod()==700) return b.weight()/b.w_btag()*b.w_bhig()/.987;
  } else if (b.type()>0 && b.type()<1000) return 1;
  
  return b.weight()/b.w_btag()*b.w_bhig();
});

const NamedFunc pico_weight_higd("pico_weight_hig_deep",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type()==-999999){ //normalize weights for TChiHZ benchmarks
    if (b.mprod()==225) return b.weight()/b.w_btag()*b.w_bhig()/.9666;
    else if (b.mprod()==400) return b.weight()/b.w_btag()*b.w_bhig()/.9705;
    else if (b.mprod()==700) return b.weight()/b.w_btag()*b.w_bhig()/.987;
  } else if (b.type()>0 && b.type()<1000) return 1;
  
  return b.weight()/b.w_btag()*b.w_bhig();
});


const NamedFunc mhig("mhig",[](const Baby &b) -> NamedFunc::ScalarType{
  float mass = -999;
  for (unsigned i(0); i<b.mc_mass()->size(); i++){
    if (b.mc_id()->at(i)==1000023) {
      mass = b.mc_mass()->at(i);
      break;
    }
  }
  return mass;
});

const NamedFunc nb_exci("nb_exci",[](const Baby &b) -> NamedFunc::ScalarType{
  int nb=0;
  for (unsigned i(0); i<b.mc_id()->size(); i++){
    if (abs(b.mc_id()->at(i))==5 && b.mc_status()->at(i)==23 
	&& abs(b.mc_mom()->at(i))!=6 && abs(b.mc_mom()->at(i))!=24 && abs(b.mc_mom()->at(i))!=23) 
      nb++;
  }
  return nb;
});

const NamedFunc nb_gs("nb_gs",[](const Baby &b) -> NamedFunc::ScalarType{
  int nb=0;
  for (unsigned i(0); i<b.mc_id()->size(); i++){
    if (abs(b.mc_id()->at(i))==5 && b.mc_status()->at(i)!=23 ) 
      nb++;
  }
  return nb;
});


// Reconstruct higgs
Double_t DeltaR(const ROOT::Math::PtEtaPhiMVector & v1, const ROOT::Math::PtEtaPhiMVector & v2)
{
   Double_t deta = v1.Eta()-v2.Eta();
   Double_t dphi = TVector2::Phi_mpi_pi(v1.Phi()-v2.Phi());
   return TMath::Sqrt( deta*deta+dphi*dphi );
}

bool greater_btag(pair<float, unsigned> bjet_1, pair<float, unsigned> bjet_2) {
  return bjet_1.first > bjet_2.first;
}

bool smaller_dm(tuple<double, double, vector<unsigned> > higgs_1, tuple<double, double, vector<unsigned> > higgs_2) {
  return get<0>(higgs_1) < get<0>(higgs_2);
}

vector<unsigned> get_higgs_bbjet_indices(vector<float> const & jet_m, vector<float> const & jet_deepcsv, vector<float> const & jet_pt, vector<float> const & jet_eta, vector<float> const & jet_phi, vector<bool> const & jet_isgood) {
  //cout<<"event start"<<endl;
  // Reconstruct resolved
  // Sort by btag
  vector<pair<float, unsigned> > bjet;
  for (unsigned ijet = 0; ijet < jet_m.size(); ijet++) {
    // Filter jets
    if (jet_pt[ijet]<=30) continue;
    if (fabs(jet_eta[ijet])>2.4) continue;
    if (!jet_isgood[ijet]) continue;
    //cout<<"ijet ["<<ijet<<"] btag: "<<jet_deepcsv[ijet]<<" pt:"<<jet_pt[ijet]<<endl;
    bjet.push_back({jet_deepcsv[ijet],ijet});
  }
  sort(bjet.begin(), bjet.end(), greater_btag);
  //for (unsigned ibjet = 0; ibjet < bjet.size(); ibjet++) {
  //  cout<<"ibjet ["<<ibjet<<"] btag: "<<bjet[ibjet].first<<" jet index: "<<bjet[ibjet].second<<endl;
  //}
  // Make three higgs candidates
  vector<vector<unsigned> > higgs_bjet_indices = {{0,1,2,3}, {0,2,1,3}, {0,3,1,2}};
  // higgs_info = [higgs_dm, higgs_am, higgs_jet_indices]
  vector<tuple<double, double, vector<unsigned> > > higgs_info;
  for (unsigned ihiggs = 0; ihiggs < higgs_bjet_indices.size(); ihiggs++) {
    //cout<<"ihiggs ["<<ihiggs<<"] bjet index: "<<higgs_bjet_indices[ihiggs][0]<<" "<<higgs_bjet_indices[ihiggs][1]<<" "<<higgs_bjet_indices[ihiggs][2]<<" "<<higgs_bjet_indices[ihiggs][3]<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] jet index: "<<bjet[higgs_bjet_indices[ihiggs][0]].second<<" "<<bjet[higgs_bjet_indices[ihiggs][1]].second<<" "<<bjet[higgs_bjet_indices[ihiggs][2]].second<<" "<<bjet[higgs_bjet_indices[ihiggs][3]].second<<endl;
    ROOT::Math::PtEtaPhiMVector higgs1_b1 (jet_pt[bjet[higgs_bjet_indices[ihiggs][0]].second], jet_eta[bjet[higgs_bjet_indices[ihiggs][0]].second], jet_phi[bjet[higgs_bjet_indices[ihiggs][0]].second], jet_m[bjet[higgs_bjet_indices[ihiggs][0]].second]);
    ROOT::Math::PtEtaPhiMVector higgs1_b2 (jet_pt[bjet[higgs_bjet_indices[ihiggs][1]].second], jet_eta[bjet[higgs_bjet_indices[ihiggs][1]].second], jet_phi[bjet[higgs_bjet_indices[ihiggs][1]].second], jet_m[bjet[higgs_bjet_indices[ihiggs][1]].second]);
    ROOT::Math::PtEtaPhiMVector higgs2_b1 (jet_pt[bjet[higgs_bjet_indices[ihiggs][2]].second], jet_eta[bjet[higgs_bjet_indices[ihiggs][2]].second], jet_phi[bjet[higgs_bjet_indices[ihiggs][2]].second], jet_m[bjet[higgs_bjet_indices[ihiggs][2]].second]);
    ROOT::Math::PtEtaPhiMVector higgs2_b2 (jet_pt[bjet[higgs_bjet_indices[ihiggs][3]].second], jet_eta[bjet[higgs_bjet_indices[ihiggs][3]].second], jet_phi[bjet[higgs_bjet_indices[ihiggs][3]].second], jet_m[bjet[higgs_bjet_indices[ihiggs][3]].second]);
    //cout<<"ihiggs ["<<ihiggs<<"] higg1_b1 e: "<<higgs1_b1.E()<<" pt: "<<higgs1_b1.Pt()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] higg1_b2 e: "<<higgs1_b2.E()<<" pt: "<<higgs1_b2.Pt()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] higg2_b1 e: "<<higgs2_b1.E()<<" pt: "<<higgs2_b1.Pt()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] higg2_b2 e: "<<higgs2_b2.E()<<" pt: "<<higgs2_b2.Pt()<<endl;

    ROOT::Math::PtEtaPhiMVector higgs1 = higgs1_b1 + higgs1_b2;
    ROOT::Math::PtEtaPhiMVector higgs2 = higgs2_b1 + higgs2_b2;
    //cout<<"ihiggs ["<<ihiggs<<"] higg1: "<<higgs1.E()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] higg2: "<<higgs2.E()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] am: "<<(higgs1.M() + higgs2.M())/2<<" dm: "<<fabs(higgs1.M()-higgs2.M())<<endl;

    vector<unsigned> higgs_indices;
    // Sort by higgs pt
    if (higgs1.Pt()>higgs2.Pt()) {
      higgs_indices = {bjet[higgs_bjet_indices[ihiggs][0]].second, bjet[higgs_bjet_indices[ihiggs][1]].second, bjet[higgs_bjet_indices[ihiggs][2]].second, bjet[higgs_bjet_indices[ihiggs][3]].second};
    } else {
      higgs_indices = {bjet[higgs_bjet_indices[ihiggs][2]].second, bjet[higgs_bjet_indices[ihiggs][3]].second, bjet[higgs_bjet_indices[ihiggs][0]].second, bjet[higgs_bjet_indices[ihiggs][1]].second};
    }
    //// Sort by dr
    //if (DeltaR(higgs1_b1, higgs1_b2) < DeltaR(higgs2_b1, higgs2_b2)) {
    //  higgs_indices = {bjet[higgs_bjet_indices[ihiggs][0]].second, bjet[higgs_bjet_indices[ihiggs][1]].second, bjet[higgs_bjet_indices[ihiggs][2]].second, bjet[higgs_bjet_indices[ihiggs][3]].second};
    //} else {
    //  higgs_indices = {bjet[higgs_bjet_indices[ihiggs][2]].second, bjet[higgs_bjet_indices[ihiggs][3]].second, bjet[higgs_bjet_indices[ihiggs][0]].second, bjet[higgs_bjet_indices[ihiggs][1]].second};
    //}

    higgs_info.push_back({fabs(higgs1.M()-higgs2.M()), (higgs1.M() + higgs2.M())/2, higgs_indices});
  }
  // Sort by dm
  sort(higgs_info.begin(), higgs_info.end(), smaller_dm);
  //cout<<"higgs [0] am: "<<get<1>(higgs_info[0])<<" dm: "<<get<0>(higgs_info[0])<<endl;
  //cout<<"higgs [1] am: "<<get<1>(higgs_info[1])<<" dm: "<<get<0>(higgs_info[1])<<endl;
  //cout<<"higgs [2] am: "<<get<1>(higgs_info[2])<<" dm: "<<get<0>(higgs_info[2])<<endl;
  //cout<<"event end"<<endl;

  return get<2>(higgs_info[0]);
}

const NamedFunc h1b1_pt("h1b1_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_pt())[bbjet_indices.at(0)];
});

const NamedFunc h1b2_pt("h1b2_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_pt())[bbjet_indices.at(1)];
});

const NamedFunc h2b1_pt("h2b1_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_pt())[bbjet_indices.at(2)];
});

const NamedFunc h2b2_pt("h2b2_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_pt())[bbjet_indices.at(3)];
});


const NamedFunc h1b1_eta("h1b1_eta",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_eta())[bbjet_indices.at(0)];
});

const NamedFunc h1b2_eta("h1b2_eta",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_eta())[bbjet_indices.at(1)];
});

const NamedFunc h2b1_eta("h2b1_eta",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_eta())[bbjet_indices.at(2)];
});

const NamedFunc h2b2_eta("h2b2_eta",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_eta())[bbjet_indices.at(3)];
});

const NamedFunc h1b1_jetid("h1b1_jetid",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_id())[bbjet_indices.at(0)];
});

const NamedFunc h1b2_jetid("h1b2_jetid",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_id())[bbjet_indices.at(1)];
});

const NamedFunc h2b1_jetid("h2b1_jetid",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_id())[bbjet_indices.at(2)];
});

const NamedFunc h2b2_jetid("h2b2_jetid",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  return (*b.jet_id())[bbjet_indices.at(3)];
});


const NamedFunc h1_dr("h1_dr",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  return DeltaR(h1b1, h1b2);
});

const NamedFunc h2_dr("h2_dr",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  return DeltaR(h2b1, h2b2);
});

const NamedFunc h1_mass("h1_mass",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  ROOT::Math::PtEtaPhiMVector higgs1 = h1b1 + h1b2;
  return higgs1.M();
});

const NamedFunc h2_mass("h2_mass",[](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  ROOT::Math::PtEtaPhiMVector higgs2 = h2b1 + h2b2;
  return higgs2.M();
});

const NamedFunc lead_signal_lepton_pt("lead_signal_lepton_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  // Search for signal electrons
  float lead_electron_pt = -1;
  for (unsigned iEl = 0; iEl < b.el_sig()->size(); ++iEl) {
    if (b.el_sig()->at(iEl)) {
      lead_electron_pt = b.el_pt()->at(iEl); 
      break;
    }
  }
  // Search for signal muons
  float lead_muon_pt = -1;
  for (unsigned iMu = 0; iMu < b.mu_sig()->size(); ++iMu) {
    if (b.mu_sig()->at(iMu)) {
      lead_muon_pt = b.mu_pt()->at(iMu); 
      break;
    }
  }
  // Warnings
  if (lead_electron_pt==-1 && lead_muon_pt==-1) {
    //cout<<"[Warning] Higfuncs::lead_signal_lepton_pt => There are no signal leptons. Returning -1. nlep: "<<b.nlep()<<" nel: "<<b.nel()<<" nmu: "<<b.nmu()<<" nvmu: "<<b.nvmu()<<endl;
    return -1;
  } else if (lead_electron_pt != -1 && lead_muon_pt != -1) {
    // Return max pt
    return max(lead_electron_pt, lead_muon_pt);
  } else if (lead_electron_pt != -1 && lead_muon_pt == -1) { // Electron case
    return lead_electron_pt;
  } else { // Muon case
    return lead_muon_pt;
  }
});

const NamedFunc lead_signal_muon_pt("lead_signal_muon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  // Search for signal muons
  float lead_muon_pt = -1;
  for (unsigned iMu = 0; iMu < b.mu_sig()->size(); ++iMu) {
    if (b.mu_sig()->at(iMu)) {
      lead_muon_pt = b.mu_pt()->at(iMu); 
      break;
    }
  }
  // Warnings
  if (lead_muon_pt==-1) {
    cout<<"[Warning] Higfuncs::lead_signal_muon_pt => There is no signal muons. Returning -1. nlep: "<<b.nlep()<<" nel: "<<b.nel()<<" nmu: "<<b.nmu()<<" nvmu: "<<b.nvmu()<<endl;
    return -1;
  } else { // Muon case
    return lead_muon_pt;
  }
});

const NamedFunc lead_signal_electron_pt("lead_signal_electron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  // Search for signal electrons
  float lead_electron_pt = -1;
  for (unsigned iEl = 0; iEl < b.el_sig()->size(); ++iEl) {
    if (b.el_sig()->at(iEl)) {
      lead_electron_pt = b.el_pt()->at(iEl); 
      break;
    }
  }
  // Warnings
  if (lead_electron_pt==-1) {
    cout<<"[Warning] Higfuncs::lead_signal_electron_pt => There is no electron leptons. Returning -1. nlep: "<<b.nlep()<<" nel: "<<b.nel()<<" nmu: "<<b.nmu()<<" nvmu: "<<b.nvmu()<<endl;
    return -1;
  } else {
    return lead_electron_pt;
  }
});

const NamedFunc pass_ecalnoisejet("pass_ecalnoisejet", [](const Baby &b) -> NamedFunc::ScalarType{
  //check top two highest pt jets in eta 2.4 to 5.0 region, if either has pt>250 and is closely aligned or anti-aligned with MET, then the event is vetoed
  if ((b.type()%100) <= (26*2)) return true; //skip 2016
  int counter = 0;
  bool r_pass_ecalnoisejet;
  bool goodjet[2] = {true, true};
  double dphi = 0.;
  for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
    if (counter >= 2) break;
    if (b.jet_pt()->at(jet_idx)>30 && fabs(b.jet_eta()->at(jet_idx))>2.4 && fabs(b.jet_eta()->at(jet_idx))<5.0) {
      dphi = fabs(TVector2::Phi_mpi_pi(b.jet_phi()->at(jet_idx)-b.met_phi()));
      if (b.jet_pt()->at(jet_idx)>250 && (dphi > 2.6 || dphi < 0.1)) goodjet[counter] = false;
      ++counter;
    }
  }
  r_pass_ecalnoisejet = goodjet[0] && goodjet[1];
  return r_pass_ecalnoisejet;
});

const NamedFunc pass_hemveto("pass_hemveto", [](const Baby &b) -> NamedFunc::ScalarType{
    //only apply for 2018 era C+D and MC
    if (abs(b.SampleType())!=2018) return true; // accept 2016 and 2017 data and mc
    if (b.SampleType()==-2018 && b.run() < 319077) return true; // accept part of 2018 data
    if (b.SampleType()==2018 && (b.event()%1961) >= 1296) return true; //accept part of 2018 mc
    bool pass_hem = true;
    for (unsigned int el_idx = 0; el_idx < b.el_pt()->size(); el_idx++) {
      if (b.el_miniso()->at(el_idx) < 0.1 && -3.0 < b.el_eta()->at(el_idx) && b.el_eta()->at(el_idx) < -1.4 && -1.57 < b.el_phi()->at(el_idx) && b.el_phi()->at(el_idx) < -0.87) {
        pass_hem = false;
      }
    }
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_pt()->at(jet_idx) > 30. && -3.2 < b.jet_eta()->at(jet_idx) && b.jet_eta()->at(jet_idx) < -1.2 && -1.77 < b.jet_phi()->at(jet_idx) && b.jet_phi()->at(jet_idx) < -0.67) {
        double dphi = fabs(TVector2::Phi_mpi_pi(b.jet_phi()->at(jet_idx)-b.met_phi()));
        if (dphi < 0.5) {
          pass_hem = false;
        }
      }
    }
    return pass_hem;
});

const NamedFunc hem_weight("hem_weight", [](const Baby & b) -> NamedFunc::ScalarType{
  if ( b.type()/1000 == 0) return 1.0; //data
  if ( abs(b.SampleType()) != 2018) return 1.0;
  if ( pass_hemveto.GetScalar(b) ) return 1.0;
  return 0.339113; //fraction of lumi before hem failure
});

const NamedFunc pass_filters("pass_filters", [](const Baby &b) -> NamedFunc::ScalarType{
  if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
  if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
  if ((b.type()/1000 != 106)  && !b.pass_cschalo_tight()) return false; //not for fastsim
  if (!b.pass_low_neutral_jet()) return false;
  if (!b.pass_htratio_dphi_tight()) return false;
  if ((b.type()/1000 == 106)  && !b.pass_jets()) return false; //back to only for fastsim
  //if (!b.pass_jets()) return false; //was modified
  if ((abs(b.SampleType())==2017 || abs(b.SampleType())==2018) && !Higfuncs::pass_ecalnoisejet.GetScalar(b)) return false; 
  if (!Higfuncs::pass_hemveto.GetScalar(b)) return false;
  return true;
});

const NamedFunc final_pass_filters = pass_filters&& "met/mht<2 && met/met_calo<2&&weight<1.5&&pass_jets";
const NamedFunc final_ttbar_pass_filters = pass_filters&& "met/met_calo<5&&weight<1.5&&pass_jets";
const NamedFunc final_zll_pass_filters = pass_filters&& "met/met_calo<5&&weight<1.5&&pass_jets";
const NamedFunc final_qcd_pass_filters = pass_filters&& "met/mht<2 && met/met_calo<2&&pass_jets";

const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double weight = 1;
  if (b.SampleType()==2016){
    return weight*35.92; // prev 35.9
  } else if (b.SampleType()==2017){
    return weight*41.53; // prev 41.5
  } else {
    return weight*59.74; // prev 59.6
  }
});

const NamedFunc final_weight = "weight"*eff_higtrig_run2*w_years*Functions::w_pileup;
const NamedFunc final_weight_notrgeff = "weight"*w_years*Functions::w_pileup;


const NamedFunc w_pileup_nosignal("w_pileup_nosignal",[](const Baby &b) -> NamedFunc::ScalarType{
  if ((b.type()/1000) == 106)  return 1.0;
  return Functions::w_pileup.GetScalar(b);
});

const NamedFunc jet_trigger = "HLT_PFJet500";

const NamedFunc met_trigger("met_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
  //can't used string-based named func because name being too long causes histograms to crash
  bool r_met_trigger = b.HLT_PFMET90_PFMHT90_IDTight()||b.HLT_PFMETNoMu90_PFMHTNoMu90_IDTight()||b.HLT_PFMET100_PFMHT100_IDTight()||b.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight()||b.HLT_PFMET110_PFMHT110_IDTight()||b.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight()||b.HLT_PFMET120_PFMHT120_IDTight()||b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight()||b.HLT_PFMET130_PFMHT130_IDTight()||b.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight()||b.HLT_PFMET140_PFMHT140_IDTight()||b.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight()||b.HLT_PFMET100_PFMHT100_IDTight_PFHT60()||b.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60()||b.HLT_PFMET110_PFMHT110_IDTight_PFHT60()||b.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60()||b.HLT_PFMET120_PFMHT120_IDTight_PFHT60()||b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60()||b.HLT_PFMET130_PFMHT130_IDTight_PFHT60()||b.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60()||b.HLT_PFMET140_PFMHT140_IDTight_PFHT60()||b.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60()||b.HLT_PFMET120_PFMHT120_IDTight_HFCleaned()||b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned()||b.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned();
  return r_met_trigger;
});

const NamedFunc el_trigger("el_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
  //can't used string-based named func because name being too long causes histograms to crash
  bool r_el_trigger = b.HLT_Ele25_WPTight_Gsf()||b.HLT_Ele27_WPTight_Gsf()||b.HLT_Ele28_WPTight_Gsf()||b.HLT_Ele32_WPTight_Gsf_L1DoubleEG()||b.HLT_Ele32_WPTight_Gsf()||b.HLT_Ele35_WPTight_Gsf()||b.HLT_Ele45_WPLoose_Gsf()||b.HLT_Ele105_CaloIdVT_GsfTrkIdT()||b.HLT_Ele115_CaloIdVT_GsfTrkIdT()||b.HLT_Ele135_CaloIdVT_GsfTrkIdT()||b.HLT_Ele145_CaloIdVT_GsfTrkIdT()||b.HLT_Ele25_eta2p1_WPTight_Gsf()||b.HLT_Ele27_eta2p1_WPTight_Gsf()||b.HLT_Ele27_eta2p1_WPLoose_Gsf()||b.HLT_Ele20_WPLoose_Gsf()||b.HLT_Ele20_eta2p1_WPLoose_Gsf()||b.HLT_Ele25_eta2p1_WPLoose_Gsf()||b.HLT_Ele15_IsoVVVL_PFHT350()||b.HLT_Ele15_IsoVVVL_PFHT400()||b.HLT_Ele15_IsoVVVL_PFHT450()||b.HLT_Ele15_IsoVVVL_PFHT600()||b.HLT_Ele50_IsoVVVL_PFHT450();
  return r_el_trigger;
});

const NamedFunc mu_trigger("mu_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
  //can't used string-based named func because name being too long causes histograms to crash
  bool r_mu_trigger = b.HLT_IsoMu20()||b.HLT_IsoMu22()||b.HLT_IsoMu24()||b.HLT_IsoMu27()||b.HLT_IsoTkMu20()||b.HLT_IsoTkMu22()||b.HLT_IsoTkMu24()||b.HLT_Mu50()||b.HLT_Mu55()||b.HLT_TkMu50()||b.HLT_IsoMu22_eta2p1()||b.HLT_IsoMu24_eta2p1()||b.HLT_Mu45_eta2p1()||b.HLT_Mu15_IsoVVVL_PFHT350()||b.HLT_Mu15_IsoVVVL_PFHT400()||b.HLT_Mu15_IsoVVVL_PFHT450()||b.HLT_Mu15_IsoVVVL_PFHT600()||b.HLT_Mu50_IsoVVVL_PFHT400()||b.HLT_Mu50_IsoVVVL_PFHT450();
  return r_mu_trigger;
});

const NamedFunc ht_trigger("ht_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
  //can't used string-based named func because name being too long causes histograms to crash
  bool r_ht_trigger = b.HLT_PFHT125()||b.HLT_PFHT200()||b.HLT_PFHT300()||b.HLT_PFHT400()||b.HLT_PFHT475()||b.HLT_PFHT600()||b.HLT_PFHT650()||b.HLT_PFHT800()||b.HLT_PFHT900()||b.HLT_PFHT180()||b.HLT_PFHT370()||b.HLT_PFHT430()||b.HLT_PFHT510()||b.HLT_PFHT590()||b.HLT_PFHT680()||b.HLT_PFHT780()||b.HLT_PFHT890()||b.HLT_PFHT1050()||b.HLT_PFHT250()||b.HLT_PFHT350();
  return r_ht_trigger;
});

const NamedFunc jetht_trigger = "HLT_PFJet500" || ht_trigger;

const NamedFunc met_trigger_v0 = "(HLT_PFMET110_PFMHT110_IDTight || HLT_PFMETNoMu110_PFMHTNoMu110_IDTight || HLT_PFMET120_PFMHT120_IDTight || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_PFMET120_PFMHT120_IDTight_PFHT60 || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)";

const NamedFunc el_trigger_v0 = "(HLT_Ele27_WPTight_Gsf || HLT_Ele35_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT)";

const NamedFunc mu_trigger_v0 = "(HLT_IsoMu24 || HLT_IsoMu27 || HLT_Mu50)";

const NamedFunc jetid_njet("jetid_njet", [](const Baby &b) -> NamedFunc::ScalarType{
  unsigned int njet = 0;
  if (b.pass_jets()) return b.njet();
  for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
    if (!b.jet_isgood()->at(jet_idx)) continue;
    if (b.jet_id()->at(jet_idx) == false) continue;
    njet += 1;
  }
  return njet;
});

const NamedFunc jetid_nb("jetid_nb", [](const Baby &b) -> NamedFunc::ScalarType{
  unsigned int nbt = 0, nbm = 0, nbl = 0;
  if (b.pass_jets()) {
    if (b.nbt() < 2) return b.nbt();
    if (b.nbm() < 3) return 2;
    if (b.nbl()==3) return 3;
    return 4;
  }
  for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
    if (!b.jet_isgood()->at(jet_idx)) continue;
    if (b.jet_id()->at(jet_idx) == false) continue;
    if (b.jet_deepcsv()->at(jet_idx) > 0.8953 && abs(b.SampleType()) == 2016) nbt += 1;
    if (b.jet_deepcsv()->at(jet_idx) > 0.8001 && abs(b.SampleType()) == 2017) nbt += 1;
    if (b.jet_deepcsv()->at(jet_idx) > 0.7527 && abs(b.SampleType()) == 2018) nbt += 1;
    if (b.jet_deepcsv()->at(jet_idx) > 0.6321 && abs(b.SampleType()) == 2016) nbm += 1;
    if (b.jet_deepcsv()->at(jet_idx) > 0.4941 && abs(b.SampleType()) == 2017) nbm += 1;
    if (b.jet_deepcsv()->at(jet_idx) > 0.4184 && abs(b.SampleType()) == 2018) nbm += 1;
    if (b.jet_deepcsv()->at(jet_idx) > 0.2217 && abs(b.SampleType()) == 2016) nbl += 1;
    if (b.jet_deepcsv()->at(jet_idx) > 0.1522 && abs(b.SampleType()) == 2017) nbl += 1;
    if (b.jet_deepcsv()->at(jet_idx) > 0.0494 && abs(b.SampleType()) == 2018) nbl += 1;
  }
  if (nbt < 2) return nbt;
  if (nbm < 3) return 2;
  if (nbl==3) return 3;
  return 4;
});

const NamedFunc jetid_low_dphi_met("jetid_low_dphi_met", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.pass_jets()) return b.low_dphi_met();
  std::vector<std::pair<float, float>> jet_pt_and_dphi;
  unsigned int njet = 0;
  for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
    if (!b.jet_isgood()->at(jet_idx)) continue;
    if (b.jet_id()->at(jet_idx) == false) continue;
    njet += 1;
    jet_pt_and_dphi.push_back(std::make_pair(b.jet_pt()->at(jet_idx), b.jet_met_dphi()->at(jet_idx)));
  }
  if (njet > 0) {
    if (jet_pt_and_dphi[0].second < 0.5) return true;
  }
  if (njet > 1) {
    if (jet_pt_and_dphi[1].second < 0.5) return true;
  }
  if (njet > 2) {
    if (jet_pt_and_dphi[2].second < 0.3) return true;
  }
  if (njet > 3) {
    if (jet_pt_and_dphi[3].second < 0.3) return true;
  }
  return false;
});

const NamedFunc jetid_hig_cand_dm("jetid_hig_cand_dm",[](const Baby &b) -> NamedFunc::ScalarType{
  //if no bad jets or if bad jet not in higgs candidate, just return normal value
  if (b.pass_jets()) {
    if (b.njet() < 4) return 999.;
    return b.hig_cand_dm()->at(0);
  }
  std::vector<bool> jet_isgood;
  bool recalculate_higgs = false;
  unsigned int njet = 0;
  for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
    jet_isgood.push_back(b.jet_isgood()->at(jet_idx) && b.jet_id()->at(jet_idx));
    if (b.jet_id()->at(jet_idx) == false) {
      if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) recalculate_higgs = true;
      continue;
    }
    if (!b.jet_isgood()->at(jet_idx)) continue;
    njet += 1;
  }
  if (njet < 4) return 999.;
  if (!recalculate_higgs) return b.hig_cand_dm()->at(0);
  //recalculate higgs
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), jet_isgood);
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  return TMath::Abs((h1b1+h1b2).M()-(h2b1+h2b2).M());
});

const NamedFunc jetid_hig_cand_am("jetid_hig_cand_am",[](const Baby &b) -> NamedFunc::ScalarType{
  //if no bad jets or if bad jet not in higgs candidate, just return normal value
  if (b.pass_jets()) {
    if (b.njet() < 4) return 999.;
    return b.hig_cand_am()->at(0);
  }
  std::vector<bool> jet_isgood;
  bool recalculate_higgs = false;
  unsigned int njet = 0;
  for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
    jet_isgood.push_back(b.jet_isgood()->at(jet_idx) && b.jet_id()->at(jet_idx));
    if (b.jet_id()->at(jet_idx) == false) {
      if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) recalculate_higgs = true;
      continue;
    }
    if (!b.jet_isgood()->at(jet_idx)) continue;
    njet += 1;
  }
  if (njet < 4) return 999.;
  if (!recalculate_higgs) return b.hig_cand_am()->at(0);
  //recalculate higgs
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), jet_isgood);
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  return ((h1b1+h1b2).M()+(h2b1+h2b2).M())/2.0;
});

const NamedFunc jetid_hig_cand_drmax("jetid_hig_cand_drmax",[](const Baby &b) -> NamedFunc::ScalarType{
  //if no bad jets or if bad jet not in higgs candidate, just return normal value
  if (b.pass_jets()) {
    if (b.njet() < 4) return 999.;
    return b.hig_cand_drmax()->at(0);
  }
  std::vector<bool> jet_isgood;
  bool recalculate_higgs = false;
  unsigned int njet = 0;
  for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
    jet_isgood.push_back(b.jet_isgood()->at(jet_idx) && b.jet_id()->at(jet_idx));
    if (b.jet_id()->at(jet_idx) == false) {
      if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) recalculate_higgs = true;
      continue;
    }
    if (!b.jet_isgood()->at(jet_idx)) continue;
    njet += 1;
  }
  if (njet < 4) return 999.;
  if (!recalculate_higgs) return b.hig_cand_drmax()->at(0);
  //recalculate higgs
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), jet_isgood);
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  float dr1 = DeltaR(h1b1, h1b2);
  float dr2 = DeltaR(h2b1, h2b2);
  if (dr1>dr2) return dr1;
  return dr2;
});

}
