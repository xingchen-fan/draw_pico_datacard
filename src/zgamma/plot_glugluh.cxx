//This script generates some tables and plots for gg->H->Zg studies

#include "core/test.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "Math/Boost.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "TColor.h"
#include "TError.h"
#include "TMath.h"

#include "core/baby.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"
#include "zgamma/apply_zg_trigeffs.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

//#include "../dnns/dnn_zg_pos.hpp"
//#include "../dnns/dnn_zg_neg.hpp"
//#include "../dnns/dnn_dy_pos.hpp"
//#include "../dnns/dnn_dy_neg.hpp"
//#include "../dnns/dnn_sig.hpp"

//costheta using unified coordinate system - see below functions
template<class coords1, class coords2>
double costheta_threevec(ROOT::Math::LorentzVector<coords1> lv1, 
                         ROOT::Math::LorentzVector<coords2> lv2) {
  ROOT::Math::XYZVector v1 = lv1.Vect();
  ROOT::Math::XYZVector v2 = lv2.Vect();
  return v1.Dot(v2)/v1.R()/v2.R();
}

//phi using unified coordinate system - see below functions
template<class coords1, class coords2, class coords3>
double phi_threevec(ROOT::Math::LorentzVector<coords1> lv, 
                    ROOT::Math::LorentzVector<coords2> lv_zref, 
                    ROOT::Math::LorentzVector<coords3> lv_xref) {
  ROOT::Math::XYZVector v_target = lv.Vect();
  ROOT::Math::XYZVector z_unit = lv_zref.Vect()/lv_zref.Vect().R();
  ROOT::Math::XYZVector y_unit = lv_zref.Vect().Cross(lv_xref.Vect());
  y_unit = y_unit/y_unit.R();
  ROOT::Math::XYZVector x_unit = y_unit.Cross(z_unit);
  v_target = v_target - v_target.Dot(z_unit)*z_unit;
  float abs_x_angle = TMath::ACos(v_target.Dot(x_unit)/v_target.R());
  float abs_y_angle = TMath::ACos(v_target.Dot(y_unit)/v_target.R());
  if (abs_y_angle < TMath::Pi()/2.0) return abs_x_angle;
  return -1.0*abs_x_angle;
}

//returns 4-vector of +z gluon in lab frame assuming the usual convention
template<class coords1, class coords2, class coords3>
ROOT::Math::PxPyPzEVector poszgluon(ROOT::Math::LorentzVector<coords1> const & p_neglep, 
                                    ROOT::Math::LorentzVector<coords2> const & p_poslep, 
                                    ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //set up boost to "zgtransverse" frame
  ROOT::Math::PtEtaPhiMVector p_zph = (p_poslep+p_neglep+p_photon);
  ROOT::Math::XYZVector boost_vec_lab_to_zgtrans = p_zph.BoostToCM();
  boost_vec_lab_to_zgtrans.SetZ(0.0);
  ROOT::Math::Boost boost_lab_to_zgtrans(boost_vec_lab_to_zgtrans);
  ROOT::Math::Boost boost_zgtrans_to_lab(boost_lab_to_zgtrans.Inverse());

  //calculate +z gluon in "zgtransverse" frame and boost back
  ROOT::Math::PtEtaPhiMVector p_zph_zgtrans(boost_lab_to_zgtrans(p_zph));
  double partz_zgtrans = p_zph_zgtrans.Pz()+p_zph_zgtrans.E();
  ROOT::Math::PxPyPzEVector p_g1(0.0,0.0,partz_zgtrans,partz_zgtrans); 
  return boost_zgtrans_to_lab(p_g1);
}

//returns cosTheta of Z in lly CM frame with +z gluon as z reference and lab +x as X reference
template<class coords1, class coords2, class coords3>
float coscaptheta(ROOT::Math::LorentzVector<coords1> const & p_neglep, 
                  ROOT::Math::LorentzVector<coords2> const & p_poslep, 
                  ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //initialize necessary 4-vectors and boosts
  ROOT::Math::PtEtaPhiMVector p_z = (p_poslep+p_neglep);
  ROOT::Math::PtEtaPhiMVector p_zph = (p_poslep+p_neglep+p_photon);
  ROOT::Math::Boost boost_lab_to_zgcm(p_zph.BoostToCM());

  //boost to CM frame and calculate
  ROOT::Math::PxPyPzEVector p_g1_cm = boost_lab_to_zgcm(poszgluon(p_neglep,p_poslep,p_photon));
  ROOT::Math::PtEtaPhiMVector p_z_cm = boost_lab_to_zgcm(p_z);
  return costheta_threevec(p_z_cm, p_g1_cm);
}

//returns Phi of Z in lly CM frame with +z gluon as z reference and lab +x as X reference
template<class coords1, class coords2, class coords3>
float capphi(ROOT::Math::LorentzVector<coords1> const & p_neglep, 
             ROOT::Math::LorentzVector<coords2> const & p_poslep, 
             ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //initialize necessary 4-vectors and boosts
  ROOT::Math::PxPyPzEVector p_xlab(1.0,0.0,0.0,0.0); 
  ROOT::Math::PtEtaPhiMVector p_z = (p_poslep+p_neglep);
  ROOT::Math::PtEtaPhiMVector p_zph = (p_poslep+p_neglep+p_photon);
  ROOT::Math::Boost boost_lab_to_zgcm(p_zph.BoostToCM());

  //boost to CM frame and calculate
  ROOT::Math::PxPyPzEVector p_g1_cm = boost_lab_to_zgcm(poszgluon(p_neglep,p_poslep,p_photon));
  ROOT::Math::PxPyPzEVector p_xlab_cm = boost_lab_to_zgcm(p_xlab);
  ROOT::Math::PtEtaPhiMVector p_z_cm = boost_lab_to_zgcm(p_z);
  return phi_threevec(p_z_cm, p_g1_cm, p_xlab_cm);
}

//returns costheta of positive lepton in ll CM frame with photon as z reference and +z gluon as x reference
template<class coords1, class coords2, class coords3>
float coslowertheta(ROOT::Math::LorentzVector<coords1> const & p_neglep, 
                    ROOT::Math::LorentzVector<coords2> const & p_poslep, 
                    ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //initialize necessary 4-vectors and boosts
  ROOT::Math::PtEtaPhiMVector p_z = (p_poslep+p_neglep);
  ROOT::Math::Boost boost_lab_to_zcm(p_z.BoostToCM());

  //boost to CM frame and calculate
  ROOT::Math::PtEtaPhiMVector p_poslep_z = boost_lab_to_zcm(p_poslep);
  ROOT::Math::PtEtaPhiMVector p_photon_z = boost_lab_to_zcm(p_photon);
  return costheta_threevec(p_poslep_z, p_photon_z);
}

//returns phi of positive lepton in ll CM frame with photon as z reference and +z gluon as x reference
template<class coords1, class coords2, class coords3>
float lowerphi(ROOT::Math::LorentzVector<coords1> const & p_neglep, 
               ROOT::Math::LorentzVector<coords2> const & p_poslep, 
               ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //initialize necessary 4-vectors and boosts
  ROOT::Math::PtEtaPhiMVector p_z = (p_poslep+p_neglep);
  ROOT::Math::Boost boost_lab_to_zcm(p_z.BoostToCM());

  //boost to CM frame and calculate
  ROOT::Math::PxPyPzEVector p_g1_z = boost_lab_to_zcm(poszgluon(p_neglep,p_poslep,p_photon));
  ROOT::Math::PtEtaPhiMVector p_poslep_z = boost_lab_to_zcm(p_poslep);
  ROOT::Math::PtEtaPhiMVector p_photon_z = boost_lab_to_zcm(p_photon);
  return phi_threevec(p_poslep_z, p_photon_z, p_g1_z);
}

ROOT::Math::PtEtaPhiMVector GetLposLV(const Baby& b) {
  if (b.ll_lepid()->at(0)==11) {
    if (b.el_charge()->at(b.ll_i1()->at(0))==1) {
      return ROOT::Math::PtEtaPhiMVector(b.el_pt()->at(b.ll_i1()->at(0)),
          b.el_eta()->at(b.ll_i1()->at(0)),b.el_phi()->at(b.ll_i1()->at(0)),0);
    }
    return ROOT::Math::PtEtaPhiMVector(b.el_pt()->at(b.ll_i2()->at(0)),
        b.el_eta()->at(b.ll_i2()->at(0)),b.el_phi()->at(b.ll_i2()->at(0)),0);
  }
  if (b.mu_charge()->at(b.ll_i1()->at(0))==1) {
    return ROOT::Math::PtEtaPhiMVector(b.mu_pt()->at(b.ll_i1()->at(0)),
        b.mu_eta()->at(b.ll_i1()->at(0)),b.mu_phi()->at(b.ll_i1()->at(0)),0);
  }
  return ROOT::Math::PtEtaPhiMVector(b.mu_pt()->at(b.ll_i2()->at(0)),
      b.mu_eta()->at(b.ll_i2()->at(0)),b.mu_phi()->at(b.ll_i2()->at(0)),0);
}

ROOT::Math::PtEtaPhiMVector GetLnegLV(const Baby& b) {
  if (b.ll_lepid()->at(0)==11) {
    if (b.el_charge()->at(b.ll_i1()->at(0))==-1) {
      return ROOT::Math::PtEtaPhiMVector(b.el_pt()->at(b.ll_i1()->at(0)),
          b.el_eta()->at(b.ll_i1()->at(0)),b.el_phi()->at(b.ll_i1()->at(0)),0);
    }
    return ROOT::Math::PtEtaPhiMVector(b.el_pt()->at(b.ll_i2()->at(0)),
        b.el_eta()->at(b.ll_i2()->at(0)),b.el_phi()->at(b.ll_i2()->at(0)),0);
  }
  if (b.mu_charge()->at(b.ll_i1()->at(0))==-1) {
    return ROOT::Math::PtEtaPhiMVector(b.mu_pt()->at(b.ll_i1()->at(0)),
        b.mu_eta()->at(b.ll_i1()->at(0)),b.mu_phi()->at(b.ll_i1()->at(0)),0);
  }
  return ROOT::Math::PtEtaPhiMVector(b.mu_pt()->at(b.ll_i2()->at(0)),
      b.mu_eta()->at(b.ll_i2()->at(0)),b.mu_phi()->at(b.ll_i2()->at(0)),0);
}

ROOT::Math::PtEtaPhiMVector GetGammaLV(const Baby& b) {
  return ROOT::Math::PtEtaPhiMVector(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0),0);
}

const NamedFunc lower_phi("lower_phi",[](const Baby &b) -> NamedFunc::ScalarType{
  return lowerphi(GetLnegLV(b), GetLposLV(b), GetGammaLV(b));
});

const NamedFunc cap_phi("cap_phi",[](const Baby &b) -> NamedFunc::ScalarType{
  return capphi(GetLnegLV(b), GetLposLV(b), GetGammaLV(b));
});

const NamedFunc cos_lower_theta("cos_lower_theta",[](const Baby &b) -> NamedFunc::ScalarType{
  return coslowertheta(GetLnegLV(b), GetLposLV(b), GetGammaLV(b));
});

const NamedFunc cos_cap_theta("cos_cap_theta",[](const Baby &b) -> NamedFunc::ScalarType{
  return coscaptheta(GetLnegLV(b), GetLposLV(b), GetGammaLV(b));
});

//const dnn_zg_pos zg_dnn_pos;
//const dnn_zg_neg zg_dnn_neg;
//const float zg_neg_frac = 0.14228788; //120<mllg<130
//const dnn_dy_pos dy_dnn_pos;
//const dnn_dy_neg dy_dnn_neg;
//const float dy_neg_frac = 0.21969835; //120<mllg<130
//const dnn_sig sig_dnn;
//
//float zg_pdf_pos(std::vector<float> vars) {
//  float dnn_score = zg_dnn_pos.evaluate(vars);
//  if (dnn_score==1.0) return 0.9999/(1.0-0.9999);
//  return dnn_score/(1.0-dnn_score)/485330.0;
//}
//
//float zg_pdf_neg(std::vector<float> vars) {
//  float dnn_score = zg_dnn_neg.evaluate(vars);
//  if (dnn_score==1.0) return 0.9999/(1.0-0.9999);
//  return dnn_score/(1.0-dnn_score)/510961.0;
//}
//
//float zg_pdf_net(std::vector<float> vars) {
//  float pos = zg_pdf_pos(vars);
//  float neg = zg_pdf_neg(vars);
//  float result = (-1.0*zg_neg_frac*neg+(1.0-zg_neg_frac)*pos)/(1.0-2.0*zg_neg_frac);
//  if (result < 0) return 0;
//  return result;
//}
//
//float dy_pdf_pos(std::vector<float> vars) {
//  float dnn_score = dy_dnn_pos.evaluate(vars);
//  if (dnn_score==1.0) return 0.9999/(1.0-0.9999);
//  return dnn_score/(1.0-dnn_score)/416964.0;
//}
//
//float dy_pdf_neg(std::vector<float> vars) {
//  float dnn_score = dy_dnn_neg.evaluate(vars);
//  if (dnn_score==1.0) return 0.9999/(1.0-0.9999);
//  return dnn_score/(1.0-dnn_score)/373446.0;
//}
//
//float dy_pdf_net(std::vector<float> vars) {
//  float pos = dy_pdf_pos(vars);
//  float neg = dy_pdf_neg(vars);
//  float result = (-1.0*dy_neg_frac*neg+(1.0-dy_neg_frac)*pos)/(1.0-2.0*dy_neg_frac);
//  if (result < 0) return 0;
//  return result;
//}
//
//float sig_pdf_net(std::vector<float> vars) {
//  float dnn_score = sig_dnn.evaluate(vars);
//  if (dnn_score==1.0) return 0.9999/(1.0-0.9999);
//  return dnn_score/(1.0-dnn_score)/399210.0; 
//}
//
//const NamedFunc dnn_likelihood_ratio("dnn_likelihood_ratio",[](const Baby &b) -> NamedFunc::ScalarType{
//    std::vector<float> vars = {b.llphoton_pt()->at(0), b.llphoton_eta()->at(0),
//        static_cast<float>(cos_cap_theta.GetScalar(b)), static_cast<float>(cap_phi.GetScalar(b)), 
//        static_cast<float>(cos_lower_theta.GetScalar(b)), static_cast<float>(lower_phi.GetScalar(b)),
//        b.photon_idmva()->at(0)};
//    float prob_bk = 0.68525*zg_pdf_net(vars)+0.31476*dy_pdf_net(vars);
//    float prob_sg = sig_pdf_net(vars);
//    return 1.0/(1.0+prob_bk/prob_sg);
//});

using namespace std;
using namespace PlotOptTypes;
using namespace ZgFunctions;
using namespace ZgUtilities;

int main() {

  //--------------------------------------------------------------------------
  //                          define NamedFuncs
  //--------------------------------------------------------------------------

  NamedFunc w_sigx20("w_sigx20",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() >= 200000 && b.type() <= 205000)
      return 20.0;
    return 1.0;
  });

  NamedFunc el1_pt = "el_pt[ll_i1[0]]";
  NamedFunc el2_pt = "el_pt[ll_i2[0]]";
  NamedFunc mu1_pt = "mu_pt[ll_i1[0]]";
  NamedFunc mu2_pt = "mu_pt[ll_i2[0]]";
  NamedFunc el1_eta = "el_eta[ll_i1[0]]";
  NamedFunc el2_eta = "el_eta[ll_i2[0]]";
  NamedFunc mu1_eta = "mu_eta[ll_i1[0]]";
  NamedFunc mu2_eta = "mu_eta[ll_i2[0]]";
  NamedFunc photon_pt = "photon_pt[0]";
  NamedFunc llphoton_cosTheta = "llphoton_cosTheta[0]";
  llphoton_cosTheta.Name("llphoton_cosCapTheta0");

  NamedFunc lep1_eta("lep1_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.ll_lepid()->at(0)==11) return b.el_eta()->at(b.ll_i1()->at(0));
    return b.mu_eta()->at(b.ll_i1()->at(0));
  });

  NamedFunc lep2_eta("lep2_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.ll_lepid()->at(0)==11) return b.el_eta()->at(b.ll_i2()->at(0));
    return b.mu_eta()->at(b.ll_i2()->at(0));
  });

  NamedFunc zg_el_cuts = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)";
  zg_el_cuts.Name("zg_el_cuts");
  NamedFunc zg_mu_cuts = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)";
  zg_mu_cuts.Name("zg_mu_cuts");
  NamedFunc higgs_window = "llphoton_m[0]>120&&llphoton_m[0]<130";

  //fix weight for Zgamma sample
  const NamedFunc w_ewkzgamma("w_ewkzgamma",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.FirstFileName().find("ZGamma2JToGamma2L2J_EWK") != std::string::npos) {
      float xs = 0.1145*1000.0; //fb
      float year_nevts = 230000.0; //2016
      if (abs(b.SampleType())==2016) year_nevts = 500000.0;
      if (abs(b.SampleType())==2017) year_nevts = 498000.0;
      return xs/year_nevts;
    }
    return b.w_lumi();
  });

  const NamedFunc tighter_baseline = NamedFunc(zg_baseline&&"photon_id80[0]&&ll_m[0]>80&&ll_m[0]<100").Name("new_baseline");

  //--------------------------------------------------------------------------
  //                          setup
  //--------------------------------------------------------------------------
  
  gErrorIgnoreLevel = 6000;
  Palette mc_colors("txt/colors_zgamma.txt","default");
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;

  bool use_local = true;
  string prod_folder("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v2/");
  std::set<int> years = {2016,2017,2018};
  std::string lumi_string = "138";
  if (use_local) {
    prod_folder = "/data"+prod_folder;
  }

  bool multiply_sig = false;
  std::string sig_name = "gg#rightarrow H#rightarrow Z#gamma";
  NamedFunc weight = w_ewkzgamma*w_years;
  if (multiply_sig) {
    sig_name = "gg#rightarrow H#rightarrow Z#gamma (x20)";
    weight = w_ewkzgamma*w_sigx20*w_years;
  }

  std::cout << "Loading samples:\n";
  auto proc_smzg      = Process::MakeShared<Baby_pico>("Z+#gamma", back, mc_colors("zgtollg"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*ZGToLLG*"}), 
                           HLT_pass_dilepton);
  auto proc_dy        = Process::MakeShared<Baby_pico>("Z+Fake Photon", back, mc_colors("dyjets"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*DYJets*"}), 
                           "stitch_dy"&&HLT_pass_dilepton);
  auto proc_tt        = Process::MakeShared<Baby_pico>("tt", back, mc_colors("tt"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*TTTo2L2Nu*"}), 
                           HLT_pass_dilepton);
  auto proc_tt1l      = Process::MakeShared<Baby_pico>("Fake lepton", back, mc_colors("fakelep"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*TTG*"}), 
                           "(ntruel+ntrumu+ntrutaul)==1"&&HLT_pass_dilepton);
  //auto proc_hzg       = Process::MakeShared<Baby_pico>(sig_name, sig, kRed,
  //                         attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
  //                         {"*HToZG*M-125*","*HToZG_M125*"}), HLT_pass_dilepton);
  auto proc_hzg       = Process::MakeShared<Baby_pico>(sig_name, sig, kRed,
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                           {"*GluGluHToZG*M-125*"}), HLT_pass_dilepton);
  auto proc_hzg_bak   = Process::MakeShared<Baby_pico>(sig_name, back, kRed,
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                           {"*GluGluHToZG*M-125*"}), HLT_pass_dilepton);

  vector<shared_ptr<Process>> procs = {proc_dy, proc_tt, proc_smzg, proc_tt1l, proc_hzg};
  vector<shared_ptr<Process>> procs_bak = {proc_dy, proc_tt, proc_smzg, proc_tt1l};
  vector<shared_ptr<Process>> procs_sig = {proc_hzg_bak};

  std::vector<std::shared_ptr<Process>> procs_loaded = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","MinimalMoreTrigs");

  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.Title(TitleType::info)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::none)
          .Bottom(BottomType::off)
          .UseCMYK(false)
          .LegendColumns(1)
          .ShowBackgroundError(false)
          .FileExtensions({"pdf"});
  PlotOpt lin_sorb = lin_lumi;
  lin_sorb.Overflow(OverflowType::both).Bottom(BottomType::sorb).FileExtensions({"pdf","root"});
  PlotOpt lin_sorb_upper = lin_lumi;
  lin_sorb_upper.Overflow(OverflowType::both).Bottom(BottomType::sorb_cut_upper).FileExtensions({"pdf","root"});
  PlotOpt lin_shapes = lin_lumi;
  lin_shapes.Stack(StackType::shapes);
  PlotOpt log_lumi2d("txt/plot_styles.txt", "Scatter");
  log_lumi2d.Title(TitleType::info)
                 .YAxis(YAxisType::log)
                 .Overflow(OverflowType::overflow)
                 .LogMinimum(0.001);
  std::vector<PlotOpt> ops = {lin_lumi, lin_shapes};
  std::vector<PlotOpt> ops_sorb = {lin_shapes,lin_sorb};
  std::vector<PlotOpt> ops_sorb_upper = {lin_shapes,lin_sorb_upper};
  std::vector<PlotOpt> ops2d = {log_lumi2d};
  std::vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","LinLumi")
                                     .Overflow(OverflowType::none).FileExtensions({"pdf","root"}),
                                     PlotOpt("txt/plot_styles.txt","Shapes")
                                     .Overflow(OverflowType::none)}; 

  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------

  //kinematic MVA
  MVAWrapper kin_bdt_reader("kinematic_bdt");
  kin_bdt_reader.SetVariable("pt_llg","llphoton_pt[0]");
  kin_bdt_reader.SetVariable("eta_llg","llphoton_eta[0]");
  kin_bdt_reader.SetVariable("cosTheta",cos_cap_theta);
  kin_bdt_reader.SetVariable("Phi",cap_phi);
  kin_bdt_reader.SetVariable("costheta",cos_lower_theta);
  kin_bdt_reader.SetVariable("phi",lower_phi);
  kin_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
  kin_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_dnnstudies_density_BDT.weights.xml");
  NamedFunc bdt_score = kin_bdt_reader.GetDiscriminant();

  //kinematic MVA
  MVAWrapper std_bdt_reader("standard_bdt");
  std_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
  std_bdt_reader.SetVariable("min_dR","photon_drmin[0]");
  std_bdt_reader.SetVariable("max_dR",photon_drmax);
  std_bdt_reader.SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
  std_bdt_reader.SetVariable("cosTheta","llphoton_cosTheta[0]");
  std_bdt_reader.SetVariable("costheta","llphoton_costheta[0]");
  std_bdt_reader.SetVariable("phi","llphoton_psi[0]");
  std_bdt_reader.SetVariable("photon_res","photon_pterr[0]/photon_pt[0]");
  std_bdt_reader.SetVariable("photon_rapidity","photon_eta[0]");
  std_bdt_reader.SetVariable("l1_rapidity",lead_lepton_eta);
  std_bdt_reader.SetVariable("l2_rapidity",sublead_lepton_eta);
  std_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_run2_BDT.weights.xml");
  NamedFunc std_bdt_score = std_bdt_reader.GetDiscriminant();

  //--------------------------------------------------------------------------
  //                         Make plots & tables
  //--------------------------------------------------------------------------

  PlotMaker pm;

  bool plot_pts = false;
  bool plot_baseline = false;
  bool plot_bdtvars = false;
  bool plot_densitystudies = true;

  if (plot_pts) {
    pm.Push<Hist1D>(Axis(40,0,80, el1_pt, "Lead Electron p_{T} [GeV]", {25}), 
        "ll_lepid[0]==11", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(40,0,80, el2_pt, "Sublead Electron p_{T} [GeV]", {15}), 
        "ll_lepid[0]==11", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(40,0,80, mu1_pt, "Lead Muon p_{T} [GeV]", {20}), 
        "ll_lepid[0]==13", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(40,0,80, mu2_pt, "Sublead Muon p_{T} [GeV]", {10}), 
        "ll_lepid[0]==13", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
  }

  if (plot_baseline) {
    NamedFunc baseline_nomll = (zg_el_cuts||zg_mu_cuts)&&
        "(photon_pt[0]/llphoton_m[0]>=15.0/110.0)&&photon_drmin[0]>0.4&&photon_idmva[0]>0.5";
    baseline_nomll.Name("baseline_nomll");
    NamedFunc baseline_noptg = (zg_el_cuts||zg_mu_cuts)&&
        "(ll_m[0]>50) && ((llphoton_m[0]+ll_m[0])>=185) && photon_drmin[0]>0.4 && photon_idmva[0]>0.5";
    baseline_noptg.Name("baseline_noptg");
    NamedFunc ptg_over_mllg = "photon_pt[0]/llphoton_m[0]";
    //mll plots
    pm.Push<Hist2D>(Axis(20, 40, 120, "ll_m[0]", "m_{ll} [GeV]", {}),
        Axis(20, 80, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        baseline_nomll, procs_bak, ops2d).Weight(weight).Tag("zgggh_bak").LuminosityTag(lumi_string);
    pm.Push<Hist2D>(Axis(20, 40, 120, "ll_m[0]", "m_{ll} [GeV]", {}),
        Axis(20, 80, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        baseline_nomll, procs_sig, ops2d).Weight(weight).Tag("zgggh_sig").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma} [GeV]"), 
        baseline_nomll, procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma} [GeV]"), 
        zg_baseline&&"photon_idmva[0]>0.5", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(20,50,130, "ll_m[0]", "m_{ll} [GeV]"), 
        baseline_nomll&&higgs_window, procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //photon pt plots
    pm.Push<Hist1D>(Axis(26,15,80, "photon_pt[0]", "Photon p_{T} [GeV]"), 
        baseline_noptg&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(26,0.115,0.667, ptg_over_mllg, "Photon p_{T}/m_{ll#gamma}"), 
        baseline_noptg&&higgs_window, procs, ops_sorb)
        .Weight("w_lumi"*w_years).Tag("zgggh").LuminosityTag(lumi_string);
  }

  if (plot_bdtvars) {
    NamedFunc photon_res = "photon_pterr[0]/photon_pt[0]";
    NamedFunc llphoton_pt_mass = "llphoton_pt[0]/llphoton_m[0]";
    NamedFunc ll_pt_mass = "ll_pt[0]/llphoton_m[0]";
    //
    pm.Push<Hist1D>(Axis(30,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA"), 
        zg_baseline&&higgs_window, procs, ops_sorb)
        .Weight("w_lumi"*w_years).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.0,0.5, photon_res, "Photon p_{T} relative uncertainty"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,0.4,3.5, "photon_drmin[0]", "Photon minimum #Delta R(l,#gamma)"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.4,5.0, photon_drmax, "Photon maximum #Delta R(l,#gamma)"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.0,4.0, llphoton_pt_mass, "p_{Tll#gamma}/m_{ll#gamma}"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,-1.0,1.0, llphoton_cosTheta, "cos #Theta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,-2.5,2.5, "photon_eta[0]", "Photon #eta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,-2.5,2.5, lep1_eta, "Lead lepton #eta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,-2.5,2.5, lep2_eta, "Sublead lepton #eta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,-1.0,1.0, "llphoton_costheta[0]", "cos #theta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,-3.14,3.14, "llphoton_psi[0]", "#phi"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,0.0,0.667, ll_pt_mass, "p_{Tll}/m_{ll#gamma}"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.0,5.0, "llphoton_dr[0]", "#Delta R(ll,#gamma)"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,0.0,80.0, "met", "p_{T}^{miss} [GeV]"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.0,80.0, "met", "p_{T}^{miss} [GeV]"), 
        zg_baseline&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
  }

  if (plot_densitystudies) {
    pm.multithreaded_ = false;
    //pm.Push<Hist1D>(Axis(30,0.0,225.0, "llphoton_pt[0]", "ll#gamma p_{T} [GeV]"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(30,-6.0,6.0, "llphoton_eta[0]", "ll#gamma #eta"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(30,-1.0,1.0, cos_cap_theta, "cos(#Theta)"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(30,-3.1416,3.1416, cap_phi, "#Phi"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(30,-1.0,1.0, cos_lower_theta, "cos(#theta)"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(30,-3.1416,3.1416, lower_phi, "#phi"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(30,0.14,1.0, "photon_idmva[0]", "Photon IDMVA"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(60,-0.4,0.4, bdt_score, "BDT Score"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(60,-0.6,0.9, std_bdt_score, "BDT Score"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //pm.Push<Hist1D>(Axis(60,0.0,1.0, dnn_likelihood_ratio, "1/(1+(p_{DNN}^{sig}/p_{DNN}^{bak})^{-1})"), 
    //    tighter_baseline&&higgs_window, procs_loaded, ops_shapes)
    //    .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
  }

  pm.min_print_ = true;
  pm.MakePlots(1.0);

  return 0;
}
