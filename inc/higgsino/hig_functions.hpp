#ifndef H_HIG_FUNCTIONS
#define H_HIG_FUNCTIONS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"
 
namespace Higfuncs{
  // count number of jets associated with true B-hadrons 
  // among the 4 jets used to build average mjj
  extern const NamedFunc ntrub;
  // nbfake = category minus ntrub
  extern const NamedFunc nbfake;

  //use both the MET-binned and HT-binned ttbar
  extern const NamedFunc stitch_htmet;
  
  // pT of leading Higgs bosons in MC
  extern const NamedFunc hig_pt1;
  extern const NamedFunc hig_pt2;

  // nominal and control region b-tag categorization
  extern const NamedFunc hig_bcat;
  extern const NamedFunc higd_bcat;
  extern const NamedFunc higd_bcat_extended;

  // alternative b-tag category options
  extern const NamedFunc higd_bcat_mmmm;
  extern const NamedFunc higd_bcat_ttll;
  extern const NamedFunc higd_bcat_tmml;

  // weight that allows subtracting the reweighted ttbar from the data histogram
  //extern const NamedFunc wgt_subtr_ttx;
  extern const NamedFunc wgt_subtr_other;
  // namedfunc to apply the weights in wgt_nb_met to all bkg. processes
  extern const NamedFunc wgt_comp;

  // calculate effect of systematics calculated for each background 
  // in the data control regions on the total bkg. kappa
  extern const NamedFunc wgt_syst_ttx;
  extern const NamedFunc wgt_syst_vjets;
  extern const NamedFunc wgt_syst_qcd;

  // NamedFunc that returns weight with scaling for composition systematic
  // note this returns the weight scaling, not a variation
  extern const NamedFunc wgt_syst_comp;
  extern const NamedFunc wgt_syst_lowptleptrig;

  // analysis trigger and its efficiency
  //NamedFunc::ScalarType trig_hig_decision(const Baby &b);
  //extern const NamedFunc trig_hig;
  extern const NamedFunc eff_higtrig;
  extern const NamedFunc eff_higtrig_run2;
  extern const NamedFunc eff_higtrig_run2_v0;
  extern const NamedFunc err_higtrig;
  extern const NamedFunc weight_higd;
  extern const NamedFunc pico_weight_higd;

  extern const NamedFunc mhig;
  extern const NamedFunc nb_exci;
  extern const NamedFunc nb_gs;

  extern const NamedFunc h1b1_pt;
  extern const NamedFunc h1b2_pt;
  extern const NamedFunc h2b1_pt;
  extern const NamedFunc h2b2_pt;
  extern const NamedFunc h1b1_eta;
  extern const NamedFunc h1b2_eta;
  extern const NamedFunc h2b1_eta;
  extern const NamedFunc h2b2_eta;
  extern const NamedFunc h1b1_jetid;
  extern const NamedFunc h1b2_jetid;
  extern const NamedFunc h2b1_jetid;
  extern const NamedFunc h2b2_jetid;
  extern const NamedFunc h1_dr;
  extern const NamedFunc h2_dr;
  extern const NamedFunc h1_mass;
  extern const NamedFunc h2_mass;

  extern const NamedFunc lead_signal_lepton_pt;
  extern const NamedFunc lead_signal_muon_pt;
  extern const NamedFunc lead_signal_electron_pt;

  //filters
  extern const NamedFunc pass_ecalnoisejet;
  extern const NamedFunc pass_hemveto;
  extern const NamedFunc hem_weight;
  extern const NamedFunc pass_filters;
  extern const NamedFunc final_pass_filters;
  extern const NamedFunc final_ttbar_pass_filters;
  extern const NamedFunc final_zll_pass_filters;
  extern const NamedFunc final_qcd_pass_filters;

  // weights
  extern const NamedFunc w_years;
  extern const NamedFunc final_weight;
  extern const NamedFunc final_weight_notrgeff;

  //triggers
  extern const NamedFunc jet_trigger;
  extern const NamedFunc met_trigger;
  extern const NamedFunc el_trigger;
  extern const NamedFunc mu_trigger;
  extern const NamedFunc ht_trigger;
  extern const NamedFunc jetht_trigger;

  extern const NamedFunc met_trigger_v0;
  extern const NamedFunc el_trigger_v0;
  extern const NamedFunc mu_trigger_v0;

  //jetid modified variables
  extern const NamedFunc jetid_njet;
  extern const NamedFunc jetid_nb;
  extern const NamedFunc jetid_low_dphi_met;
  extern const NamedFunc jetid_hig_cand_dm;
  extern const NamedFunc jetid_hig_cand_am;
  extern const NamedFunc jetid_hig_cand_drmax;
}

#endif
