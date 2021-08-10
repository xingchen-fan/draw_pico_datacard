#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"
#include "TMath.h"
#include "Math/Vector4D.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

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

vector<unsigned> get_higgs_bbjet_indices_islep(vector<float> const & jet_m, vector<float> const & jet_deepcsv, vector<float> const & jet_pt, vector<float> const & jet_eta, vector<float> const & jet_phi, vector<bool> const & jet_islep) {
  //cout<<"event start"<<endl;
  // Reconstruct resolved
  // Sort by btag
  vector<pair<float, unsigned> > bjet;
  for (unsigned ijet = 0; ijet < jet_m.size(); ijet++) {
    // Filter jets
    if (jet_pt[ijet]<=30) continue;
    if (fabs(jet_eta[ijet])>2.4) continue;
    if (jet_islep[ijet]) continue;
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

extern const NamedFunc higgs_h1b1_pt("higgs_h1b1_pt", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return (*b.jet_pt())[bbjet_indices.at(0)];
});
extern const NamedFunc higgs_h1b2_pt("higgs_h1b2_pt", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return (*b.jet_pt())[bbjet_indices.at(1)];
});
extern const NamedFunc higgs_h2b1_pt("higgs_h2b1_pt", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return (*b.jet_pt())[bbjet_indices.at(2)];
});
extern const NamedFunc higgs_h2b2_pt("higgs_h2b2_pt", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return (*b.jet_pt())[bbjet_indices.at(3)];
});

const NamedFunc boost_low_dphi("boost_low_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  bool low_dphi = false;
  if ((*b.jet_met_dphi()).size()>3) {
    low_dphi = (*b.jet_met_dphi())[0]<0.5 || (*b.jet_met_dphi())[1]<0.5 ||
                          (*b.jet_met_dphi())[2]<0.3 || (*b.jet_met_dphi())[3]<0.3;
  }
  if ((*b.jet_met_dphi()).size()==3) {
    low_dphi = (*b.jet_met_dphi())[0]<0.5 || (*b.jet_met_dphi())[1]<0.5 ||
                          (*b.jet_met_dphi())[2]<0.3;
  }
  return low_dphi;
});
extern const NamedFunc higgs_h1_dr("higgs_h1_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  return DeltaR(h1b1, h1b2);
});
extern const NamedFunc higgs_h2_dr("higgs_h2_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  return DeltaR(h2b1, h2b2);
});
extern const NamedFunc higgs_min_dr("higgs_min_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  float h1_dr = DeltaR(h1b1, h1b2);
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  float h2_dr = DeltaR(h2b1, h2b2);
  return min(h1_dr, h2_dr);
});
extern const NamedFunc higgs_h1_mass("higgs_h1_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  ROOT::Math::PtEtaPhiMVector higgs1 = h1b1 + h1b2;
  return higgs1.M();
});
extern const NamedFunc higgs_h2_mass("higgs_h2_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  ROOT::Math::PtEtaPhiMVector higgs2 = h2b1 + h2b2;

  //cout<<"higg0 am: "<<(*b.hig_cand_am())[0]<<endl;
  //ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  //ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  //ROOT::Math::PtEtaPhiMVector h1 = h1b1 + h1b2;
  //ROOT::Math::PtEtaPhiMVector h2 = h2b1 + h2b2;
  //cout<<h1.M()<<" "<<h2.M()<<" "<<(h1.M()+h2.M())/2<<endl;
  return higgs2.M();
});
extern const NamedFunc higgs_average_mass("higgs_average_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  ROOT::Math::PtEtaPhiMVector higgs1 = h1b1 + h1b2;
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  ROOT::Math::PtEtaPhiMVector higgs2 = h2b1 + h2b2;
  return (higgs1.M()+higgs2.M())/2;
});

extern const NamedFunc higgs_h1b1_bcsv("higgs_h1b1_bcsv", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return (*b.jet_deepcsv())[bbjet_indices.at(0)];
});
extern const NamedFunc higgs_h1b2_bcsv("higgs_h1b2_bcsv", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return (*b.jet_deepcsv())[bbjet_indices.at(1)];
});
extern const NamedFunc higgs_h1_bcsv_diff("higgs_h1_bcsv_diff", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  float h1b1_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(0)];
  float h1b2_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(1)];
  return fabs(h1b1_bcsv - h1b2_bcsv);
});
extern const NamedFunc higgs_h2b1_bcsv("higgs_h2b1_bcsv", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return (*b.jet_deepcsv())[bbjet_indices.at(2)];
});
extern const NamedFunc higgs_h2b2_bcsv("higgs_h2b2_bcsv", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return (*b.jet_deepcsv())[bbjet_indices.at(3)];
});
extern const NamedFunc higgs_h2_bcsv_diff("higgs_h2_bcsv_diff", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  float h2b1_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(2)];
  float h2b2_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(3)];
  return fabs(h2b1_bcsv - h2b2_bcsv);
});
extern const NamedFunc higgs_btag_good("higgs_btag_good", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  float h1b1_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(0)];
  float h1b2_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(1)];
  float h2b1_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(2)];
  float h2b2_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(3)];
  return h1b1_bcsv>0.8953&&h1b2_bcsv>0.2217&&h2b1_bcsv>0.8953&&h2b2_bcsv>0.2217;
});
extern const NamedFunc higgs_h1b1_pflavor("higgs_h1b1_pflavor", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return abs((*b.jet_pflavor())[bbjet_indices.at(0)]);
});
extern const NamedFunc higgs_h1b2_pflavor("higgs_h1b2_pflavor", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return abs((*b.jet_pflavor())[bbjet_indices.at(1)]);
});
extern const NamedFunc higgs_h2b1_pflavor("higgs_h2b1_pflavor", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return abs((*b.jet_pflavor())[bbjet_indices.at(2)]);
});
extern const NamedFunc higgs_h2b2_pflavor("higgs_h2b2_pflavor", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return abs((*b.jet_pflavor())[bbjet_indices.at(3)]);
});

extern const NamedFunc higgs_h1b2_jet_mass("higgs_h1b2_jet_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  // Combine h1b2 with other jets
  vector<ROOT::Math::PtEtaPhiMVector> w_cand;
  w_cand.reserve((*b.jet_m()).size()-1);
  for (unsigned ijet = 0; ijet < (*b.jet_m()).size(); ijet++) {
    if (ijet == bbjet_indices.at(0)) continue;
    if (ijet == bbjet_indices.at(1)) continue;
    if (ijet == bbjet_indices.at(2)) continue;
    if (ijet == bbjet_indices.at(3)) continue;
    ROOT::Math::PtEtaPhiMVector jet ((*b.jet_pt())[ijet], (*b.jet_eta())[ijet], (*b.jet_phi())[ijet], (*b.jet_m())[ijet]);
    w_cand.push_back(jet+h1b2);
  }
  //// Select smallest R
  //vector<ROOT::Math::PtEtaPhiMVector> w_cand(1);
  //float smallestR = 10;
  //for (unsigned ijet = 0; ijet < (*b.jet_m()).size(); ijet++) {
  //  if (ijet == bbjet_indices.at(1)) continue;
  //  ROOT::Math::PtEtaPhiMVector jet ((*b.jet_pt())[ijet], (*b.jet_eta())[ijet], (*b.jet_phi())[ijet], (*b.jet_m())[ijet]);
  //  if (smallestR > DeltaR(jet, h1b2)) {
  //    smallestR = DeltaR(jet, h1b2);
  //    w_cand[0] = jet+h1b2;
  //  }
  //}
  // Return mass
  vector<double> w_mass;
  w_mass.reserve((*b.jet_m()).size()-1);
  if (w_cand.size()==0) return 9999;
  else return w_cand[0].M();
});
extern const NamedFunc higgs_h1b1_jet_mass("higgs_h1b1_jet_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  int bbjet_index = 0;
  ROOT::Math::PtEtaPhiMVector bjet ((*b.jet_pt())[bbjet_indices.at(bbjet_index)], (*b.jet_eta())[bbjet_indices.at(bbjet_index)], (*b.jet_phi())[bbjet_indices.at(bbjet_index)], (*b.jet_m())[bbjet_indices.at(bbjet_index)]);
  // Combine with other jets
  vector<ROOT::Math::PtEtaPhiMVector> w_cand;
  w_cand.reserve((*b.jet_m()).size()-1);
  for (unsigned ijet = 0; ijet < (*b.jet_m()).size(); ijet++) {
    if (ijet == bbjet_indices.at(bbjet_index)) continue;
    if (ijet == bbjet_indices.at(0)) continue;
    if (ijet == bbjet_indices.at(1)) continue;
    if (ijet == bbjet_indices.at(2)) continue;
    if (ijet == bbjet_indices.at(3)) continue;
    ROOT::Math::PtEtaPhiMVector jet ((*b.jet_pt())[ijet], (*b.jet_eta())[ijet], (*b.jet_phi())[ijet], (*b.jet_m())[ijet]);
    w_cand.push_back(jet+bjet);
  }
  // Return mass
  vector<double> w_mass;
  w_mass.reserve((*b.jet_m()).size()-1);
  for (unsigned iw = 0; iw < w_cand.size(); ++iw) {
    w_mass.push_back(w_cand[iw].M());
  }
  return w_mass[0];
});
extern const NamedFunc higgs_h2b1_jet_mass("higgs_h2b1_jet_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  int bbjet_index = 2;
  ROOT::Math::PtEtaPhiMVector bjet ((*b.jet_pt())[bbjet_indices.at(bbjet_index)], (*b.jet_eta())[bbjet_indices.at(bbjet_index)], (*b.jet_phi())[bbjet_indices.at(bbjet_index)], (*b.jet_m())[bbjet_indices.at(bbjet_index)]);
  // Combine with other jets
  vector<ROOT::Math::PtEtaPhiMVector> w_cand;
  w_cand.reserve((*b.jet_m()).size()-1);
  for (unsigned ijet = 0; ijet < (*b.jet_m()).size(); ijet++) {
    if (ijet == bbjet_indices.at(bbjet_index)) continue;
    if (ijet == bbjet_indices.at(0)) continue;
    if (ijet == bbjet_indices.at(1)) continue;
    if (ijet == bbjet_indices.at(2)) continue;
    if (ijet == bbjet_indices.at(3)) continue;
    ROOT::Math::PtEtaPhiMVector jet ((*b.jet_pt())[ijet], (*b.jet_eta())[ijet], (*b.jet_phi())[ijet], (*b.jet_m())[ijet]);
    w_cand.push_back(jet+bjet);
  }
  // Return mass
  vector<double> w_mass;
  w_mass.reserve((*b.jet_m()).size()-1);
  for (unsigned iw = 0; iw < w_cand.size(); ++iw) {
    w_mass.push_back(w_cand[iw].M());
  }
  return w_mass[0];
});
extern const NamedFunc higgs_h2b2_jet_mass("higgs_h2b2_jet_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  int bbjet_index = 3;
  ROOT::Math::PtEtaPhiMVector bjet ((*b.jet_pt())[bbjet_indices.at(bbjet_index)], (*b.jet_eta())[bbjet_indices.at(bbjet_index)], (*b.jet_phi())[bbjet_indices.at(bbjet_index)], (*b.jet_m())[bbjet_indices.at(bbjet_index)]);
  // Combine with other jets
  vector<ROOT::Math::PtEtaPhiMVector> w_cand;
  w_cand.reserve((*b.jet_m()).size()-1);
  for (unsigned ijet = 0; ijet < (*b.jet_m()).size(); ijet++) {
    if (ijet == bbjet_indices.at(bbjet_index)) continue;
    if (ijet == bbjet_indices.at(0)) continue;
    if (ijet == bbjet_indices.at(1)) continue;
    if (ijet == bbjet_indices.at(2)) continue;
    if (ijet == bbjet_indices.at(3)) continue;
    ROOT::Math::PtEtaPhiMVector jet ((*b.jet_pt())[ijet], (*b.jet_eta())[ijet], (*b.jet_phi())[ijet], (*b.jet_m())[ijet]);
    w_cand.push_back(jet+bjet);
  }
  // Return mass
  vector<double> w_mass;
  w_mass.reserve((*b.jet_m()).size()-1);
  for (unsigned iw = 0; iw < w_cand.size(); ++iw) {
    w_mass.push_back(w_cand[iw].M());
  }
  return w_mass[0];
});

const NamedFunc min_jet_dphi("min_jet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  for (size_t ijet = 0; ijet < (*b.jet_phi()).size(); ++ijet) {
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.jet_phi())[ijet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});

extern const NamedFunc jet0_p("jet0_p", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.jet_phi()->size() < 1) return -9999;
  else return ROOT::Math::PtEtaPhiMVector ((*b.jet_pt())[0], (*b.jet_eta())[0], (*b.jet_phi())[0], (*b.jet_m())[0]).P();
});

extern const NamedFunc average_btag("average_btag", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
  return ((*b.jet_deepcsv())[bbjet_indices.at(0)]+(*b.jet_deepcsv())[bbjet_indices.at(1)]+(*b.jet_deepcsv())[bbjet_indices.at(2)]+(*b.jet_deepcsv())[bbjet_indices.at(3)])/4;
});

const NamedFunc basecut("basecut", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.met()>150&&"ntk==0&&!low_dphi_met&&nvlep==0";
});

// TODO debug
vector<int> getJetTypes(vector<float> const & jet_deepcsv, vector<bool> const & jet_islep,
                         vector<float> const & jet_pt, vector<float> const & jet_eta, vector<float> const & jet_phi, vector<float> const & jet_m, 
                         vector<int> const & mc_id, vector<int> const & mc_mom, vector<int> const & mc_momidx,
                         vector<float> const & mc_pt, vector<float> const & mc_eta, vector<float> const & mc_phi,vector<float> const & mc_mass
                        ) {
  // Find jets that are used to reconstruct higgs
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(jet_m, jet_deepcsv, jet_pt, jet_eta, jet_phi, jet_islep);
  // TODO search for particles from tau
  // Find (duscb) quarks that are from ttbar's W and b quark from ttbar
  // ID: 0: unknown, (1:d, 2:u, 3:s, 4:c, 5:b) from W or tau, 6:b from t
  // ttbarQuarks = [(quark id, quark 4vector)]
  vector<tuple<int, ROOT::Math::PtEtaPhiMVector> > ttbarQuarks; 
  for (unsigned iMc = 0; iMc < (mc_pt).size(); ++iMc) {
    //int grandmom_org = (mc_mom)[(mc_momidx)[iMc]];
    //cout<<"id: "<<mc_id[iMc]<<" mom_id: "<<(mc_mom)[iMc]<<" grandmom_id: "<<grandmom_org<<endl;
    // b quark from ttbar
    if (abs((mc_mom)[iMc])==6 && abs((mc_id)[iMc])==5) {
      ttbarQuarks.push_back({6, ROOT::Math::PtEtaPhiMVector((mc_pt)[iMc], (mc_eta)[iMc], (mc_phi)[iMc], (mc_mass)[iMc])});
      //cout<<"store from t"<<endl;
    }
    // quarks from W
    if (abs((mc_mom)[iMc])==24) {
      int grandmom = (mc_mom)[(mc_momidx)[iMc]];
      if (abs(grandmom) == 6) {
        //cout<<"store from w"<<endl;
        ttbarQuarks.push_back({abs((mc_id)[iMc]), ROOT::Math::PtEtaPhiMVector((mc_pt)[iMc], (mc_eta)[iMc], (mc_phi)[iMc], (mc_mass)[iMc])});
      }
    }
  }
  //// Print quarks from ttbar
  //for (unsigned iQuark = 0; iQuark < ttbarQuarks.size(); ++iQuark) {
  //  cout<<"id: "<<get<0>(ttbarQuarks[iQuark])<<" pt: "<<get<1>(ttbarQuarks[iQuark]).Pt()<<endl;
  //}

  //// Find closest quark to jet
  //// closestQuarkToJet = [(id, deltaR)]
  //vector<tuple<int,int> > closestQuarkToJet;
  //for (unsigned hJet = 0; hJet < bbjet_indices.size(); ++hJet) {
  //  int ijet = bbjet_indices[hJet];
  //  ROOT::Math::PtEtaPhiMVector jet ((jet_pt)[ijet], (jet_eta)[ijet], (jet_phi)[ijet], (jet_m)[ijet]);
  //  float minR = 1000;
  //  int minQuarkId = -1;
  //  for (unsigned iQuark = 0; iQuark < ttbarQuarks.size(); ++iQuark) {
  //    int quarkId = get<0>(ttbarQuarks[iQuark]);
  //    ROOT::Math::PtEtaPhiMVector & quark = get<1>(ttbarQuarks[iQuark]);
  //    float deltaR = DeltaR(jet, quark);
  //    if (minR>deltaR) {
  //      minR = deltaR;
  //      minQuarkId = quarkId;
  //    }
  //  }
  //  closestQuarkToJet.push_back({minQuarkId, minR});
  //}
  //vector<int> quarks;
  //for (unsigned iJet = 0; iJet < closestQuarkToJet.size(); ++iJet) {
  //  int quarkId = get<0>(closestQuarkToJet[iJet]);
  //  quarks.push_back(quarkId);
  //};

  // Match quarks to jets
  //matchQuarkToJet[quark index] = (jet index, deltaR)
  map<int, tuple<int, float> > matchQuarkToJet;
  for (unsigned iQuark = 0; iQuark < ttbarQuarks.size(); ++iQuark) {
    //int quarkId = get<0>(ttbarQuarks[iQuark]);
    ROOT::Math::PtEtaPhiMVector const & quark = get<1>(ttbarQuarks[iQuark]);
    int minJetIndex = -1;
    float minR = 1000;
    for (unsigned hJet = 0; hJet < bbjet_indices.size(); ++hJet) {
      int ijet = bbjet_indices[hJet];
      ROOT::Math::PtEtaPhiMVector jet ((jet_pt)[ijet], (jet_eta)[ijet], (jet_phi)[ijet], (jet_m)[ijet]);
      float deltaR = DeltaR(jet, quark);
      if (minR>deltaR) {
        minR = deltaR;
        minJetIndex = hJet;
      }
    }
    matchQuarkToJet[iQuark] = {minJetIndex, minR};
  }
  vector<int> quarks(4,-1);
  for (auto it = matchQuarkToJet.begin(); it != matchQuarkToJet.end(); ++it) {
    int quarkIndex = (*it).first;
    int quarkId = get<0>(ttbarQuarks[quarkIndex]);
    int jetIndex = get<0>((*it).second);
    float minR = get<1>((*it).second);
    if (minR>0.4) continue;
    quarks[jetIndex] = quarkId;
  }
  //cout<<quarks[0]<<" "<<quarks[1]<<" "<<quarks[2]<<" "<<quarks[3]<<endl;

  return quarks;
}

const NamedFunc jetType_h1b1("jetType_h1b1", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<int> jetTypes = getJetTypes(*b.jet_deepcsv(), *b.jet_islep(),
               *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_m(), 
               *b.mc_id(), *b.mc_mom(), *b.mc_momidx(),
               *b.mc_pt(), *b.mc_eta(), *b.mc_phi(), *b.mc_mass());
  return jetTypes[0];
});
const NamedFunc jetType_h1b2("jetType_h1b2", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<int> jetTypes = getJetTypes(*b.jet_deepcsv(), *b.jet_islep(),
               *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_m(), 
               *b.mc_id(), *b.mc_mom(), *b.mc_momidx(),
               *b.mc_pt(), *b.mc_eta(), *b.mc_phi(), *b.mc_mass());
  return jetTypes[1];
});
const NamedFunc jetType_h2b1("jetType_h2b1", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<int> jetTypes = getJetTypes(*b.jet_deepcsv(), *b.jet_islep(),
               *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_m(), 
               *b.mc_id(), *b.mc_mom(), *b.mc_momidx(),
               *b.mc_pt(), *b.mc_eta(), *b.mc_phi(), *b.mc_mass());
  return jetTypes[2];
});
const NamedFunc jetType_h2b2("jetType_h2b2", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<int> jetTypes = getJetTypes(*b.jet_deepcsv(), *b.jet_islep(),
               *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_m(), 
               *b.mc_id(), *b.mc_mom(), *b.mc_momidx(),
               *b.mc_pt(), *b.mc_eta(), *b.mc_phi(), *b.mc_mass());
  return jetTypes[3];
});

//const NamedFunc eventType("eventType", [](const Baby &b) -> NamedFunc::VectorType{
//  // Find jets that are used to reconstruct higgs
//  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices_islep(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_islep());
//  // TODO search for particles from tau
//  // Find (duscb) quarks that are from ttbar's W and b quark from ttbar
//  // ID: 0: unknown, (1:d, 2:u, 3:s, 4:c, 5:b) from W or tau, 6:b from t
//  // ttbarQuarks = [(quark id, quark 4vector)]
//  vector<tuple<int, ROOT::Math::PtEtaPhiMVector> > ttbarQuarks; 
//  for (unsigned iMc = 0; iMc < (*b.mc_pt()).size(); ++iMc) {
//    // b quark from ttbar
//    if (abs((*b.mc_mom())[iMc])==6 && abs((*b.mc_id())[iMc])==5) {
//      ttbarQuarks.push_back({6, ROOT::Math::PtEtaPhiMVector((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], (*b.mc_mass())[iMc])});
//    }
//    // quarks from W
//    if (abs((*b.mc_mom())[iMc])==24) {
//      int grandmom = (*b.mc_mom())[(*b.mc_momidx())[iMc]];
//      if (abs(grandmom) == 6) {
//        ttbarQuarks.push_back({abs((*b.mc_id())[iMc]), ROOT::Math::PtEtaPhiMVector((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], (*b.mc_mass())[iMc])});
//      }
//    }
//  }
//  // Find closest quark to jet
//  // closestQuarkToJet = [(id, deltaR)]
//  vector<tuple<int,int> > closestQuarkToJet;
//  for (unsigned hJet = 0; hJet < bbjet_indices.size(); ++hJet) {
//    int ijet = bbjet_indices[hJet];
//    ROOT::Math::PtEtaPhiMVector jet ((*b.jet_pt())[ijet], (*b.jet_eta())[ijet], (*b.jet_phi())[ijet], (*b.jet_m())[ijet]);
//    float minR = 1000;
//    int minQuarkId = -1;
//    for (unsigned iQuark = 0; iQuark < ttbarQuarks.size(); ++iQuark) {
//      int quarkId = get<0>(ttbarQuarks[iQuark]);
//      ROOT::Math::PtEtaPhiMVector & quark = get<1>(ttbarQuarks[iQuark]);
//      float deltaR = DeltaR(jet, quark);
//      if (minR>deltaR) {
//        minR = deltaR;
//        minQuarkId = quarkId;
//      }
//    }
//    closestQuarkToJet.push_back({minQuarkId, minR});
//  }
//  // Tag event (xxxx) MSB is h1b1, LSB is h2b2. x is 0: unknown, 1:d, 2:u, 3:s, 4:c, 5:b, 6:b from t
//  //int eventTag = 0;
//  //for (unsigned iJet = 0; iJet < closestQuarkToJet.size(); ++iJet) {
//  //  int quarkId = get<0>(closestQuarkToJet[iJet]);
//  //  float minR = get<1>(closestQuarkToJet[iJet]);
//  //  if (minR>0.4) quarkId = 0;
//  //  eventTag += quarkId * pow(10,(closestQuarkToJet.size()-iJet-1));
//  //}
//  //return eventTag;
//
//  vector<double> quarks;
//  for (unsigned iJet = 0; iJet < closestQuarkToJet.size(); ++iJet) {
//    int quarkId = get<0>(closestQuarkToJet[iJet]);
//    quarks.push_back(quarkId);
//  };
//  return quarks;
//});

const NamedFunc nhhjets_true("nhhjets_true", [](const Baby &b) -> NamedFunc::ScalarType{
  // Find mc b quarks
  vector<int> bIndices;
  for (unsigned iMc=0; iMc<b.mc_pt()->size(); iMc++){
    if (abs(b.mc_mom()->at(iMc))==25&&abs(b.mc_id()->at(iMc))==5) bIndices.push_back(iMc);
  }
  // Count jets near b quarks
  int nJetsCloseTob = 0;
  for (unsigned iJet = 0; iJet < b.jet_pt()->size(); iJet++) {
    bool mcCloseToJet = false;
    // Jet cut
    if ((*b.jet_pt())[iJet]<=30) continue;
    if (fabs((*b.jet_eta())[iJet])>2.4) continue;
    // Find number of jets close to b quarks
    ROOT::Math::PtEtaPhiMVector jet ((*b.jet_pt())[iJet], (*b.jet_eta())[iJet], (*b.jet_phi())[iJet], (*b.jet_m())[iJet]);
    for (unsigned bIndex = 0; bIndex < bIndices.size(); ++bIndex) {
      int iMc = bIndices[bIndex];
      ROOT::Math::PtEtaPhiMVector bQuark ((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], 4);
      if (DeltaR(jet, bQuark)<0.4) mcCloseToJet = true;
    }
    if (mcCloseToJet) nJetsCloseTob++;
  }
  return nJetsCloseTob;
});

const NamedFunc nbbclose_true("nbbclose_true", [](const Baby &b) -> NamedFunc::ScalarType{
  // Find mc b quarks
  vector<int> bIndices;
  for (unsigned iMc=0; iMc<b.mc_pt()->size(); iMc++){
    if (abs(b.mc_mom()->at(iMc))==25&&abs(b.mc_id()->at(iMc))==5) bIndices.push_back(iMc);
  }
  int nMcCloseMc = 0;
  // Count b quarks near b quarks
  for (unsigned bIndex = 0; bIndex < bIndices.size(); ++bIndex) {
    bool mcCloseToMc = false;
    int iMc = bIndices[bIndex];
    ROOT::Math::PtEtaPhiMVector bQuark ((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], 4);
    // Find number of bquarks close to b quarks
    for (unsigned b2Index = 0; b2Index < bIndices.size(); ++b2Index) {
      if (b2Index == bIndex) continue;
      int iMc2 = bIndices[b2Index];
      ROOT::Math::PtEtaPhiMVector b2Quark ((*b.mc_pt())[iMc2], (*b.mc_eta())[iMc2], (*b.mc_phi())[iMc2], 4);
      if (DeltaR(b2Quark, bQuark)<0.8) mcCloseToMc = true;
    }
    if (mcCloseToMc) nMcCloseMc++;
  }
  return nMcCloseMc;
});

const NamedFunc nbout_true("nbout_true", [](const Baby &b) -> NamedFunc::ScalarType{
  // Find mc b quarks
  vector<int> bIndices;
  for (unsigned iMc=0; iMc<b.mc_pt()->size(); iMc++){
    if (abs(b.mc_mom()->at(iMc))==25&&abs(b.mc_id()->at(iMc))==5) bIndices.push_back(iMc);
  }
  int nbout = 0;
  for (unsigned bIndex = 0; bIndex < bIndices.size(); ++bIndex) {
    int iMc = bIndices[bIndex];
    ROOT::Math::PtEtaPhiMVector bQuark ((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], 4);
    if (bQuark.Pt()<=30 || fabs(bQuark.Eta())>2.4) nbout++;
  }
  return nbout;
});

const NamedFunc nb_ptout_true("nb_ptout_true", [](const Baby &b) -> NamedFunc::ScalarType{
  // Find mc b quarks
  vector<int> bIndices;
  for (unsigned iMc=0; iMc<b.mc_pt()->size(); iMc++){
    if (abs(b.mc_mom()->at(iMc))==25&&abs(b.mc_id()->at(iMc))==5) bIndices.push_back(iMc);
  }
  int nbout = 0;
  for (unsigned bIndex = 0; bIndex < bIndices.size(); ++bIndex) {
    int iMc = bIndices[bIndex];
    ROOT::Math::PtEtaPhiMVector bQuark ((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], 4);
    if (bQuark.Pt()<=30) nbout++;
  }
  return nbout;
});

const NamedFunc nb_etaout_true("nb_etaout_true", [](const Baby &b) -> NamedFunc::ScalarType{
  // Find mc b quarks
  vector<int> bIndices;
  for (unsigned iMc=0; iMc<b.mc_pt()->size(); iMc++){
    if (abs(b.mc_mom()->at(iMc))==25&&abs(b.mc_id()->at(iMc))==5) bIndices.push_back(iMc);
  }
  int nbout = 0;
  for (unsigned bIndex = 0; bIndex < bIndices.size(); ++bIndex) {
    int iMc = bIndices[bIndex];
    ROOT::Math::PtEtaPhiMVector bQuark ((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], 4);
    if (fabs(bQuark.Eta())>2.4) nbout++;
  }
  return nbout;
});

const NamedFunc b_pt("b_pt", [](const Baby &b) -> NamedFunc::VectorType{
  // Find mc b quarks
  vector<int> bIndices;
  for (unsigned iMc=0; iMc<b.mc_pt()->size(); iMc++){
    if (abs(b.mc_mom()->at(iMc))==25&&abs(b.mc_id()->at(iMc))==5) bIndices.push_back(iMc);
  }
  vector<double> pts;
  for (unsigned bIndex = 0; bIndex < bIndices.size(); ++bIndex) {
    int iMc = bIndices[bIndex];
    ROOT::Math::PtEtaPhiMVector bQuark ((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], 4);
    pts.push_back(bQuark.Pt());
  }
  return pts;
});

const NamedFunc b_eta("b_eta", [](const Baby &b) -> NamedFunc::VectorType{
  // Find mc b quarks
  vector<int> bIndices;
  for (unsigned iMc=0; iMc<b.mc_pt()->size(); iMc++){
    if (abs(b.mc_mom()->at(iMc))==25&&abs(b.mc_id()->at(iMc))==5) bIndices.push_back(iMc);
  }
  vector<double> pts;
  for (unsigned bIndex = 0; bIndex < bIndices.size(); ++bIndex) {
    int iMc = bIndices[bIndex];
    ROOT::Math::PtEtaPhiMVector bQuark ((*b.mc_pt())[iMc], (*b.mc_eta())[iMc], (*b.mc_phi())[iMc], 4);
    pts.push_back(bQuark.Eta());
  }
  return pts;
});
  
namespace{
  bool single_thread = true;
  //bool do_twiki = true;
  double luminosity = 1.;
  //double luminosity = 137;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  PlotOpt style("txt/plot_styles.txt", "Scatter");
  //vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};
  vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};

  set<int> years;
  //years = {2016, 2017, 2018};
  years = {2016};
  //years = {2017};

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if((Contains(hostname, "cms") || Contains(hostname, "compute-")))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  //string mc_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  //string mc_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_humboldt/";
  string mc_base_folder = "/net/cms25/cms25r5/jbkim/pico/NanoAODv5/higgsino_humboldt/";
  // TODO
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  //string mc_skim_folder = "mc/unskimmed/";
  //string mc_skim_folder = "mc/merged_higmc_met150/";
  string sig_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";

  NamedFunc base_filters = HigUtilities::pass_2016; //since pass_fsjets is not quite usable...
  //NamedFunc wgt = "w_lumi*w_isr"*w_years*Higfuncs::eff_higtrig;
  //NamedFunc wgt = "w_lumi*w_isr*137."*Higfuncs::eff_higtrig;
  NamedFunc wgt = "w_lumi*w_isr"*Higfuncs::eff_higtrig;
  if (years.size()==1 && *years.begin()==2016) wgt *= "137.";
  else wgt *= w_years;

  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_SingleLept*",
                                  "*TTJets_DiLept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                  "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"
                                 });
  //mctags["ttonly"]     = set<string>({"*TTJets_*Lept*"});
  //mctags["topx"]     = set<string>({"*_TTZ*.root", "*_TTW*.root",
  //                                   "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  //mctags["qcd"]     = set<string>({
  //                                 "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               //"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               //"*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("tt", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),base_filters&&"stitch")); 

  //procs.push_back(Process::MakeShared<Baby_pico>("tt+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["topx"]),base_filters&&"stitch"));
  //procs.push_back(Process::MakeShared<Baby_pico>("ttonly", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["ttonly"]),base_filters&&"stitch"));

  bool plot_btag_comparison = false;
  bool plot_btag_relation = false;
  bool plot_btag_category = false;
  bool plot_ABCD = false;
  bool plot_dphi = false;

  //vector<string> sigm = {"200", "600", "950"}; 
  vector<string> sigm = {"225", "400", "700"}; 
  //vector<string> sigm = {"950"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need sigm.size() >= sig_colors.size()
  if (plot_dphi||plot_ABCD) {
    for (unsigned isig(0); isig<sigm.size(); isig++){
      procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
        sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), base_filters));
    }
  }

  string baseline = "ntk==0&&!low_dphi_met&&nvlep==0&&met>150";
  //string baseline = "ntk==0&&!low_dphi_met&&nvlep==0&&met>50";
  //string baseline = "!low_dphi_met&&nvlep==0&&met>150";
  //string baseline = "ntk==0&&!low_dphi_met&&nlep==2&&met>150";

  string baseline_fourjet = "njet>=4&&njet<=5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&nbm>=2"; 
  //NamedFunc baseline_fourjet = "njet>=4&&njet<=5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40"&&higgs_h1b2_pt>80&&higgs_h2b2_pt>80; 
  //NamedFunc baseline_fourjet = "njet>=4&&njet<=6&&hig_cand_am[0]<200&&hig_cand_dm[0]<60"; 
  //string baseline_fourjet = "njet>=4&&njet<=5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&nbm>=2&&hig_cand_am[0]<60"; 
  //string baseline_fourjet = "njet>=4&&njet<=5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&nbm>=2&&hig_cand_am[0]<150&&hig_cand_am[0]>100"; 
  //string baseline_fourjet = "njet==4&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&nbm>=2"; 
  //string baseline_fourjet = "njet==5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&nbm>=2"; 
  //string drcut = "hig_cand_drmax[0]>=0 && hig_cand_drmax[0]<2.2";

  //string drcut = "hig_cand_drmax[0]>=1.1 && hig_cand_drmax[0]<2.2";
  string drcut = "hig_cand_drmax[0]>=0 && hig_cand_drmax[0]<2.2";
  //string drcut = "hig_cand_drmax[0]>=0.8 && hig_cand_drmax[0]<1.2";

  // TTML
  string bcut_2b = "(nbt==2&&nbm==2)";
  string bcut_3b = "(nbt>=2&&nbm==3&&nbl==3)";
  string bcut_4b = "(nbt>=2&&nbm>=3&&nbl>=4)";
  //// MMMM
  //string bcut_2b = "(nbm==2)";
  //string bcut_3b = "(nbm==3)";
  //string bcut_4b = "(nbm>=4)";
  // TTMM
  //string bcut_2b = "(nbt==2&&nbm==2)";
  //string bcut_3b = "(nbt>=2&&nbm==3)";
  //string bcut_4b = "(nbt>=2&&nbm>=4)";
  //// TMMM
  //string bcut_2b = "(nbt>=1&&nbm==2)";
  //string bcut_3b = "(nbt>=1&&nbm==3)";
  //string bcut_4b = "(nbt>=1&&nbm>=4)";
  //// TMML
  //string bcut_2b = "(nbt>=1&&nbm==2)";
  //string bcut_3b = "(nbt>=1&&nbm==3&&nbl==3)";
  //string bcut_4b = "(nbt>=1&&nbm>=3&&nbl>=4)";
  //// TMLL ver1
  //string bcut_2b = "(nbt>=1&&nbm==2&&(nbl==2||nbl==3))";
  //string bcut_3b = "(nbt>=1&&nbm==3&&nbl==3)";
  //string bcut_4b = "(nbt>=1&&nbm>=2&&nbl>=4)";
  //// TMLL ver2
  //string bcut_2b = "(nbt>=1&&nbm>=2&&nbl==2)";
  //string bcut_3b = "(nbt>=1&&nbm>=2&&nbl==3)";
  //string bcut_4b = "(nbt>=1&&nbm>=2&&nbl>=4)";
  ////MMML
  //string bcut_2b = "(nbm==2)";
  //string bcut_3b = "(nbm==3&&nbl==3)";
  //string bcut_4b = "(nbm>=3&&nbl>=4)";
  ////LLLL (Actually MMLL because of baseline nbm>=2 cut)
  //string bcut_2b = "(nbl==2)";
  //string bcut_3b = "(nbl==3)";
  //string bcut_4b = "(nbl>=4)";
  ////TTTT
  //string bcut_2b = "(nbt==2)";
  //string bcut_3b = "(nbt==3)";
  //string bcut_4b = "(nbt>=4)";
  //// TTML
  //string bcut_2b = "(nbt==2&&nbm==2&&nbl>=4)";
  //string bcut_3b = "(nbt>=2&&nbm==3&&nbl==3)";
  //string bcut_4b = "(nbt>=2&&nbm>=3&&nbl>=4)";
  //// TTLL ver1
  //string bcut_2b = "(nbt==2&&nbm==2&&(nbl==2||nbl==3))";
  //string bcut_3b = "(nbt>=2&&nbm==3&&nbl==3)";
  //string bcut_4b = "(nbt>=2&&nbm>=2&&nbl>=4)";
  //// TTLL ver2
  //string bcut_2b = "(nbt>=2&&nbm>=2&&nbl==2)";
  //string bcut_3b = "(nbt>=2&&nbm>=2&&nbl==3)";
  //string bcut_4b = "(nbt>=2&&nbm>=2&&nbl>=4)";

  //baseline_fourjet += "&&("+bcut_3b+"||"+bcut_4b+")";
  //baseline_fourjet += "&&"+bcut_4b;
  //baseline_fourjet += "&&"+bcut_2b;

  NamedFunc extra_cut = "1";
  //NamedFunc extra_cut = higgs_h2b2_pt>80;
  //NamedFunc extra_cut = higgs_h1b2_pt>80&&higgs_h2b2_pt>80;
  //NamedFunc extra_cut = higgs_h2b2_pt>80;
  //NamedFunc extra_cut = higgs_h2b2_pt>80;
  //NamedFunc extra_cut = "hig_cand_am[0]>100&&hig_cand_am[0]<140";
  //NamedFunc extra_cut = "hig_cand_am[0]>100";
  //NamedFunc extra_cut = "hig_cand_am[0]<80";

  PlotMaker pm;
  bool plot_jet_study = true;
  if (plot_jet_study) {
    vector<shared_ptr<Process> > procs_signal;
    for (unsigned isig(0); isig<sigm.size(); isig++){
      procs_signal.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
        sig_colors[isig], attach_folder(sig_base_folder, years, "SMS-TChiHH_2D/unskimmed/", {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), "1"));
    }
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, "njet", "njets", {}),
      "1", procs_signal, plt_lin).Weight("weight*35.9");
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, nhhjets_true, "true hh jets", {}),
      "met>150", procs_signal, plt_lin).Weight("weight*35.9");
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, nbbclose_true, "Number of close b quarks", {}),
      "met>150", procs_signal, plt_lin).Weight("weight*35.9"*Higfuncs::eff_higtrig);
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, nbout_true, "Number b quarks out of acceptance", {}),
      "met>150", procs_signal, plt_lin).Weight("weight*35.9"*Higfuncs::eff_higtrig);
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, nb_ptout_true, "Number b quarks out of pt acceptance", {}),
      "met>150", procs_signal, plt_lin).Weight("weight*35.9"*Higfuncs::eff_higtrig);
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, nb_etaout_true, "Number b quarks out of eta acceptance", {}),
      "met>150", procs_signal, plt_lin).Weight("weight*35.9"*Higfuncs::eff_higtrig);
    pm.Push<Hist1D>(Axis(50, 0, 500, b_pt, "b quarks pt", {}),
      "met>150", procs_signal, plt_shapes_info).Weight("weight*35.9"*Higfuncs::eff_higtrig);
    pm.Push<Hist1D>(Axis(30, -3, 3, b_eta, "b quarks eta", {}),
      "met>150", procs_signal, plt_shapes_info).Weight("weight*35.9"*Higfuncs::eff_higtrig);
  }

  if (plot_btag_comparison) {
    // 2b, 3b, 4b comparison
    vector<shared_ptr<Process> > procs_btag;
    procs_btag.push_back(Process::MakeShared<Baby_pico>("2b", Process::Type::background,sig_colors[0],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut));
    procs_btag.push_back(Process::MakeShared<Baby_pico>("3b", Process::Type::background,sig_colors[1],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut));
    procs_btag.push_back(Process::MakeShared<Baby_pico>("4b", Process::Type::background,sig_colors[2],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut));
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      1, procs_btag, plt_shapes).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "subleading mass [GeV]", {100,140}),
      1, procs_btag, plt_shapes).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      1, procs_btag, plt_shapes).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "dr", {1.1,2.2}),
      1, procs_btag, plt_shapes).Weight(wgt);

    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);

    pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "subleading mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "subleading mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "subleading mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);

    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&bcut_3b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);

    //pm.Push<Hist1D>(Axis(22, 0, 2.2, higgs_h1_dr, "dr", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(22, 0, 2.2, higgs_h1_dr, "dr", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(22, 0, 2.2, higgs_h1_dr, "dr", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);

    pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "dr", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "dr", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "dr", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);

    //pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
  }

  if (plot_btag_relation) {
    // b tag relation
    pm.Push<Hist2D>(Axis(20, 0, 1, higgs_h1b1_bcsv, "higgs1 b1 btag", {0.2217, 0.6321, 0.8953}),
      Axis(15, 0, 1, higgs_h1b2_bcsv, "higgs1 b2 btag", {0.2217, 0.6321, 0.8953}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_2D).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h1_bcsv_diff, "Btag diff of leading Higg", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist2D>(Axis(20, 0, 1, higgs_h2b1_bcsv, "higgs2 b1 btag", {0.2217, 0.6321, 0.8953}),
      Axis(15, 0, 1, higgs_h2b2_bcsv, "higgs2 b2 btag", {0.2217, 0.6321, 0.8953}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_2D).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h2_bcsv_diff, "Btag diff of sub-leading Higg", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    // b tags
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h1b1_bcsv, "Leading pt Higg, leading Btag jet, Btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h1b2_bcsv, "Leading pt Higg, subleading Btag jet, Btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h2b1_bcsv, "SubLeading pt Higg, leading Btag jet, Btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h2b2_bcsv, "SubLeading pt Higg, subleading Btag jet, Btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
  }

  if (plot_btag_category) {
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&higgs_btag_good&&extra_cut, procs, plt_lin).Weight(wgt);
  }

  if (plot_ABCD) {

    //// 2b, 3b, 4b comparison
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);


    vector<shared_ptr<Process> > procs_btag;
    procs_btag.push_back(Process::MakeShared<Baby_pico>("2b", Process::Type::background,sig_colors[0],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut));
    procs_btag.push_back(Process::MakeShared<Baby_pico>("3b", Process::Type::background,sig_colors[1],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut));
    procs_btag.push_back(Process::MakeShared<Baby_pico>("4b", Process::Type::background,sig_colors[2],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut));

    // Observables
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "Average Higgs mass [GeV]", {100,140}),
      1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsAm");
    pm.Push<Hist1D>(Axis(15, 0, 300, "isr_tru_pt", "ISR true pT [GeV]", {}),
      1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsIsrTruPt");
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average Higgs mass [GeV]", {100,140}),
    //  higgs_h1b2_pflavor>=5&&higgs_h2b2_pflavor>=5, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsAm");
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h1b2, "jetType h1b2", {}),
      1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsh1b2");
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h2b2, "jetType h2b2", {}),
      1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsh2b2");
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, "nisr", "nisr", {}),
      1, procs_btag, plt_shapes).Weight(wgt).Tag("nisr");

    // Important variables
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h1b1, "jetType h1b1", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h1b2, "jetType h1b2", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h2b1, "jetType h2b1", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h2b2, "jetType h2b2", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    //// low mass
    //pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h1b1, "jetType h1b1", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut&&"hig_cand_am[0]<80", procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h1b2, "jetType h1b2", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut&&"hig_cand_am[0]<80", procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h2b1, "jetType h2b1", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut&&"hig_cand_am[0]<80", procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h2b2, "jetType h2b2", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut&&"hig_cand_am[0]<80", procs, plt_lin).Weight(wgt);
    //// high mass
    //pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h1b1, "jetType h1b1", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut&&"hig_cand_am[0]>100", procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h1b2, "jetType h1b2", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut&&"hig_cand_am[0]>100", procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h2b1, "jetType h2b1", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut&&"hig_cand_am[0]>100", procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, -0.5, 7.5, jetType_h2b2, "jetType h2b2", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut&&"hig_cand_am[0]>100", procs, plt_lin).Weight(wgt);

    pm.Push<Hist1D>(Axis(18, 150, 600, "met", "MET [GeV]", {200, 300, 450}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 0, 150, "hig_cand_dm[0]", "Diff mass [GeV]", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(10, 0, 10, "njet", "Number of jets", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(10, 0, 10, "ntk", "Number of iso tk", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 0, 300, higgs_h1b1_pt, "pT of Higgs 1, jet 1 [GeV]", {80}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 0, 300, higgs_h1b2_pt, "pT of Higgs 1, jet 2 [GeV]", {80}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 0, 300, higgs_h2b1_pt, "pT of Higgs 2, jet 1 [GeV]", {80}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 0, 300, higgs_h2b2_pt, "pT of Higgs 2, jet 2 [GeV]", {80}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h1b1_bcsv, "Btag of Higgs 1, jet 1", {0.2217, 0.6321, 0.8953}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h1b2_bcsv, "Btag of Higgs 1, jet 2", {0.2217, 0.6321, 0.8953}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h2b1_bcsv, "Btag of Higgs 2, jet 1", {0.2217, 0.6321, 0.8953}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h2b2_bcsv, "Btag of Higgs 2, jet 2", {0.2217, 0.6321, 0.8953}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 0, 300, "isr_tru_pt", "ISR true pT [GeV]", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);

    // lead sublead higgs
    //pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h1_mass, "Leading Higgs mass [GeV]", {100,140}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsLeadM");
    //pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "Subleading Higgs mass [GeV]", {100,140}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsSubM");
    //pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h1_dr, "Leading Higgs dr", {1.1,2.2}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsLeadM");
    //pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "Subleading Higgs dr", {1.1,2.2}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsSubM");

    //pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h1b1_pt, "pT of Higg 1, jet 1", {100./0.8,140./0.8}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("btag");
    //pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h1b2_pt, "pT of Higg 1, jet 2", {100./0.8,140./0.8}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("btag");
    //pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2b1_pt, "pT of Higg 2, jet 1", {100./0.9,140./0.9}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("btag");
    //pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2b2_pt, "pT of Higg 2, jet 2", {100./0.9,140./0.9}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("btag");

    //// TODO none
    // flavor
    //pm.Push<Hist1D>(Axis(22, 0, 21.5, higgs_h1b1_pflavor, "pflavor of Higg 1, jet 1", {}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("btag");
    //pm.Push<Hist1D>(Axis(22, 0, 21.5, higgs_h1b2_pflavor, "pflavor of Higg 1, jet 2", {}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("btag");
    //pm.Push<Hist1D>(Axis(22, 0, 21.5, higgs_h2b1_pflavor, "pflavor of Higg 2, jet 1", {}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("btag");
    //pm.Push<Hist1D>(Axis(22, 0, 21.5, higgs_h2b2_pflavor, "pflavor of Higg 2, jet 2", {}),
    //  1, procs_btag, plt_shapes).Weight(wgt).Tag("btag");

    //pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h1_dr, "Leading Higgs dr", {1.1}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "Subleading Higgs dr", {1.1}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h1_mass, "Leading Higgs mass [GeV]", {100, 140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "Subleading Higgs mass [GeV]", {100, 140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);

    //pm.Push<Hist1D>(Axis(100, 0, 500, higgs_h1b1_jet_mass, "h1b1+jet mass [GeV]", {80.379}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(100, 0, 500, higgs_h1b2_jet_mass, "h1b2+jet mass [GeV]", {80.379}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(100, 0, 500, higgs_h2b1_jet_mass, "h2b1+jet mass [GeV]", {80.379}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(100, 0, 500, higgs_h2b2_jet_mass, "h2b2+jet mass [GeV]", {80.379}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);


    //vector<shared_ptr<Process> > procs_h1b1_btag;
    //NamedFunc h1b1_b0 = higgs_h1b1_bcsv<=0.2217;
    //NamedFunc h1b1_bL = higgs_h1b1_bcsv>0.2217;
    //NamedFunc h1b1_bM = higgs_h1b1_bcsv>0.6321;
    //NamedFunc h1b1_bT = higgs_h1b1_bcsv<0.8953;
    //procs_h1b1_btag.push_back(Process::MakeShared<Baby_pico>("0b", Process::Type::background,sig_colors[0],
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
    //                base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&h1b1_b0&&extra_cut));
    //procs_h1b1_btag.push_back(Process::MakeShared<Baby_pico>("bL", Process::Type::background,sig_colors[1],
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
    //                base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&h1b1_bL&&extra_cut));
    //procs_h1b1_btag.push_back(Process::MakeShared<Baby_pico>("bM", Process::Type::background,sig_colors[2],
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
    //                base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&h1b1_bM&&extra_cut));
    //procs_h1b1_btag.push_back(Process::MakeShared<Baby_pico>("bT", Process::Type::background,sig_colors[3],
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
    //                base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&h1b1_bT&&extra_cut));
    //pm.Push<Hist1D>(Axis(20, 30, 400, higgs_h1b1_pt, "pT of Higg 1, jet 1", {}),
    //  1, procs_h1b1_btag, plt_shapes).Weight(wgt).Tag("h1b1_btag");


    //// Average btag vs average higgs mass: 
    //// - Doesn't show strong relation bewteen average btag and hig mass
    //pm.Push<Hist2D>(Axis(20, 0, 1, average_btag, "average_btag", {0.2217, 0.6321, 0.8953}),
    //  Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_2D).Weight(wgt);
    ////// Below two plots don't show much
    //////pm.Push<Hist1D>(Axis(20, 0, 1, average_btag, "average btag", {}),
    //////  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    //////pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
    //////  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);

    //// See how mass shape changes depending on average btag
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.2217&&average_btag<=0.6321, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.6321&&average_btag<=0.8953, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.8953, procs, plt_lin).Weight(wgt);
    //// See how mass shape changes depending on average btag
    //vector<shared_ptr<Process> > procs_average_btag;
    //procs_average_btag.push_back(Process::MakeShared<Baby_pico>("low b", Process::Type::background,sig_colors[0],
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
    //                base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.2217&&average_btag<=0.6321));
    //procs_average_btag.push_back(Process::MakeShared<Baby_pico>("mid b", Process::Type::background,sig_colors[1],
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
    //                base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.6321&&average_btag<=0.8953));
    //procs_average_btag.push_back(Process::MakeShared<Baby_pico>("high b", Process::Type::background,sig_colors[2],
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
    //                base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.8953));
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
    //  1, procs_average_btag, plt_shapes).Weight(wgt).Tag("btagVsAm");

    //// Average b tag vs dr max
    //pm.Push<Hist2D>(Axis(20, 0, 1, average_btag, "average_btag", {0.2217, 0.6321, 0.8953}),
    //  Axis(16, 0, 3.2, "hig_cand_drmax[0]", "Max dela R_{bb}", {1.1,2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_2D).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 0, 1, average_btag, "average btag", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    //// See how drmax changes depending on average btag
    //pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.2217&&average_btag<=0.6321, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.6321&&average_btag<=0.8953, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut&&average_btag>0.8953, procs, plt_lin).Weight(wgt);
    //// See how mass shape changes depending on average btag
    //pm.Push<Hist1D>(Axis(15, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
    //  1, procs_average_btag, plt_shapes).Weight(wgt).Tag("btagVsDrmax");
  }

  if (plot_dphi) {
    pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {0.5}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {0.5}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {0.5}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b&&extra_cut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {0.5}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b&&extra_cut, procs, plt_lin).Weight(wgt);
  }

  pm.multithreaded_ = true;
  //pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(luminosity);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(false){
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

vector<unsigned> higidx(const Baby &b){
  vector<unsigned> idx;
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_id()->at(i)==25) idx.push_back(i);
    if (idx.size()>1) break;
  }
  return idx;
}

// vector<unsigned> bidx(const Baby &b){
//   vector<unsigned> idx;
//   for (unsigned i(0); i<b.mc_pt()->size(); i++){
//     if (b.mc_id()->at(i)==25) higidx.push_back(i);
//     if (higidx.size()>1) break;
//   }
//   return idx;
// }
