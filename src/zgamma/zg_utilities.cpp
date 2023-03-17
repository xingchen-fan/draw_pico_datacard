#include <algorithm>
#include <stdlib.h>
#include <regex>
#include "zgamma/zg_utilities.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"
#include "core/baby.hpp"

namespace ZgUtilities {
  using std::string;
  using std::to_string;
  using std::cout;
  using std::endl;
  // Returns negative lepton 4-momentum
  TLorentzVector AssignL1(const Baby &b, bool gen) {
    TLorentzVector l1;
    if(gen) {
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == 11 || 
           b.mc_id()->at(i) == 13
           //|| b.mc_id()->at(i) == 15
           ) {
            if(b.mc_mom()->at(i) == 23) {
              l1.SetPtEtaPhiM(b.mc_pt()->at(i),
                              b.mc_eta()->at(i),
                              b.mc_phi()->at(i),
                              b.mc_mass()->at(i));
              break;
            }
          }
    }
    else{
      int il(-1);
      if(b.ll_lepid()->at(0) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
        else                                        il = b.ll_i2()->at(0);
        l1.SetPtEtaPhiM(b.el_pt() ->at(il), 
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.00511);
      }
      else if(b.ll_lepid()->at(0) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
        else                                        il = b.ll_i2()->at(0);
        l1.SetPtEtaPhiM(b.mu_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
      }
    }
    return l1;
  }

/*  TLorentzVector AssignJJ(const Baby &b, bool gen){
    TLorentzVector jj;
    double  jj_m = 100000;
    int r_idx = -1;

    for(size_t i; i < b.dijet_m()->size(); i++){
      if( fabs( jj_m - 91.2 ) < fabs( (b.dijet_m() -> at(i)) - 91.2)){
        jj_m = b.dijet_m()->at(i);
        r_idx = i;
      }

    }
   
   jj.SetPtEtaPhiM(b.dijet_pt() ->at(r_idx),
                   b.dijet_eta()->at(r_idx),
                   b.dijet_phi()->at(r_idx),
   return jj;    
  }
*/
  // Returns positive lepton 4-momentum
  TLorentzVector AssignL2(const Baby &b, bool gen) {
    TLorentzVector l2;
    if(gen) {
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == -11 || 
           b.mc_id()->at(i) == -13
           //|| b.mc_id()->at(i) == -15
           ) {
            if(b.mc_mom()->at(i) == 23) {
              //cout<<"AssignL2 Z index: "<<b.mc_momidx()->at(i)<<endl;
              //cout<<"AssignL2 Z mom index: "<<b.mc_momidx()->at(b.mc_momidx()->at(i))<<endl;
              //cout<<"AssignL2 Z mom pid: "<<b.mc_id()->at(b.mc_momidx()->at(b.mc_momidx()->at(i)))<<endl;
              //cout<<"AssignL2 Z mom pt: "<<b.mc_pt()->at(b.mc_momidx()->at(b.mc_momidx()->at(i)))<<endl;
              l2.SetPtEtaPhiM(b.mc_pt()->at(i),
                              b.mc_eta()->at(i),
                              b.mc_phi()->at(i),
                              b.mc_mass()->at(i));
              break;
            }
          }
    }
    else{
      int il(-1);
      if(b.ll_lepid()->at(0) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
        else                                        il = b.ll_i1()->at(0);
        l2.SetPtEtaPhiM(b.el_pt() ->at(il), 
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.00511);
      }
      else if(b.ll_lepid()->at(0) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
        else                                        il = b.ll_i1()->at(0);
        l2.SetPtEtaPhiM(b.mu_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
      }
    }
    return l2;
  }

  // Returns Z 4-momentum
  TLorentzVector AssignZ(const Baby &b, bool gen) {
    TLorentzVector ll;
    if(gen) 
        ll = AssignL1(b,gen) + AssignL2(b,gen);
    else
      ll.SetPtEtaPhiM(b.ll_pt() ->at(0),
                      b.ll_eta()->at(0),
                      b.ll_phi()->at(0),
                      b.ll_m()  ->at(0));
    return ll;
  }

  // Returns photon 4-momentum
  TLorentzVector AssignGamma(const Baby &b, bool gen) {
    TLorentzVector gamma;
    bool FoundGamma(false);
    if(gen) {
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == 22 && b.mc_pt()->at(i) > 5) {
          if (!FoundGamma) gamma.SetPtEtaPhiM(b.mc_pt()->at(i),
                                              b.mc_eta()->at(i),
                                              b.mc_phi()->at(i),
                                              b.mc_mass()->at(i));
          FoundGamma = true;
          if(b.mc_mom()->at(i) == 25) {
            //cout<<"Higgs pt:"<<b.mc_pt()->at(b.mc_momidx()->at(i))<<endl;
            //cout<<"AssignGamma H index: "<<b.mc_momidx()->at(i)<<endl;
            //cout<<"AssignGamma H pt: "<<b.mc_pt()->at(b.mc_momidx()->at(i))<<endl;
            gamma.SetPtEtaPhiM(b.mc_pt()->at(i),
                               b.mc_eta()->at(i),
                               b.mc_phi()->at(i),
                               b.mc_mass()->at(i));
            break;
          }
      }
    } else {
      gamma.SetPtEtaPhiM(b.photon_pt() ->at(0),
                         b.photon_eta()->at(0),
                         b.photon_phi()->at(0), 0);
    }
    return gamma;
  }

  TLorentzVector FindGamma(const Baby &b, bool gen) {
    TLorentzVector gamma;

    if(gen){
    			gamma.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0);
    }

    if(!gen){
      for(unsigned int i = 0; i < b.photon_pt()->size(); i++) {
	  	  if(b.photon_pt() -> at(i) > 10 && fabs(b.photon_eta() -> at(i) < 2.6)){ 
         			gamma.SetPtEtaPhiM(b.photon_pt()->at(i), b.photon_eta()->at(i), b.photon_phi()->at(i), 0);
              break;

        }
	    }
      //gamma.SetPtEtaPhiM(0,0,0,0); -- See if this helps seg faults
    }
    return gamma;

  }
     
  double Findlly(const Baby &b) { 
    TLorentzVector photon,ll;
    if(b.nllphoton() > 0){
      return (b.llphoton_m() -> at(0));
    } else {
      photon = FindGamma(b);
      ll = AssignZ(b);
    }
      return (ll + photon).M();
  }



  // Returns Higgs 4-momentum
  TLorentzVector AssignH(const Baby &b, bool gen) {
    TLorentzVector h, y;
    if(gen) {
      bool FoundH(false);
      for(size_t i = 0; i < b.mc_id()->size(); i++) {
        if(b.mc_id()->at(i) == 25) {
          //cout<<"AssignH index:"<<i<<endl;
          //cout<<"AssignH pt:"<<b.mc_pt()->at(i)<<endl;
          FoundH = true;
          h.SetPtEtaPhiM(b.mc_pt()->at(i),
                         b.mc_eta()->at(i),
                         b.mc_phi()->at(i),
                         b.mc_mass()->at(i));
        }
      }
      if(!FoundH) h = AssignZ(b,gen) + AssignGamma(b,gen);
    }
    else {  
	 TLorentzVector z = AssignZ(b);
	 TLorentzVector photon;

   if(b.llphoton_m() -> size() > 0){
  	 //for(unsigned int i = 0; i < b.llphoton_pt()->size(); i++) {
	
        h.SetPtEtaPhiM(b.llphoton_pt()->at(0),
                       b.llphoton_eta()->at(0),
                       b.llphoton_phi()->at(0),
                       b.llphoton_m()->at(0));
      //}
   } else if ( (b.photon_pt() -> size() > 0) && (b.ll_m() -> size() > 0) && (b.llphoton_m() -> size() == 0) ){
        h = AssignZ(b,gen) + AssignGamma(b,gen);

   } else {
        h.SetPtEtaPhiM(0,0,0,0);
   }
  
  }
   return h;  
  }

  //
  // Variables used for defining kinematic angles presented in https://arxiv.org/pdf/1108.2274.pdf
  //

  // Returns 4-momentum of q1 (quark from gluon-gluon fusion)
  //  Defined in Equation 4
  TLorentzVector AssignQ1(const Baby &b, bool gen) {
    TLorentzVector h = AssignH(b,gen);
    TVector3 htran = h.BoostVector();
    htran.SetZ(0);
    h.Boost(-1*htran);
    TLorentzVector k1;
    double pz, E;
    pz = h.Pz() + h.E();
    E  = h.E()  + h.Pz();
    k1.SetPxPyPzE(0,0,pz/2,E/2);
    k1.Boost(htran);
    return k1;
  }


  // Returns 4-momentum of q2 (quark from gluon-gluon fusion)
  //  Defined in Equation 5
  TLorentzVector AssignQ2(const Baby &b, bool gen) {
    TLorentzVector k2;
    TLorentzVector h = AssignH(b,gen);
    TVector3 htran = h.BoostVector();
    htran.SetZ(0);
    h.Boost(-1*htran);
    double pz, E;
    pz = h.Pz() - h.E();
    E  = h.E()  - h.Pz();
    k2.SetPxPyPzE(0,0,pz/2,E/2);
    k2.Boost(htran);
    return k2;
  }
  
  // Returns magnitude of Z candidate 3-momentum 
  //  Defined in Equation 7
  double lambdaZ(const Baby &b, bool gen) {
    TLorentzVector P = AssignH(b, gen);
    TLorentzVector l1 = AssignL1(b,gen);
    TLorentzVector l2 = AssignL2(b,gen);
    TLorentzVector Z = AssignZ(b, gen);
    double M = P.M(), mll = Z.M();
    return sqrt(pow(P.Dot(Z)/M,2)-pow(mll,2));
  }

  // Cosine of angle between lepton 1 and parent Z in Higgs frame 
  //  Defined in Equation 13
  double cos_theta(const Baby &b, bool gen) {
    TLorentzVector P =  AssignH(b,gen);
    TLorentzVector l1 = AssignL1(b,gen);
    TLorentzVector l2 = AssignL2(b,gen);
    double M = P.M();
    double lZ = lambdaZ(b,gen);
    double ctheta = P.Dot(l1-l2)/(M*lZ);
    if(ctheta > 1) ctheta = 0.999;
    if(ctheta <-1) ctheta = -0.999;
    return ctheta;
  }

  // Cosine of angle between incoming quarks and outgoing Zs in higgs frame 
  //  Defined in Equation 8
  double cos_Theta(const Baby &b, bool gen) {
    TLorentzVector H  = AssignH(b,gen);
    TLorentzVector Z  = AssignZ(b,gen);
    TLorentzVector q1 = AssignQ1(b,gen);
    TLorentzVector q2 = AssignQ2(b,gen);
    //cout<<"a H: "<<H.Pt()<<" Z: "<<Z.Pt()<<endl;
    double M = H.M();
    double lZ = lambdaZ(b,gen);
    double cosTheta = Z.Dot(q1-q2)/(M*lZ);
    //cout<<"a q1: "<<q1.Pt()<<" q2: "<<q2.Pt()<<" M: "<<M<<" lZ:"<<lZ<<" cosT: "<<cosTheta<<endl;
    if(abs(cosTheta) > 1.01) cout << "ERROR: cTheta = " << cosTheta <<  endl;
    return cosTheta;
  }

  // Angle of the Z decay plane from the z-axis (defined in Equation 1) in the higgs frame
  //  Defined in Equation 21+22
  double Getphi(const Baby &b, bool gen) {
    TVector3 l1 = AssignL1(b, gen).Vect();
    TVector3 l2 = AssignL2(b, gen).Vect();
    TVector3 q1 = AssignQ1(b, gen).Vect();
    TVector3 Z  = AssignZ( b, gen).Vect();
    double cosphi, sinphi;
    cosphi = -1*l1.Cross(l2).Dot(q1.Cross(Z))/l1.Cross(l2).Mag()/q1.Cross(Z).Mag();
    sinphi = -1*l1.Cross(l2).Dot(q1)/l1.Cross(l2).Mag()/q1.Mag();
    double phi(0);
    if(abs(cosphi) > 1.01) cout << "ERROR: cphi = " << cosphi <<  endl;
    if(cosphi > 1) cosphi = 1;
    if(cosphi < -1) cosphi = -1;
    if(sinphi < 0) phi = -1*acos(cosphi);
    else           phi = acos(cosphi);
    return phi;
  }
}







