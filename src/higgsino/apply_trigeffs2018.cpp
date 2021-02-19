#include <vector>
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/apply_trigeffs2018.hpp"

namespace Higfuncs{

const NamedFunc get_0l_trigeff2018("get_0l_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.0519126; errup = 0.0266522; errdown = 0.0266522;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.102473; errup = 0.0506841; errdown = 0.0506841;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.111702; errup = 0.0336122; errdown = 0.0336122;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.111702; errup = 0.0336122; errdown = 0.0336122;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 180) {eff = 0.177866; errup = 0.0458729; errdown = 0.0458729;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 190) {eff = 0.221311; errup = 0.0614236; errdown = 0.0614236;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 200) {eff = 0.272727; errup = 0.0747522; errdown = 0.0747522;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 9999) {eff = 0.336957; errup = 0.0757904; errdown = 0.0757904;}
  else if (ht> 200 && ht<= 300 && met> 150 && met<= 155) {eff = 0.200972; errup = 0.0942924; errdown = 0.0942924;}
  else if (ht> 200 && ht<= 300 && met> 155 && met<= 160) {eff = 0.259058; errup = 0.121196; errdown = 0.121196;}
  else if (ht> 200 && ht<= 300 && met> 160 && met<= 165) {eff = 0.283133; errup = 0.0653851; errdown = 0.0653851;}
  else if (ht> 200 && ht<= 300 && met> 165 && met<= 170) {eff = 0.331126; errup = 0.076019; errdown = 0.076019;}
  else if (ht> 200 && ht<= 300 && met> 170 && met<= 175) {eff = 0.451128; errup = 0.102174; errdown = 0.102174;}
  else if (ht> 200 && ht<= 300 && met> 175 && met<= 180) {eff = 0.41875; errup = 0.0960246; errdown = 0.0960246;}
  else if (ht> 200 && ht<= 300 && met> 180 && met<= 185) {eff = 0.515504; errup = 0.117364; errdown = 0.117364;}
  else if (ht> 200 && ht<= 300 && met> 185 && met<= 190) {eff = 0.579439; errup = 0.1316; errdown = 0.1316;}
  else if (ht> 200 && ht<= 300 && met> 190 && met<= 195) {eff = 0.568807; errup = 0.129293; errdown = 0.129293;}
  else if (ht> 200 && ht<= 300 && met> 195 && met<= 200) {eff = 0.670968; errup = 0.152051; errdown = 0.152051;}
  else if (ht> 200 && ht<= 300 && met> 200 && met<= 210) {eff = 0.73516; errup = 0.129121; errdown = 0.129121;}
  else if (ht> 200 && ht<= 300 && met> 210 && met<= 220) {eff = 0.774038; errup = 0.135416; errdown = 0.135416;}
  else if (ht> 200 && ht<= 300 && met> 220 && met<= 230) {eff = 0.867089; errup = 0.132911; errdown = 0.150618;}
  else if (ht> 200 && ht<= 300 && met> 230 && met<= 240) {eff = 0.87156; errup = 0.0909439; errdown = 0.0909439;}
  else if (ht> 200 && ht<= 300 && met> 240 && met<= 250) {eff = 0.886076; errup = 0.093621; errdown = 0.093621;}
  else if (ht> 200 && ht<= 300 && met> 250 && met<= 275) {eff = 0.935484; errup = 0.0645161; errdown = 0.0808035;}
  else if (ht> 200 && ht<= 300 && met> 275 && met<= 9999) {eff = 0.972222; errup = 0.0277778; errdown = 0.0842695;}
  else if (ht> 300 && ht<= 400 && met> 150 && met<= 155) {eff = 0.337349; errup = 0.163848; errdown = 0.163848;}
  else if (ht> 300 && ht<= 400 && met> 155 && met<= 160) {eff = 0.360269; errup = 0.175003; errdown = 0.175003;}
  else if (ht> 300 && ht<= 400 && met> 160 && met<= 165) {eff = 0.422925; errup = 0.0939919; errdown = 0.0939919;}
  else if (ht> 300 && ht<= 400 && met> 165 && met<= 170) {eff = 0.466926; errup = 0.102767; errdown = 0.102767;}
  else if (ht> 300 && ht<= 400 && met> 170 && met<= 175) {eff = 0.553398; errup = 0.121137; errdown = 0.121137;}
  else if (ht> 300 && ht<= 400 && met> 175 && met<= 180) {eff = 0.569832; errup = 0.125124; errdown = 0.125124;}
  else if (ht> 300 && ht<= 400 && met> 180 && met<= 185) {eff = 0.701493; errup = 0.134049; errdown = 0.134049;}
  else if (ht> 300 && ht<= 400 && met> 185 && met<= 190) {eff = 0.700535; errup = 0.134175; errdown = 0.134175;}
  else if (ht> 300 && ht<= 400 && met> 190 && met<= 195) {eff = 0.724138; errup = 0.13934; errdown = 0.13934;}
  else if (ht> 300 && ht<= 400 && met> 195 && met<= 200) {eff = 0.826667; errup = 0.156406; errdown = 0.156406;}
  else if (ht> 300 && ht<= 400 && met> 200 && met<= 210) {eff = 0.862069; errup = 0.0880559; errdown = 0.0880559;}
  else if (ht> 300 && ht<= 400 && met> 210 && met<= 220) {eff = 0.867816; errup = 0.0890126; errdown = 0.0890126;}
  else if (ht> 300 && ht<= 400 && met> 220 && met<= 230) {eff = 0.940541; errup = 0.0594595; errdown = 0.0939934;}
  else if (ht> 300 && ht<= 400 && met> 230 && met<= 240) {eff = 0.950413; errup = 0.0401099; errdown = 0.0401099;}
  else if (ht> 300 && ht<= 400 && met> 240 && met<= 250) {eff = 0.982609; errup = 0.0173913; errdown = 0.0381041;}
  else if (ht> 300 && ht<= 400 && met> 250 && met<= 275) {eff = 0.985437; errup = 0.0145631; errdown = 0.0236853;}
  else if (ht> 300 && ht<= 400 && met> 275 && met<= 300) {eff = 0.985507; errup = 0.0144928; errdown = 0.0243905;}
  else if (ht> 300 && ht<= 400 && met> 300 && met<= 9999) {eff = 0.99; errup = 0.01; errdown = 0.0169057;}
  else if (ht> 400 && ht<= 600 && met> 150 && met<= 155) {eff = 0.381271; errup = 0.107837; errdown = 0.107837;}
  else if (ht> 400 && ht<= 600 && met> 155 && met<= 160) {eff = 0.43619; errup = 0.123176; errdown = 0.123176;}
  else if (ht> 400 && ht<= 600 && met> 160 && met<= 165) {eff = 0.446581; errup = 0.0874369; errdown = 0.0874369;}
  else if (ht> 400 && ht<= 600 && met> 165 && met<= 170) {eff = 0.558411; errup = 0.108185; errdown = 0.108185;}
  else if (ht> 400 && ht<= 600 && met> 170 && met<= 175) {eff = 0.581047; errup = 0.112496; errdown = 0.112496;}
  else if (ht> 400 && ht<= 600 && met> 175 && met<= 180) {eff = 0.677778; errup = 0.130386; errdown = 0.130386;}
  else if (ht> 400 && ht<= 600 && met> 180 && met<= 185) {eff = 0.722071; errup = 0.0808773; errdown = 0.0808773;}
  else if (ht> 400 && ht<= 600 && met> 185 && met<= 190) {eff = 0.780899; errup = 0.0865531; errdown = 0.0865531;}
  else if (ht> 400 && ht<= 600 && met> 190 && met<= 195) {eff = 0.780488; errup = 0.0871803; errdown = 0.0871803;}
  else if (ht> 400 && ht<= 600 && met> 195 && met<= 200) {eff = 0.857143; errup = 0.0943144; errdown = 0.0943144;}
  else if (ht> 400 && ht<= 600 && met> 200 && met<= 210) {eff = 0.858427; errup = 0.0715949; errdown = 0.0715949;}
  else if (ht> 400 && ht<= 600 && met> 210 && met<= 220) {eff = 0.925926; errup = 0.0740741; errdown = 0.0762577;}
  else if (ht> 400 && ht<= 600 && met> 220 && met<= 230) {eff = 0.955844; errup = 0.0441558; errdown = 0.0782704;}
  else if (ht> 400 && ht<= 600 && met> 230 && met<= 240) {eff = 0.978723; errup = 0.0212766; errdown = 0.0309743;}
  else if (ht> 400 && ht<= 600 && met> 240 && met<= 250) {eff = 0.971631; errup = 0.0283688; errdown = 0.0313196;}
  else if (ht> 400 && ht<= 600 && met> 250 && met<= 275) {eff = 0.994444; errup = 0.00555556; errdown = 0.022177;}
  else if (ht> 400 && ht<= 600 && met> 275 && met<= 300) {eff = 0.987981; errup = 0.0120192; errdown = 0.0224476;}
  else if (ht> 400 && ht<= 600 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00358039;}
  else if (ht> 600 && ht<= 950 && met> 150 && met<= 155) {eff = 0.400491; errup = 0.0744788; errdown = 0.0744788;}
  else if (ht> 600 && ht<= 950 && met> 155 && met<= 160) {eff = 0.504021; errup = 0.0923124; errdown = 0.0923124;}
  else if (ht> 600 && ht<= 950 && met> 160 && met<= 165) {eff = 0.507003; errup = 0.103808; errdown = 0.103808;}
  else if (ht> 600 && ht<= 950 && met> 165 && met<= 170) {eff = 0.551155; errup = 0.1128; errdown = 0.1128;}
  else if (ht> 600 && ht<= 950 && met> 170 && met<= 175) {eff = 0.640927; errup = 0.130348; errdown = 0.130348;}
  else if (ht> 600 && ht<= 950 && met> 175 && met<= 180) {eff = 0.743191; errup = 0.149643; errdown = 0.149643;}
  else if (ht> 600 && ht<= 950 && met> 180 && met<= 185) {eff = 0.725806; errup = 0.0445408; errdown = 0.0445408;}
  else if (ht> 600 && ht<= 950 && met> 185 && met<= 190) {eff = 0.835498; errup = 0.046481; errdown = 0.046481;}
  else if (ht> 600 && ht<= 950 && met> 190 && met<= 195) {eff = 0.797297; errup = 0.046407; errdown = 0.046407;}
  else if (ht> 600 && ht<= 950 && met> 195 && met<= 200) {eff = 0.850299; errup = 0.0488229; errdown = 0.0488229;}
  else if (ht> 600 && ht<= 950 && met> 200 && met<= 210) {eff = 0.886628; errup = 0.0654922; errdown = 0.0654922;}
  else if (ht> 600 && ht<= 950 && met> 210 && met<= 220) {eff = 0.926357; errup = 0.068027; errdown = 0.068027;}
  else if (ht> 600 && ht<= 950 && met> 220 && met<= 230) {eff = 0.933333; errup = 0.0666667; errdown = 0.0684725;}
  else if (ht> 600 && ht<= 950 && met> 230 && met<= 240) {eff = 0.982759; errup = 0.0172414; errdown = 0.0310059;}
  else if (ht> 600 && ht<= 950 && met> 240 && met<= 250) {eff = 0.980583; errup = 0.0194175; errdown = 0.0312543;}
  else if (ht> 600 && ht<= 950 && met> 250 && met<= 275) {eff = 0.989333; errup = 0.0106667; errdown = 0.0226867;}
  else if (ht> 600 && ht<= 950 && met> 275 && met<= 300) {eff = 0.996753; errup = 0.00324675; errdown = 0.0224584;}
  else if (ht> 600 && ht<= 950 && met> 300 && met<= 9999) {eff = 0.996324; errup = 0.00367647; errdown = 0.00711638;}
  else if (ht> 950 && ht<= 9999 && met> 150 && met<= 160) {eff = 0.493151; errup = 0.110839; errdown = 0.110839;}
  else if (ht> 950 && ht<= 9999 && met> 160 && met<= 170) {eff = 0.644444; errup = 0.103253; errdown = 0.103253;}
  else if (ht> 950 && ht<= 9999 && met> 170 && met<= 180) {eff = 0.681416; errup = 0.109285; errdown = 0.109285;}
  else if (ht> 950 && ht<= 9999 && met> 180 && met<= 190) {eff = 0.863158; errup = 0.100786; errdown = 0.100786;}
  else if (ht> 950 && ht<= 9999 && met> 190 && met<= 200) {eff = 0.847826; errup = 0.100015; errdown = 0.100015;}
  else if (ht> 950 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.903614; errup = 0.0574661; errdown = 0.0574661;}
  else if (ht> 950 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.904762; errup = 0.0602203; errdown = 0.0602203;}
  else if (ht> 950 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.953125; errup = 0.046875; errdown = 0.0566106;}
  else if (ht> 950 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.979167; errup = 0.0208333; errdown = 0.0411332;}
  else if (ht> 950 && ht<= 9999 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.0363517;}
  else if (ht> 950 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.965909; errup = 0.0291201; errdown = 0.0291201;}
  else if (ht> 950 && ht<= 9999 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.0225349;}
  else if (ht> 950 && ht<= 9999 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.021753;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_fakemet_trigeff2018("get_0l_fakemet_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 350 && met> 150 && met<= 160) {eff = 0.078; errup = 0.011993; errdown = 0.011993;}
  else if (ht> 0 && ht<= 350 && met> 160 && met<= 170) {eff = 0.102837; errup = 0.0180878; errdown = 0.0180878;}
  else if (ht> 0 && ht<= 350 && met> 170 && met<= 180) {eff = 0.186207; errup = 0.0323274; errdown = 0.0323274;}
  else if (ht> 0 && ht<= 350 && met> 180 && met<= 190) {eff = 0.253731; errup = 0.0531615; errdown = 0.0531615;}
  else if (ht> 0 && ht<= 350 && met> 190 && met<= 200) {eff = 0.363636; errup = 0.0725204; errdown = 0.0725204;}
  else if (ht> 0 && ht<= 350 && met> 200 && met<= 225) {eff = 0.428571; errup = 0.0763604; errdown = 0.0763604;}
  else if (ht> 0 && ht<= 350 && met> 225 && met<= 250) {eff = 0.529412; errup = 0.121058; errdown = 0.121058;}
  else if (ht> 0 && ht<= 350 && met> 250 && met<= 9999) {eff = 0.823529; errup = 0.0924594; errdown = 0.0924594;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.153846; errup = 0.0114786; errdown = 0.0114786;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.164521; errup = 0.014023; errdown = 0.014023;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.227273; errup = 0.0178692; errdown = 0.0178692;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.239374; errup = 0.0201823; errdown = 0.0201823;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.285294; errup = 0.024489; errdown = 0.024489;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.309353; errup = 0.0277225; errdown = 0.0277225;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.402655; errup = 0.0326231; errdown = 0.0326231;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.411111; errup = 0.0366741; errdown = 0.0366741;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.503876; errup = 0.0440212; errdown = 0.0440212;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.473684; errup = 0.0512278; errdown = 0.0512278;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.649351; errup = 0.0384517; errdown = 0.0384517;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.666667; errup = 0.0488824; errdown = 0.0488824;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.860759; errup = 0.0389502; errdown = 0.0389502;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.722222; errup = 0.0609519; errdown = 0.0609519;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.837209; errup = 0.0562986; errdown = 0.0562986;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 300) {eff = 0.93617; errup = 0.025213; errdown = 0.025213;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.973684; errup = 0.0259672; errdown = 0.0259672;}
  else if (ht> 450 && ht<= 550 && met> 150 && met<= 155) {eff = 0.185145; errup = 0.0090772; errdown = 0.0090772;}
  else if (ht> 450 && ht<= 550 && met> 155 && met<= 160) {eff = 0.224373; errup = 0.0111693; errdown = 0.0111693;}
  else if (ht> 450 && ht<= 550 && met> 160 && met<= 165) {eff = 0.253235; errup = 0.0132203; errdown = 0.0132203;}
  else if (ht> 450 && ht<= 550 && met> 165 && met<= 170) {eff = 0.296253; errup = 0.0156247; errdown = 0.0156247;}
  else if (ht> 450 && ht<= 550 && met> 170 && met<= 175) {eff = 0.351145; errup = 0.0186508; errdown = 0.0186508;}
  else if (ht> 450 && ht<= 550 && met> 175 && met<= 180) {eff = 0.397363; errup = 0.0212361; errdown = 0.0212361;}
  else if (ht> 450 && ht<= 550 && met> 180 && met<= 185) {eff = 0.430233; errup = 0.0238763; errdown = 0.0238763;}
  else if (ht> 450 && ht<= 550 && met> 185 && met<= 190) {eff = 0.512894; errup = 0.0267555; errdown = 0.0267555;}
  else if (ht> 450 && ht<= 550 && met> 190 && met<= 195) {eff = 0.553191; errup = 0.0296056; errdown = 0.0296056;}
  else if (ht> 450 && ht<= 550 && met> 195 && met<= 200) {eff = 0.571429; errup = 0.033065; errdown = 0.033065;}
  else if (ht> 450 && ht<= 550 && met> 200 && met<= 210) {eff = 0.626263; errup = 0.0243116; errdown = 0.0243116;}
  else if (ht> 450 && ht<= 550 && met> 210 && met<= 220) {eff = 0.731707; errup = 0.0282492; errdown = 0.0282492;}
  else if (ht> 450 && ht<= 550 && met> 220 && met<= 230) {eff = 0.797753; errup = 0.0301069; errdown = 0.0301069;}
  else if (ht> 450 && ht<= 550 && met> 230 && met<= 240) {eff = 0.758865; errup = 0.0360249; errdown = 0.0360249;}
  else if (ht> 450 && ht<= 550 && met> 240 && met<= 250) {eff = 0.824074; errup = 0.0366384; errdown = 0.0366384;}
  else if (ht> 450 && ht<= 550 && met> 250 && met<= 300) {eff = 0.90458; errup = 0.0181507; errdown = 0.0181507;}
  else if (ht> 450 && ht<= 550 && met> 300 && met<= 400) {eff = 0.948387; errup = 0.0177708; errdown = 0.0177708;}
  else if (ht> 450 && ht<= 550 && met> 400 && met<= 9999) {eff = 0.959459; errup = 0.0229267; errdown = 0.0229267;}
  else if (ht> 550 && ht<= 650 && met> 150 && met<= 155) {eff = 0.2286; errup = 0.00942298; errdown = 0.00942298;}
  else if (ht> 550 && ht<= 650 && met> 155 && met<= 160) {eff = 0.262344; errup = 0.0107295; errdown = 0.0107295;}
  else if (ht> 550 && ht<= 650 && met> 160 && met<= 165) {eff = 0.332286; errup = 0.0132019; errdown = 0.0132019;}
  else if (ht> 550 && ht<= 650 && met> 165 && met<= 170) {eff = 0.356405; errup = 0.0153936; errdown = 0.0153936;}
  else if (ht> 550 && ht<= 650 && met> 170 && met<= 175) {eff = 0.418685; errup = 0.0167548; errdown = 0.0167548;}
  else if (ht> 550 && ht<= 650 && met> 175 && met<= 180) {eff = 0.471642; errup = 0.0192856; errdown = 0.0192856;}
  else if (ht> 550 && ht<= 650 && met> 180 && met<= 185) {eff = 0.472222; errup = 0.0222374; errdown = 0.0222374;}
  else if (ht> 550 && ht<= 650 && met> 185 && met<= 190) {eff = 0.543033; errup = 0.02255; errdown = 0.02255;}
  else if (ht> 550 && ht<= 650 && met> 190 && met<= 195) {eff = 0.584699; errup = 0.0257577; errdown = 0.0257577;}
  else if (ht> 550 && ht<= 650 && met> 195 && met<= 200) {eff = 0.635328; errup = 0.0256919; errdown = 0.0256919;}
  else if (ht> 550 && ht<= 650 && met> 200 && met<= 210) {eff = 0.669922; errup = 0.0207819; errdown = 0.0207819;}
  else if (ht> 550 && ht<= 650 && met> 210 && met<= 220) {eff = 0.737113; errup = 0.0223478; errdown = 0.0223478;}
  else if (ht> 550 && ht<= 650 && met> 220 && met<= 230) {eff = 0.744409; errup = 0.0246551; errdown = 0.0246551;}
  else if (ht> 550 && ht<= 650 && met> 230 && met<= 240) {eff = 0.791304; errup = 0.0267957; errdown = 0.0267957;}
  else if (ht> 550 && ht<= 650 && met> 240 && met<= 250) {eff = 0.828125; errup = 0.0272272; errdown = 0.0272272;}
  else if (ht> 550 && ht<= 650 && met> 250 && met<= 275) {eff = 0.881459; errup = 0.0178212; errdown = 0.0178212;}
  else if (ht> 550 && ht<= 650 && met> 275 && met<= 300) {eff = 0.88587; errup = 0.023441; errdown = 0.023441;}
  else if (ht> 550 && ht<= 650 && met> 300 && met<= 400) {eff = 0.887029; errup = 0.014479; errdown = 0.014479;}
  else if (ht> 550 && ht<= 650 && met> 400 && met<= 9999) {eff = 0.977528; errup = 0.00473689; errdown = 0.00473689;}
  else if (ht> 650 && ht<= 800 && met> 150 && met<= 155) {eff = 0.312199; errup = 0.00737961; errdown = 0.00737961;}
  else if (ht> 650 && ht<= 800 && met> 155 && met<= 160) {eff = 0.355075; errup = 0.00835433; errdown = 0.00835433;}
  else if (ht> 650 && ht<= 800 && met> 160 && met<= 165) {eff = 0.397912; errup = 0.00962519; errdown = 0.00962519;}
  else if (ht> 650 && ht<= 800 && met> 165 && met<= 170) {eff = 0.446444; errup = 0.0104157; errdown = 0.0104157;}
  else if (ht> 650 && ht<= 800 && met> 170 && met<= 175) {eff = 0.505056; errup = 0.0115341; errdown = 0.0115341;}
  else if (ht> 650 && ht<= 800 && met> 175 && met<= 180) {eff = 0.540892; errup = 0.012404; errdown = 0.012404;}
  else if (ht> 650 && ht<= 800 && met> 180 && met<= 185) {eff = 0.593127; errup = 0.0128787; errdown = 0.0128787;}
  else if (ht> 650 && ht<= 800 && met> 185 && met<= 190) {eff = 0.609113; errup = 0.0137958; errdown = 0.0137958;}
  else if (ht> 650 && ht<= 800 && met> 190 && met<= 195) {eff = 0.628141; errup = 0.0153217; errdown = 0.0153217;}
  else if (ht> 650 && ht<= 800 && met> 195 && met<= 200) {eff = 0.692389; errup = 0.0150048; errdown = 0.0150048;}
  else if (ht> 650 && ht<= 800 && met> 200 && met<= 210) {eff = 0.703607; errup = 0.011694; errdown = 0.011694;}
  else if (ht> 650 && ht<= 800 && met> 210 && met<= 220) {eff = 0.790055; errup = 0.0114418; errdown = 0.0114418;}
  else if (ht> 650 && ht<= 800 && met> 220 && met<= 230) {eff = 0.805634; errup = 0.0121256; errdown = 0.0121256;}
  else if (ht> 650 && ht<= 800 && met> 230 && met<= 240) {eff = 0.862187; errup = 0.0116332; errdown = 0.0116332;}
  else if (ht> 650 && ht<= 800 && met> 240 && met<= 250) {eff = 0.868639; errup = 0.0116205; errdown = 0.0116205;}
  else if (ht> 650 && ht<= 800 && met> 250 && met<= 275) {eff = 0.878536; errup = 0.00769318; errdown = 0.00769318;}
  else if (ht> 650 && ht<= 800 && met> 275 && met<= 300) {eff = 0.922261; errup = 0.00711815; errdown = 0.00711815;}
  else if (ht> 650 && ht<= 800 && met> 300 && met<= 350) {eff = 0.934977; errup = 0.00508305; errdown = 0.00508305;}
  else if (ht> 650 && ht<= 800 && met> 350 && met<= 400) {eff = 0.95888; errup = 0.00484746; errdown = 0.00484746;}
  else if (ht> 650 && ht<= 800 && met> 400 && met<= 450) {eff = 0.969673; errup = 0.00484455; errdown = 0.00484455;}
  else if (ht> 650 && ht<= 800 && met> 450 && met<= 500) {eff = 0.979543; errup = 0.00491061; errdown = 0.00491061;}
  else if (ht> 650 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.986716; errup = 0.00311023; errdown = 0.00311023;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.539823; errup = 0.0034145; errdown = 0.0034145;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.571429; errup = 0.00358377; errdown = 0.00358377;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.61682; errup = 0.00383997; errdown = 0.00383997;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.65604; errup = 0.00399833; errdown = 0.00399833;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.675896; errup = 0.00424876; errdown = 0.00424876;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.719656; errup = 0.00434348; errdown = 0.00434348;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.735982; errup = 0.00455557; errdown = 0.00455557;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.765251; errup = 0.00464919; errdown = 0.00464919;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.773482; errup = 0.00482465; errdown = 0.00482465;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.806076; errup = 0.00486078; errdown = 0.00486078;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.833808; errup = 0.00351213; errdown = 0.00351213;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.858087; errup = 0.00366153; errdown = 0.00366153;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.874377; errup = 0.00379521; errdown = 0.00379521;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.891151; errup = 0.00394335; errdown = 0.00394335;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.909908; errup = 0.0039598; errdown = 0.0039598;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.920343; errup = 0.00270517; errdown = 0.00270517;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.940531; errup = 0.00282431; errdown = 0.00282431;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.943467; errup = 0.0025356; errdown = 0.0025356;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.952579; errup = 0.00332293; errdown = 0.00332293;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.968186; errup = 0.00382436; errdown = 0.00382436;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.969146; errup = 0.00486385; errdown = 0.00486385;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.987481; errup = 0.00248803; errdown = 0.00248803;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.432701; errup = 0.00130445; errdown = 0.00130445;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.468416; errup = 0.00141951; errdown = 0.00141951;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.500201; errup = 0.00153028; errdown = 0.00153028;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.534339; errup = 0.00164201; errdown = 0.00164201;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.565809; errup = 0.00176191; errdown = 0.00176191;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.597802; errup = 0.00187686; errdown = 0.00187686;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.624353; errup = 0.00198579; errdown = 0.00198579;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.657567; errup = 0.00209913; errdown = 0.00209913;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.683468; errup = 0.00220487; errdown = 0.00220487;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.707921; errup = 0.00231018; errdown = 0.00231018;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.736016; errup = 0.00174943; errdown = 0.00174943;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.77612; errup = 0.0018808; errdown = 0.0018808;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.81054; errup = 0.00199175; errdown = 0.00199175;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.831745; errup = 0.00215975; errdown = 0.00215975;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.855244; errup = 0.00226767; errdown = 0.00226767;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.882622; errup = 0.00156784; errdown = 0.00156784;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.90554; errup = 0.00178883; errdown = 0.00178883;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.927802; errup = 0.00151072; errdown = 0.00151072;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.941449; errup = 0.00195205; errdown = 0.00195205;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.944989; errup = 0.00260629; errdown = 0.00260629;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.954409; errup = 0.00317334; errdown = 0.00317334;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.960595; errup = 0.00225241; errdown = 0.00225241;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1el_trigeff2018("get_1el_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.0097629; errup = 0.00521891; errdown = 0.00359356;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.238559; errup = 0.0140446; errdown = 0.0135171;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.523018; errup = 0.00945585; errdown = 0.00947208;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.017301; errup = 0.00729105; errdown = 0.00535648;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.297573; errup = 0.0172358; errdown = 0.0166989;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.583857; errup = 0.0110342; errdown = 0.011117;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.260109; errup = 0.0153141; errdown = 0.0147704;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.444776; errup = 0.0162215; errdown = 0.0161081;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.650049; errup = 0.00869262; errdown = 0.00879159;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.387642; errup = 0.0066101; errdown = 0.00656955;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.596085; errup = 0.00792443; errdown = 0.00797407;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.751802; errup = 0.00372287; errdown = 0.00376028;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.386105; errup = 0.0171245; errdown = 0.0168562;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.643836; errup = 0.0259838; errdown = 0.0268136;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.754565; errup = 0.0118009; errdown = 0.0121843;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.761468; errup = 0.00860319; errdown = 0.00881662;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.800239; errup = 0.0140459; errdown = 0.0147947;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.87316; errup = 0.00563279; errdown = 0.0058494;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.00446429; errup = 0.0101903; errdown = 0.00369336;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.267442; errup = 0.0230545; errdown = 0.0219156;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.551826; errup = 0.0123281; errdown = 0.0123906;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.0408163; errup = 0.0195138; errdown = 0.0140021;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.307937; errup = 0.0282018; errdown = 0.0269105;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.622134; errup = 0.0139227; errdown = 0.0141213;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.33518; errup = 0.0266756; errdown = 0.0257133;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.524887; errup = 0.024796; errdown = 0.0249138;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.689604; errup = 0.010443; errdown = 0.0106352;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.490962; errup = 0.0112954; errdown = 0.0112864;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.663012; errup = 0.0115434; errdown = 0.0117363;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.785698; errup = 0.0044274; errdown = 0.00449406;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.470199; errup = 0.0304274; errdown = 0.0302191;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.636364; errup = 0.0394636; errdown = 0.0412433;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.793506; errup = 0.0148321; errdown = 0.0156277;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.790724; errup = 0.0138995; errdown = 0.014584;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.872; errup = 0.0175109; errdown = 0.0196448;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.897516; errup = 0.00642761; errdown = 0.00679249;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.0361446; errup = 0.0339112; errdown = 0.0196074;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.278607; errup = 0.0352149; errdown = 0.0328462;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.65252; errup = 0.0144558; errdown = 0.0147338;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.045977; errup = 0.0348613; errdown = 0.0218756;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.409639; errup = 0.0416348; errdown = 0.0404647;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.704675; errup = 0.0157244; errdown = 0.0162081;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.453988; errup = 0.0422296; errdown = 0.041623;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.61597; errup = 0.0313617; errdown = 0.0322973;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.763399; errup = 0.0110119; errdown = 0.0113664;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.626178; errup = 0.0160347; errdown = 0.0163076;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.76572; errup = 0.0137106; errdown = 0.0142695;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.825713; errup = 0.00481489; errdown = 0.00492038;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.610687; errup = 0.0453152; errdown = 0.0471475;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.775; errup = 0.0487807; errdown = 0.0565118;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.859083; errup = 0.0145352; errdown = 0.0158254;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.850498; errup = 0.0209401; errdown = 0.0234521;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.859459; errup = 0.0260778; errdown = 0.0303603;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.914036; errup = 0.00731012; errdown = 0.00789175;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.111111; errup = 0.0791717; errdown = 0.0524058;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.3125; errup = 0.0543931; errdown = 0.0500862;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.710778; errup = 0.0171191; errdown = 0.0177172;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0; errup = 0.0683597; errdown = 0;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.363636; errup = 0.0583279; errdown = 0.0549108;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.755848; errup = 0.0167559; errdown = 0.0175366;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.588235; errup = 0.0577422; errdown = 0.0600294;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.702479; errup = 0.0436783; errdown = 0.0473382;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.807843; errup = 0.012509; errdown = 0.0131362;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.705495; errup = 0.0219711; errdown = 0.0229195;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.808732; errup = 0.018281; errdown = 0.0196374;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.868353; errup = 0.0052035; errdown = 0.00537996;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.865385; errup = 0.0482807; errdown = 0.0649611;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.942308; errup = 0.0312343; errdown = 0.0529489;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.876147; errup = 0.0159935; errdown = 0.0178418;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.898438; errup = 0.0270367; errdown = 0.0341365;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.9; errup = 0.0338535; errdown = 0.0456142;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.933071; errup = 0.00788272; errdown = 0.00878826;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.1875; errup = 0.149399; errdown = 0.100212;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.583333; errup = 0.0918854; errdown = 0.0971962;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.799163; errup = 0.0186952; errdown = 0.0200167;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.238095; errup = 0.128988; errdown = 0.0987131;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.673077; errup = 0.0701338; errdown = 0.0776765;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.842784; errup = 0.0188114; errdown = 0.0206985;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.522727; errup = 0.0849595; errdown = 0.0861304;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.786885; errup = 0.0548156; errdown = 0.0655726;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.887324; errup = 0.0119924; errdown = 0.0131472;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.865546; errup = 0.0225156; errdown = 0.0258705;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.885621; errup = 0.0184512; errdown = 0.0211925;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.908734; errup = 0.00523456; errdown = 0.00551032;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.857143; errup = 0.0670785; errdown = 0.0986253;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.9; errup = 0.0643201; errdown = 0.116971;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.910448; errup = 0.0176206; errdown = 0.0210035;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.915254; errup = 0.0361154; errdown = 0.0532642;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.940299; errup = 0.0283546; errdown = 0.0446914;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.953344; errup = 0.00833639; errdown = 0.00987704;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.5; errup = 0.195182; errdown = 0.195182;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.666667; errup = 0.106107; errdown = 0.122517;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.871875; errup = 0.0189792; errdown = 0.021493;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.666667; errup = 0.17521; errdown = 0.221361;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.772727; errup = 0.0944237; errdown = 0.12452;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.9375; errup = 0.014321; errdown = 0.0177218;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.764706; errup = 0.108959; errdown = 0.147312;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.913043; errup = 0.0411491; errdown = 0.0634308;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.924335; errup = 0.0120474; errdown = 0.013929;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.894309; errup = 0.028102; errdown = 0.0354144;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.94; errup = 0.0168114; errdown = 0.0218501;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.952295; errup = 0.00453493; errdown = 0.00496221;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.142229;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.15411;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.969543; errup = 0.0120157; errdown = 0.0177485;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.969697; errup = 0.0250817; errdown = 0.0662602;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.972973; errup = 0.0223689; errdown = 0.0594217;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.970297; errup = 0.00752031; errdown = 0.00962462;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.666667; errup = 0.277375; errdown = 0.414535;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.8; errup = 0.0835235; errdown = 0.112668;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.952381; errup = 0.0116187; errdown = 0.0146509;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.777778; errup = 0.142118; errdown = 0.221429;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.8; errup = 0.0835235; errdown = 0.112668;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.983444; errup = 0.00713312; errdown = 0.0110449;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.9; errup = 0.0643201; errdown = 0.116971;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.918919; errup = 0.0438002; errdown = 0.0726265;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.965649; errup = 0.00793934; errdown = 0.00992767;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.957983; errup = 0.0180304; errdown = 0.027424;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.960199; errup = 0.0136568; errdown = 0.0190433;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.984898; errup = 0.00245839; errdown = 0.00288698;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.168149;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.95; errup = 0.0413995; errdown = 0.105764;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.975845; errup = 0.0103945; errdown = 0.0160098;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0576587;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0472931;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.994718; errup = 0.00287314; errdown = 0.0051109;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.15411;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0879414;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.991758; errup = 0.00448218; errdown = 0.00795189;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.264229;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 0.947368; errup = 0.0435805; errdown = 0.110836;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 0.988201; errup = 0.00563868; errdown = 0.00923117;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.184992;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 0.979592; errup = 0.0168888; errdown = 0.0453679;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 0.995208; errup = 0.00260705; errdown = 0.00463962;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0200277;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0104058;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.996068; errup = 0.00111791; errdown = 0.00148999;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.308024;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0802771;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 0.993266; errup = 0.00434837; errdown = 0.00881245;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.102638;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0271039;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.99867; errup = 0.00110009; errdown = 0.00305118;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1mu_trigeff2018("get_1mu_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.0509554; errup = 0.0103616; errdown = 0.00880306;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.406528; errup = 0.0160299; errdown = 0.0158392;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.920266; errup = 0.00495963; errdown = 0.00524822;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.296203; errup = 0.024754; errdown = 0.0236669;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.57037; errup = 0.019671; errdown = 0.0198874;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.918153; errup = 0.00631369; errdown = 0.00677068;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.45394; errup = 0.0171845; errdown = 0.017079;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.734914; errup = 0.0120119; errdown = 0.0123595;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.94219; errup = 0.00339353; errdown = 0.0035852;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.904807; errup = 0.00908061; errdown = 0.00988317;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.939845; errup = 0.00607334; errdown = 0.00667259;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.953995; errup = 0.00248988; errdown = 0.00262123;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.0708333; errup = 0.0206706; errdown = 0.0166547;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.462222; errup = 0.0246733; errdown = 0.0244977;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.905621; errup = 0.00668024; errdown = 0.00711475;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.358382; errup = 0.0401117; errdown = 0.0383512;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.66548; errup = 0.0292488; errdown = 0.0304982;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.924117; errup = 0.00730077; errdown = 0.00797078;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.5629; errup = 0.0238039; errdown = 0.0240843;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.770013; errup = 0.015271; errdown = 0.0159857;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.950936; errup = 0.0036475; errdown = 0.00391321;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.923333; errup = 0.010941; errdown = 0.0124589;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.946875; errup = 0.00726638; errdown = 0.00826383;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.95772; errup = 0.00278904; errdown = 0.00297018;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.105263; errup = 0.0415156; errdown = 0.0318852;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.44898; errup = 0.0339661; errdown = 0.033524;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.924188; errup = 0.0071534; errdown = 0.0077968;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.410256; errup = 0.0629437; errdown = 0.0603993;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.71978; errup = 0.0346421; errdown = 0.0372475;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.93621; errup = 0.00752273; errdown = 0.00839199;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.681159; errup = 0.0291135; errdown = 0.0305079;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.833676; errup = 0.0171643; errdown = 0.0186162;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.953016; errup = 0.00387315; errdown = 0.00418804;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.918539; errup = 0.0146216; errdown = 0.0171965;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.956954; errup = 0.00826908; errdown = 0.00993237;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.965588; errup = 0.00275628; errdown = 0.00297726;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.111111; errup = 0.0791717; errdown = 0.0524058;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.519685; errup = 0.0478966; errdown = 0.0482324;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.938907; errup = 0.00787908; errdown = 0.00888368;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.68; errup = 0.102272; errdown = 0.119276;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.788462; errup = 0.0415895; errdown = 0.0477639;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.947566; errup = 0.00790478; errdown = 0.00911009;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.703125; errup = 0.042389; errdown = 0.0458523;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.882759; errup = 0.0191698; errdown = 0.0220444;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.965564; errup = 0.00386237; errdown = 0.00430287;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.976303; errup = 0.0101982; errdown = 0.0157124;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.971354; errup = 0.00844407; errdown = 0.01127;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.974014; errup = 0.00275157; errdown = 0.00305009;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.1; errup = 0.116971; errdown = 0.0643201;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.627119; errup = 0.0683561; errdown = 0.0731904;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.973875; errup = 0.00605301; errdown = 0.00759168;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.75; errup = 0.115499; errdown = 0.153966;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.869565; errup = 0.0504938; errdown = 0.0697763;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.97757; errup = 0.00634716; errdown = 0.00839086;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.923077; errup = 0.032824; errdown = 0.0486878;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.923858; errup = 0.01902; errdown = 0.0238853;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.980064; errup = 0.00353709; errdown = 0.0042104;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.963964; errup = 0.0171678; errdown = 0.0275761;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.995305; errup = 0.00388411; errdown = 0.0107125;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.984236; errup = 0.00249933; errdown = 0.00292235;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.875; errup = 0.103637; errdown = 0.23225;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.723404; errup = 0.0696186; errdown = 0.0805167;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.987923; errup = 0.00520699; errdown = 0.00808753;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.941176; errup = 0.048713; errdown = 0.12258;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.962963; errup = 0.0306592; errdown = 0.080075;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.980978; errup = 0.00698888; errdown = 0.0100953;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.971429; errup = 0.0236478; errdown = 0.0626552;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.957447; errup = 0.0202555; errdown = 0.0323679;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.981584; errup = 0.00406048; errdown = 0.00504708;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.967742; errup = 0.0208084; errdown = 0.0409676;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.986667; errup = 0.00860748; errdown = 0.0173148;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.992337; errup = 0.00201909; errdown = 0.00263025;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.9; errup = 0.082873; errdown = 0.194135;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.9; errup = 0.0539222; errdown = 0.0877974;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.989035; errup = 0.00472821; errdown = 0.00734954;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.857143; errup = 0.11848; errdown = 0.257124;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0659133;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.994819; errup = 0.00334597; errdown = 0.00679283;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.971429; errup = 0.0236478; errdown = 0.0626552;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.98913; errup = 0.00899357; errdown = 0.0245495;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.989756; errup = 0.00279678; errdown = 0.00367595;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.967213; errup = 0.0211491; errdown = 0.0416131;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0133482;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.99708; errup = 0.00115753; errdown = 0.0017398;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0769247;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00350724;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0802771;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 0.995204; errup = 0.00309728; errdown = 0.00629067;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0636358;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.016449;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 0.998025; errup = 0.00107467; errdown = 0.00191738;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0439098;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.01428;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.998493; errup = 0.000720896; errdown = 0.00118964;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2el_trigeff2018("get_2el_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.897959; errup = 0.0249664; errdown = 0.0249664;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.950413; errup = 0.0139551; errdown = 0.0139551;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.968085; errup = 0.0104672; errdown = 0.0104672;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.95203; errup = 0.0129816; errdown = 0.0129816;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.92926; errup = 0.0145385; errdown = 0.0145385;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.953846; errup = 0.0130124; errdown = 0.0130124;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.962963; errup = 0.0114932; errdown = 0.0114932;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.960937; errup = 0.012109; errdown = 0.012109;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.952381; errup = 0.0140117; errdown = 0.0140117;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.961353; errup = 0.0133973; errdown = 0.0133973;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.975248; errup = 0.0109318; errdown = 0.0109318;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.978022; errup = 0.0108676; errdown = 0.0108676;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.979899; errup = 0.00994873; errdown = 0.00994873;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.987578; errup = 0.00872921; errdown = 0.00872921;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.996761; errup = 0.000784379; errdown = 0.000784379;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2mu_trigeff2018("get_2mu_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.984496; errup = 0.00769161; errdown = 0.00769161;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.960106; errup = 0.0100929; errdown = 0.0100929;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.988916; errup = 0.000322656; errdown = 0.000322656;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

}
