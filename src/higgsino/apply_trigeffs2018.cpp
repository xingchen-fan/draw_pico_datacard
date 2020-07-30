#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/apply_trigeffs2018.hpp"

namespace Higfuncs{

const NamedFunc get_0l_trigeff2018("get_0l_trigeff2018", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.113095; errup = 0.00624072; errdown = 0.00596183;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.188875; errup = 0.00345936; errdown = 0.00341194;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.294183; errup = 0.0160298; errdown = 0.0155528;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.298851; errup = 0.031016; errdown = 0.0293744;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.376238; errup = 0.0542182; errdown = 0.0515331;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.135683; errup = 0.00757278; errdown = 0.00724204;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.228229; errup = 0.00396375; errdown = 0.00391636;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.365462; errup = 0.0184634; errdown = 0.0180898;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.381295; errup = 0.0313246; errdown = 0.03042;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.44898; errup = 0.0556081; errdown = 0.0544687;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.151927; errup = 0.00903262; errdown = 0.00862502;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.277572; errup = 0.00449048; errdown = 0.00444667;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.403631; errup = 0.01916; errdown = 0.0188808;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.423077; errup = 0.0328108; errdown = 0.0321835;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.4375; errup = 0.0621526; errdown = 0.0604283;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.19427; errup = 0.0110263; errdown = 0.010585;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.327308; errup = 0.00496953; errdown = 0.00493151;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.474398; errup = 0.0201556; errdown = 0.0200756;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.55336; errup = 0.0329386; errdown = 0.033386;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.492958; errup = 0.0659664; errdown = 0.0657465;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.211255; errup = 0.012699; errdown = 0.0121812;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.38853; errup = 0.00545869; errdown = 0.00543118;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.481229; errup = 0.0215077; errdown = 0.0214411;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.54717; errup = 0.0362169; errdown = 0.0366912;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.492537; errup = 0.0680949; errdown = 0.0678474;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.283391; errup = 0.0161913; errdown = 0.0156695;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.425364; errup = 0.00595615; errdown = 0.00593482;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.560892; errup = 0.0212838; errdown = 0.0215011;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.610526; errup = 0.0372811; errdown = 0.0385269;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.577778; errup = 0.0563216; errdown = 0.0582206;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.313754; errup = 0.0185449; errdown = 0.0179899;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.48824; errup = 0.00636054; errdown = 0.00635678;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.62549; errup = 0.0221331; errdown = 0.0226471;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.670213; errup = 0.0358777; errdown = 0.0378223;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.666667; errup = 0.0654592; errdown = 0.0717057;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.369388; errup = 0.0230792; errdown = 0.0225219;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.513148; errup = 0.00679507; errdown = 0.00679987;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.689362; errup = 0.0219642; errdown = 0.0228091;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.689266; errup = 0.0363549; errdown = 0.0386589;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.690141; errup = 0.0584835; errdown = 0.0644708;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.405896; errup = 0.0247128; errdown = 0.0242659;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.591733; errup = 0.00709589; errdown = 0.00713379;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.746237; errup = 0.0206791; errdown = 0.0217931;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.723077; errup = 0.0410611; errdown = 0.0448107;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.779661; errup = 0.0565389; errdown = 0.0673831;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.407104; errup = 0.0272753; errdown = 0.0267416;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.627949; errup = 0.00736333; errdown = 0.00742228;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.764423; errup = 0.021315; errdown = 0.0226577;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.816327; errup = 0.0328827; errdown = 0.0376241;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.854545; errup = 0.0487151; errdown = 0.06388;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.480638; errup = 0.0250011; errdown = 0.0249089;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.68834; errup = 0.00546682; errdown = 0.00551914;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.806452; errup = 0.0141336; errdown = 0.0149273;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.832714; errup = 0.023253; errdown = 0.0259226;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.811966; errup = 0.037304; errdown = 0.0432346;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.592466; errup = 0.0300791; errdown = 0.0307487;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.764847; errup = 0.00559549; errdown = 0.0056878;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.862416; errup = 0.0143005; errdown = 0.0155877;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.88785; errup = 0.0218923; errdown = 0.025898;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.857143; errup = 0.0349593; errdown = 0.0426621;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.645631; errup = 0.0349094; errdown = 0.0364216;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.82982; errup = 0.00563359; errdown = 0.00578284;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.90566; errup = 0.0122204; errdown = 0.0137092;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.871921; errup = 0.023876; errdown = 0.0279125;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.908046; errup = 0.0311873; errdown = 0.0422192;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.78169; errup = 0.0359018; errdown = 0.0402526;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.858189; errup = 0.00578039; errdown = 0.00597919;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.940397; errup = 0.011173; errdown = 0.0133076;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.941176; errup = 0.0180294; errdown = 0.024042;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.936709; errup = 0.0270672; errdown = 0.0405476;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.817308; errup = 0.0391317; errdown = 0.0459649;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.904392; errup = 0.00543304; errdown = 0.00571476;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.945626; errup = 0.0110594; errdown = 0.0133891;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.946667; errup = 0.0182461; errdown = 0.0252511;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.932432; errup = 0.028876; errdown = 0.0431238;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.825688; errup = 0.0374493; errdown = 0.0441394;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.944826; errup = 0.00325468; errdown = 0.00344007;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.97401; errup = 0.00558092; errdown = 0.00688314;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.973988; errup = 0.00845747; errdown = 0.0116538;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.956522; errup = 0.0171125; errdown = 0.0250664;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.797101; errup = 0.0504208; errdown = 0.0602257;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.97125; errup = 0.00311311; errdown = 0.00345772;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.982759; errup = 0.0050968; errdown = 0.00683869;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.983122; errup = 0.00806029; errdown = 0.0131435;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.968085; errup = 0.0173196; errdown = 0.0300713;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 9999) {eff = 0.897959; errup = 0.0433605; errdown = 0.0631366;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.98731; errup = 0.00212537; errdown = 0.00250879;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.988889; errup = 0.0036249; errdown = 0.00503323;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.993399; errup = 0.00426229; errdown = 0.00863929;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0118835;}
  return eff;
});

const NamedFunc get_0l_fakemet_trigeff2018("get_0l_fakemet_trigeff2018", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.783069; errup = 0.0136199; errdown = 0.0142422;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.76465; errup = 0.00349924; errdown = 0.00353528;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.56146; errup = 0.0108663; errdown = 0.0109243;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.624842; errup = 0.00339922; errdown = 0.00341149;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.466008; errup = 0.00230628; errdown = 0.00230483;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.803213; errup = 0.0127787; errdown = 0.0134115;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.772851; errup = 0.00347832; errdown = 0.00351597;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.619071; errup = 0.0110367; errdown = 0.0111582;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.666379; errup = 0.00348106; errdown = 0.00349915;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.509226; errup = 0.0025102; errdown = 0.00251067;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.823529; errup = 0.0121963; errdown = 0.0128678;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.777045; errup = 0.0035589; errdown = 0.00359946;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.641411; errup = 0.0117358; errdown = 0.0119034;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.700949; errup = 0.00358966; errdown = 0.00361433;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.540696; errup = 0.00271578; errdown = 0.00271818;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.797505; errup = 0.0126284; errdown = 0.0132211;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.78555; errup = 0.00355616; errdown = 0.0035991;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.665031; errup = 0.0118898; errdown = 0.0120975;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.723216; errup = 0.0036945; errdown = 0.00372492;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.584946; errup = 0.00287986; errdown = 0.00288564;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.790654; errup = 0.0126188; errdown = 0.0131822;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.787925; errup = 0.00359856; errdown = 0.00364326;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.687737; errup = 0.0118419; errdown = 0.0120856;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.756919; errup = 0.00375845; errdown = 0.00379792;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.614853; errup = 0.00306386; errdown = 0.00307294;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.802715; errup = 0.012136; errdown = 0.0127043;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.793769; errup = 0.00363703; errdown = 0.0036846;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.72912; errup = 0.0123851; errdown = 0.0127404;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.777236; errup = 0.00383485; errdown = 0.00388201;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.651016; errup = 0.00324028; errdown = 0.0032542;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.810811; errup = 0.011716; errdown = 0.0122781;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.797007; errup = 0.00369439; errdown = 0.00374461;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.731029; errup = 0.0123471; errdown = 0.0127048;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.80504; errup = 0.00381664; errdown = 0.00387342;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.678762; errup = 0.0034433; errdown = 0.00346269;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.81283; errup = 0.0117133; errdown = 0.0122836;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.797124; errup = 0.00377134; errdown = 0.00382372;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.756178; errup = 0.0125133; errdown = 0.0129493;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.81854; errup = 0.00389061; errdown = 0.00395578;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.718694; errup = 0.00354248; errdown = 0.0035696;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.817032; errup = 0.0113684; errdown = 0.0119228;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.802714; errup = 0.00377713; errdown = 0.00383181;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.781277; errup = 0.0124058; errdown = 0.0129153;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.840313; errup = 0.00388819; errdown = 0.00396528;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.741856; errup = 0.00368455; errdown = 0.00371883;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.825083; errup = 0.0110439; errdown = 0.0116004;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.80674; errup = 0.00382727; errdown = 0.00388508;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.815399; errup = 0.011973; errdown = 0.0125808;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.857696; errup = 0.00386119; errdown = 0.00394919;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.762597; errup = 0.00384232; errdown = 0.00388518;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.822276; errup = 0.00773276; errdown = 0.00799867;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.804838; errup = 0.00276043; errdown = 0.00279007;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.834778; errup = 0.0084645; errdown = 0.00881621;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.879252; errup = 0.00270911; errdown = 0.0027619;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.800344; errup = 0.00281079; errdown = 0.00284054;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.814829; errup = 0.00756296; errdown = 0.00780333;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.814466; errup = 0.00280749; errdown = 0.00284038;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.857221; errup = 0.00825338; errdown = 0.00865716;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.895166; errup = 0.00278813; errdown = 0.00285418;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.843622; errup = 0.00287179; errdown = 0.00291493;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.807988; errup = 0.00760664; errdown = 0.00783779;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.814443; errup = 0.00290636; errdown = 0.0029416;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.865638; errup = 0.00822363; errdown = 0.00865626;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.918065; errup = 0.00264745; errdown = 0.00272631;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.872208; errup = 0.00293961; errdown = 0.0029977;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.828114; errup = 0.00701563; errdown = 0.0072444;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.812111; errup = 0.00304631; errdown = 0.00308437;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.912081; errup = 0.00696799; errdown = 0.00748204;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.928379; errup = 0.00268497; errdown = 0.00277924;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.889586; errup = 0.00306626; errdown = 0.00314151;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.828745; errup = 0.00718099; errdown = 0.0074219;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.808351; errup = 0.0031814; errdown = 0.00322179;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.906933; errup = 0.00730866; errdown = 0.0078387;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.940203; errup = 0.0026612; errdown = 0.00277416;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.911654; errup = 0.00309177; errdown = 0.00319074;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.82974; errup = 0.00445576; errdown = 0.00454893;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.805113; errup = 0.00221366; errdown = 0.00223275;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.931136; errup = 0.00416028; errdown = 0.00439877;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.951112; errup = 0.00170434; errdown = 0.00176156;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.930972; errup = 0.00201305; errdown = 0.00206803;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.839123; errup = 0.00470632; errdown = 0.00481832;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.808875; errup = 0.00254223; errdown = 0.00256811;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.950015; errup = 0.00378067; errdown = 0.00406072;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.963379; errup = 0.00172179; errdown = 0.00180134;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.94979; errup = 0.0021136; errdown = 0.0021994;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 350) {eff = 0.85927; errup = 0.0037519; errdown = 0.00383613;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.807335; errup = 0.00232967; errdown = 0.00235116;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.961112; errup = 0.00252941; errdown = 0.00269201;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.968023; errup = 0.00143163; errdown = 0.0014949;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.965564; errup = 0.00163386; errdown = 0.00171025;}
  else if (ht> 0 && ht<= 200 && met> 350 && met<= 400) {eff = 0.876007; errup = 0.00517695; errdown = 0.00536483;}
  else if (ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.807657; errup = 0.00344952; errdown = 0.00349677;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.972668; errup = 0.00240389; errdown = 0.0026182;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.975655; errup = 0.00170611; errdown = 0.00182667;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.974032; errup = 0.00197332; errdown = 0.00212467;}
  else if (ht> 0 && ht<= 200 && met> 400 && met<= 450) {eff = 0.896953; errup = 0.00720465; errdown = 0.00766138;}
  else if (ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.832974; errup = 0.00483036; errdown = 0.00494271;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.981407; errup = 0.00237668; errdown = 0.00269446;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.983543; errup = 0.00192253; errdown = 0.0021567;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.972542; errup = 0.00281052; errdown = 0.0031042;}
  else if (ht> 0 && ht<= 200 && met> 450 && met<= 500) {eff = 0.925876; errup = 0.00968138; errdown = 0.0109088;}
  else if (ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.839387; errup = 0.00698255; errdown = 0.00723044;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.983294; errup = 0.00279405; errdown = 0.00329469;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.984362; errup = 0.0024794; errdown = 0.00289914;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.982897; errup = 0.00290135; errdown = 0.00342911;}
  else if (ht> 0 && ht<= 200 && met> 500 && met<= 9999) {eff = 0.918994; errup = 0.0145423; errdown = 0.0171055;}
  else if (ht> 200 && ht<= 600 && met> 500 && met<= 9999) {eff = 0.842346; errup = 0.0101654; errdown = 0.0107064;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.984349; errup = 0.00354204; errdown = 0.00443195;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.987952; errup = 0.00288307; errdown = 0.00365891;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.983651; errup = 0.00379901; errdown = 0.0047816;}
  return eff;
});

const NamedFunc get_1el_trigeff2018("get_1el_trigeff2018", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = b.el_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (el_pt> 20 && el_pt<= 25 && met> 150 && met<= 110) {eff = 0.0613195; errup = 0.00423119; errdown = 0.00398277;}
  else if (el_pt> 25 && el_pt<= 30 && met> 150 && met<= 110) {eff = 0.065407; errup = 0.00505952; errdown = 0.00473156;}
  else if (el_pt> 30 && el_pt<= 110 && met> 150 && met<= 110) {eff = 0.503611; errup = 0.0033815; errdown = 0.00338183;}
  else if (el_pt> 110 && el_pt<= 120 && met> 150 && met<= 110) {eff = 0.565506; errup = 0.0118904; errdown = 0.0119646;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 110) {eff = 0.833754; errup = 0.00543067; errdown = 0.00557368;}
  else if (el_pt> 20 && el_pt<= 25 && met> 155 && met<= 120) {eff = 0.294118; errup = 0.035356; errdown = 0.0331892;}
  else if (el_pt> 25 && el_pt<= 30 && met> 155 && met<= 120) {eff = 0.239766; errup = 0.0371352; errdown = 0.0338313;}
  else if (el_pt> 30 && el_pt<= 110 && met> 155 && met<= 120) {eff = 0.638531; errup = 0.0124216; errdown = 0.0126048;}
  else if (el_pt> 110 && el_pt<= 120 && met> 155 && met<= 120) {eff = 0.680851; errup = 0.0412595; errdown = 0.0440438;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 155 && met<= 120) {eff = 0.846939; errup = 0.0185027; errdown = 0.020394;}
  else if (el_pt> 20 && el_pt<= 25 && met> 160 && met<= 130) {eff = 0.34375; errup = 0.0376574; errdown = 0.0359119;}
  else if (el_pt> 25 && el_pt<= 30 && met> 160 && met<= 130) {eff = 0.320261; errup = 0.0421267; errdown = 0.0395842;}
  else if (el_pt> 30 && el_pt<= 110 && met> 160 && met<= 130) {eff = 0.6659; errup = 0.0133043; errdown = 0.0135659;}
  else if (el_pt> 110 && el_pt<= 120 && met> 160 && met<= 130) {eff = 0.744526; errup = 0.0388247; errdown = 0.0427257;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 160 && met<= 130) {eff = 0.870769; errup = 0.0189029; errdown = 0.0213686;}
  else if (el_pt> 20 && el_pt<= 25 && met> 165 && met<= 140) {eff = 0.493506; errup = 0.0434277; errdown = 0.0433369;}
  else if (el_pt> 25 && el_pt<= 30 && met> 165 && met<= 140) {eff = 0.433333; errup = 0.0497845; errdown = 0.0485776;}
  else if (el_pt> 30 && el_pt<= 110 && met> 165 && met<= 140) {eff = 0.72155; errup = 0.0129487; errdown = 0.0133175;}
  else if (el_pt> 110 && el_pt<= 120 && met> 165 && met<= 140) {eff = 0.724771; errup = 0.0448899; errdown = 0.0494292;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 165 && met<= 140) {eff = 0.847896; errup = 0.0208172; errdown = 0.0232425;}
  else if (el_pt> 20 && el_pt<= 25 && met> 170 && met<= 150) {eff = 0.529851; errup = 0.046418; errdown = 0.0468998;}
  else if (el_pt> 25 && el_pt<= 30 && met> 170 && met<= 150) {eff = 0.482759; errup = 0.0591637; errdown = 0.0587281;}
  else if (el_pt> 30 && el_pt<= 110 && met> 170 && met<= 150) {eff = 0.782133; errup = 0.0126483; errdown = 0.0131812;}
  else if (el_pt> 110 && el_pt<= 120 && met> 170 && met<= 150) {eff = 0.804348; errup = 0.0428773; errdown = 0.0503096;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 170 && met<= 150) {eff = 0.861199; errup = 0.0197526; errdown = 0.0222118;}
  else if (el_pt> 20 && el_pt<= 25 && met> 175 && met<= 160) {eff = 0.574074; errup = 0.0512024; errdown = 0.0526986;}
  else if (el_pt> 25 && el_pt<= 30 && met> 175 && met<= 160) {eff = 0.663043; errup = 0.0524528; errdown = 0.0563602;}
  else if (el_pt> 30 && el_pt<= 110 && met> 175 && met<= 160) {eff = 0.817901; errup = 0.0125485; errdown = 0.0132294;}
  else if (el_pt> 110 && el_pt<= 120 && met> 175 && met<= 160) {eff = 0.844828; errup = 0.0344943; errdown = 0.0411643;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 160) {eff = 0.919708; errup = 0.0165531; errdown = 0.0199497;}
  else if (el_pt> 20 && el_pt<= 25 && met> 180 && met<= 170) {eff = 0.654545; errup = 0.0481134; errdown = 0.0511868;}
  else if (el_pt> 25 && el_pt<= 30 && met> 180 && met<= 170) {eff = 0.792683; errup = 0.0465837; errdown = 0.0546294;}
  else if (el_pt> 30 && el_pt<= 110 && met> 180 && met<= 170) {eff = 0.852136; errup = 0.011889; errdown = 0.0126959;}
  else if (el_pt> 110 && el_pt<= 120 && met> 180 && met<= 170) {eff = 0.90099; errup = 0.0300397; errdown = 0.0392408;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 180 && met<= 170) {eff = 0.948339; errup = 0.0134475; errdown = 0.0172017;}
  else if (el_pt> 20 && el_pt<= 25 && met> 185 && met<= 180) {eff = 0.771084; errup = 0.0481984; errdown = 0.0555253;}
  else if (el_pt> 25 && el_pt<= 30 && met> 185 && met<= 180) {eff = 0.851351; errup = 0.0424361; errdown = 0.053383;}
  else if (el_pt> 30 && el_pt<= 110 && met> 185 && met<= 180) {eff = 0.897801; errup = 0.0109971; errdown = 0.0120867;}
  else if (el_pt> 110 && el_pt<= 120 && met> 185 && met<= 180) {eff = 0.885057; errup = 0.0347292; errdown = 0.0449884;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 185 && met<= 180) {eff = 0.953586; errup = 0.0136201; errdown = 0.0180267;}
  else if (el_pt> 20 && el_pt<= 25 && met> 190 && met<= 190) {eff = 0.804878; errup = 0.0454168; errdown = 0.0538207;}
  else if (el_pt> 25 && el_pt<= 30 && met> 190 && met<= 190) {eff = 0.848485; errup = 0.0453331; errdown = 0.0575786;}
  else if (el_pt> 30 && el_pt<= 110 && met> 190 && met<= 190) {eff = 0.919376; errup = 0.00988977; errdown = 0.0110519;}
  else if (el_pt> 110 && el_pt<= 120 && met> 190 && met<= 190) {eff = 0.927083; errup = 0.0265003; errdown = 0.0370657;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 190 && met<= 190) {eff = 0.979167; errup = 0.00896988; errdown = 0.0138476;}
  else if (el_pt> 20 && el_pt<= 25 && met> 195 && met<= 200) {eff = 0.783333; errup = 0.0556641; errdown = 0.0664662;}
  else if (el_pt> 25 && el_pt<= 30 && met> 195 && met<= 200) {eff = 0.888889; errup = 0.0400561; errdown = 0.0547111;}
  else if (el_pt> 30 && el_pt<= 110 && met> 195 && met<= 200) {eff = 0.931298; errup = 0.00994385; errdown = 0.0113613;}
  else if (el_pt> 110 && el_pt<= 120 && met> 195 && met<= 200) {eff = 1; errup = 0; errdown = 0.0287997;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 195 && met<= 200) {eff = 0.969565; errup = 0.011157; errdown = 0.0160084;}
  else if (el_pt> 20 && el_pt<= 25 && met> 200 && met<= 210) {eff = 0.880597; errup = 0.0402397; errdown = 0.0536015;}
  else if (el_pt> 25 && el_pt<= 30 && met> 200 && met<= 210) {eff = 0.945455; errup = 0.0295393; errdown = 0.0502232;}
  else if (el_pt> 30 && el_pt<= 110 && met> 200 && met<= 210) {eff = 0.951299; errup = 0.00869555; errdown = 0.0102965;}
  else if (el_pt> 110 && el_pt<= 120 && met> 200 && met<= 210) {eff = 0.970149; errup = 0.0192575; errdown = 0.0380185;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 200 && met<= 210) {eff = 0.981395; errup = 0.00888312; errdown = 0.0144657;}
  else if (el_pt> 20 && el_pt<= 25 && met> 210 && met<= 9999) {eff = 0.956044; errup = 0.0150689; errdown = 0.0209638;}
  else if (el_pt> 25 && el_pt<= 30 && met> 210 && met<= 9999) {eff = 0.927273; errup = 0.0203074; errdown = 0.0262392;}
  else if (el_pt> 30 && el_pt<= 110 && met> 210 && met<= 9999) {eff = 0.973943; errup = 0.003533; errdown = 0.00403137;}
  else if (el_pt> 110 && el_pt<= 120 && met> 210 && met<= 9999) {eff = 0.980237; errup = 0.0085104; errdown = 0.0131481;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 210 && met<= 9999) {eff = 0.991837; errup = 0.00323336; errdown = 0.00484391;}
  return eff;
});

const NamedFunc get_1mu_trigeff2018("get_1mu_trigeff2018", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = b.mu_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 110) {eff = 0.134868; errup = 0.00589481; errdown = 0.00568981;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 110) {eff = 0.491734; errup = 0.0099961; errdown = 0.00998962;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 110) {eff = 0.676738; errup = 0.00566981; errdown = 0.00572157;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 110) {eff = 0.932745; errup = 0.00254411; errdown = 0.00263477;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 155 && met<= 120) {eff = 0.392308; errup = 0.0325795; errdown = 0.0317007;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 155 && met<= 120) {eff = 0.738889; errup = 0.0339931; errdown = 0.0368625;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 155 && met<= 120) {eff = 0.83156; errup = 0.016017; errdown = 0.0172575;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 155 && met<= 120) {eff = 0.966785; errup = 0.00617036; errdown = 0.00738353;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 160 && met<= 130) {eff = 0.404545; errup = 0.0357413; errdown = 0.0348166;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 160 && met<= 130) {eff = 0.784884; errup = 0.0323551; errdown = 0.0359617;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 160 && met<= 130) {eff = 0.874187; errup = 0.015658; errdown = 0.0173931;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 160 && met<= 130) {eff = 0.965426; errup = 0.00666058; errdown = 0.00802104;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 165 && met<= 140) {eff = 0.588571; errup = 0.0393922; errdown = 0.0404783;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 165 && met<= 140) {eff = 0.806818; errup = 0.0306292; errdown = 0.0344349;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 165 && met<= 140) {eff = 0.90583; errup = 0.01397; errdown = 0.0159349;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 165 && met<= 140) {eff = 0.981162; errup = 0.00533576; errdown = 0.00706507;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 170 && met<= 150) {eff = 0.67; errup = 0.0347494; errdown = 0.0365714;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 170 && met<= 150) {eff = 0.856061; errup = 0.0312558; errdown = 0.0372898;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 170 && met<= 150) {eff = 0.936842; errup = 0.01254; errdown = 0.0150776;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 170 && met<= 150) {eff = 0.976705; errup = 0.0061118; errdown = 0.00791103;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 160) {eff = 0.732824; errup = 0.040383; errdown = 0.0442694;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 160) {eff = 0.813559; errup = 0.0370102; errdown = 0.0429202;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 160) {eff = 0.94864; errup = 0.0121479; errdown = 0.0151907;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 160) {eff = 0.989111; errup = 0.00431103; errdown = 0.00644732;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 180 && met<= 170) {eff = 0.753846; errup = 0.0393436; errdown = 0.0436253;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 180 && met<= 170) {eff = 0.931373; errup = 0.0249636; errdown = 0.0350079;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 180 && met<= 170) {eff = 0.934132; errup = 0.013645; errdown = 0.0165261;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 180 && met<= 170) {eff = 0.986275; errup = 0.0050482; errdown = 0.00731462;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 185 && met<= 180) {eff = 0.810526; errup = 0.0416155; errdown = 0.0489577;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 185 && met<= 180) {eff = 0.964706; errup = 0.0191475; errdown = 0.0331419;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 185 && met<= 180) {eff = 0.957529; errup = 0.0124757; errdown = 0.0165429;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 185 && met<= 180) {eff = 0.99115; errup = 0.00423058; errdown = 0.00694183;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 190 && met<= 190) {eff = 0.934066; errup = 0.0258435; errdown = 0.0373058;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 190 && met<= 190) {eff = 0.988372; errup = 0.00962116; errdown = 0.0262292;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 190 && met<= 190) {eff = 0.959459; errup = 0.011418; errdown = 0.0149723;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 190 && met<= 190) {eff = 0.995556; errup = 0.00287019; errdown = 0.00583174;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 195 && met<= 200) {eff = 0.894737; errup = 0.0355919; errdown = 0.0478087;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 195 && met<= 200) {eff = 0.956522; errup = 0.0235699; errdown = 0.0404889;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 195 && met<= 200) {eff = 0.995745; errup = 0.00352047; errdown = 0.0097167;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 195 && met<= 200) {eff = 0.99723; errup = 0.00229166; errdown = 0.00634082;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 200 && met<= 210) {eff = 0.96; errup = 0.0190465; errdown = 0.0304977;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 200 && met<= 210) {eff = 0.972603; errup = 0.0176765; errdown = 0.0349951;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 200 && met<= 210) {eff = 0.990868; errup = 0.00589653; errdown = 0.0119178;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 200 && met<= 210) {eff = 0.994667; errup = 0.0034441; errdown = 0.00699085;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 210 && met<= 9999) {eff = 0.954545; errup = 0.013342; errdown = 0.0176666;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 210 && met<= 9999) {eff = 0.978947; errup = 0.0100488; errdown = 0.0163327;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 210 && met<= 9999) {eff = 0.992614; errup = 0.00318649; errdown = 0.00496534;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 210 && met<= 9999) {eff = 0.997575; errup = 0.00131961; errdown = 0.00235344;}
  return eff;
});

const NamedFunc get_2el_trigeff2018("get_2el_trigeff2018", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = b.el_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.832117; errup = 0.0319327; errdown = 0.0319327;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.868687; errup = 0.0240023; errdown = 0.0240023;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.929825; errup = 0.0169171; errdown = 0.0169171;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.881279; errup = 0.0218574; errdown = 0.0218574;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.821577; errup = 0.0246627; errdown = 0.0246627;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.865922; errup = 0.0254678; errdown = 0.0254678;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.883621; errup = 0.0210536; errdown = 0.0210536;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.848168; errup = 0.0259661; errdown = 0.0259661;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.853659; errup = 0.0275997; errdown = 0.0275997;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.832258; errup = 0.0300112; errdown = 0.0300112;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.923077; errup = 0.0222833; errdown = 0.0222833;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.897638; errup = 0.0268979; errdown = 0.0268979;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.899329; errup = 0.0246501; errdown = 0.0246501;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.901961; errup = 0.0294438; errdown = 0.0294438;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.991383; errup = 0.00137385; errdown = 0.00137385;}
  return eff;
});

const NamedFunc get_2mu_trigeff2018("get_2mu_trigeff2018", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = b.mu_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.965812; errup = 0.0118789; errdown = 0.0118789;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.935754; errup = 0.0129587; errdown = 0.0129587;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.984612; errup = 0.000361107; errdown = 0.000361107;}
  return eff;
});

}
