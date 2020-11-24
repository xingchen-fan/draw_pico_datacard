#include <vector>
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/hig_utilities.hpp"
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
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.202703; errup = 0.0577474; errdown = 0.0486331;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.243401; errup = 0.00516343; errdown = 0.00509116;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.289245; errup = 0.00570872; errdown = 0.00564355;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.330854; errup = 0.00511006; errdown = 0.0050709;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.408515; errup = 0.0027316; errdown = 0.002726;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.297872; errup = 0.0814059; errdown = 0.0715756;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.30406; errup = 0.00619349; errdown = 0.00612423;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.335865; errup = 0.00655117; errdown = 0.00648949;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.368558; errup = 0.00574175; errdown = 0.00570519;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.449374; errup = 0.00301222; errdown = 0.00300853;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.351351; errup = 0.0952483; errdown = 0.0860018;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.347403; errup = 0.00698388; errdown = 0.00691979;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.389544; errup = 0.00740883; errdown = 0.00735894;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.40968; errup = 0.00627806; errdown = 0.00624912;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.479847; errup = 0.0032744; errdown = 0.00327268;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.34375; errup = 0.103529; errdown = 0.0921925;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.381182; errup = 0.0081055; errdown = 0.00804092;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.442154; errup = 0.00811501; errdown = 0.00808473;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.463081; errup = 0.00702641; errdown = 0.00701198;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.521223; errup = 0.00352075; errdown = 0.00352284;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.375; errup = 0.122914; errdown = 0.110663;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.469179; errup = 0.00908694; errdown = 0.00906691;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.477634; errup = 0.00881021; errdown = 0.00879655;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.489059; errup = 0.00750549; errdown = 0.00750063;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.556443; errup = 0.00377733; errdown = 0.00378382;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.466667; errup = 0.15827; errdown = 0.152935;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.496869; errup = 0.0097792; errdown = 0.00977685;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.554628; errup = 0.00962846; errdown = 0.00966887;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.526151; errup = 0.00817664; errdown = 0.00819047;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.591429; errup = 0.00405927; errdown = 0.00407168;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.647059; errup = 0.129922; errdown = 0.150722;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.562444; errup = 0.0106692; errdown = 0.0107261;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.584582; errup = 0.0103731; errdown = 0.010447;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.570557; errup = 0.00855781; errdown = 0.00859949;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.615664; errup = 0.00435457; errdown = 0.00437302;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.714286; errup = 0.13123; errdown = 0.168877;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.621965; errup = 0.0115951; errdown = 0.0117329;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.626667; errup = 0.0111406; errdown = 0.0112735;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.577571; errup = 0.00930451; errdown = 0.00935887;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.653083; errup = 0.00454038; errdown = 0.00456813;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.636364; errup = 0.164762; errdown = 0.195231;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.641304; errup = 0.0119969; errdown = 0.0121718;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.657925; errup = 0.011646; errdown = 0.0118347;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.612853; errup = 0.00979111; errdown = 0.00988133;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.688074; errup = 0.00471269; errdown = 0.00475151;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.875; errup = 0.103637; errdown = 0.23225;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.684962; errup = 0.0129659; errdown = 0.0132521;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.682345; errup = 0.0120156; errdown = 0.0122569;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.651388; errup = 0.010154; errdown = 0.0102904;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.708353; errup = 0.0049779; errdown = 0.00502779;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.791667; errup = 0.0868677; errdown = 0.116379;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.740741; errup = 0.00960556; errdown = 0.0098368;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.735731; errup = 0.00881812; errdown = 0.00900652;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.686861; errup = 0.00769244; errdown = 0.00779483;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.739948; errup = 0.00374889; errdown = 0.00378393;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.75; errup = 0.132707; errdown = 0.1849;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.803825; errup = 0.00997623; errdown = 0.0103627;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.78953; errup = 0.00957715; errdown = 0.00989845;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.714932; errup = 0.00820568; errdown = 0.00834736;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.78435; errup = 0.00396082; errdown = 0.00401366;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.875; errup = 0.103637; errdown = 0.23225;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.833603; errup = 0.010707; errdown = 0.0112662;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.825459; errup = 0.00982864; errdown = 0.01027;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.745188; errup = 0.00901608; errdown = 0.00922602;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.824385; errup = 0.00409059; errdown = 0.0041659;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.23126;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.85963; errup = 0.0115907; errdown = 0.01241;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.850206; errup = 0.0103481; errdown = 0.0109475;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.757143; errup = 0.0098043; errdown = 0.0100735;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.835609; errup = 0.00441926; errdown = 0.00451523;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.308024;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.902013; errup = 0.0109914; errdown = 0.0121363;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.865828; errup = 0.0111558; errdown = 0.0119586;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.784155; errup = 0.0104025; errdown = 0.0107677;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.866879; errup = 0.00454197; errdown = 0.00467439;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.940419; errup = 0.00674673; errdown = 0.00749814;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.896878; errup = 0.00722917; errdown = 0.00768865;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.813946; errup = 0.00718428; errdown = 0.00739968;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.889713; errup = 0.00306427; errdown = 0.00313953;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.959155; errup = 0.00743723; errdown = 0.00885222;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.934353; errup = 0.00746562; errdown = 0.00829373;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.82008; errup = 0.00864714; errdown = 0.00897445;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.917388; errup = 0.00326247; errdown = 0.00338145;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 350) {eff = 1; errup = 0; errdown = 0.841345;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.963964; errup = 0.00721984; errdown = 0.00875676;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.925249; errup = 0.00762407; errdown = 0.00836884;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.818534; errup = 0.00807502; errdown = 0.00835695;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.940378; errup = 0.00263708; errdown = 0.00274835;}
  else if (ht> 0 && ht<= 200 && met> 350 && met<= 400) {eff = 1.25382e-308; errup = 5.56355e-307; errdown = 1.78249e-306;}
  else if (ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.983471; errup = 0.0078941; errdown = 0.012876;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.892; errup = 0.0140398; errdown = 0.0157188;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.831216; errup = 0.0114216; errdown = 0.0120465;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.958167; errup = 0.00307583; errdown = 0.00329941;}
  else if (ht> 0 && ht<= 200 && met> 400 && met<= 450) {eff = 9.53547e-322; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.969388; errup = 0.0166146; errdown = 0.0288818;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.870833; errup = 0.022025; errdown = 0.0254039;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.876686; errup = 0.0146189; errdown = 0.0161656;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.956439; errup = 0.00445297; errdown = 0.00490815;}
  else if (ht> 0 && ht<= 200 && met> 450 && met<= 500) {eff = 1.0; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.961538; errup = 0.0318392; errdown = 0.0829559;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.858268; errup = 0.0316515; errdown = 0.0379788;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.864151; errup = 0.0214229; errdown = 0.0244106;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.970889; errup = 0.0047152; errdown = 0.00551772;}
  else if (ht> 0 && ht<= 200 && met> 500 && met<= 9999) {eff = 1.0; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 500 && met<= 9999) {eff = 0.9; errup = 0.082873; errdown = 0.194135;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.892857; errup = 0.0416676; errdown = 0.0585112;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.910828; errup = 0.0229624; errdown = 0.028904;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.975714; errup = 0.00579072; errdown = 0.00731565;}
  return eff;
});

const NamedFunc get_1el_trigeff2018("get_1el_trigeff2018", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 20 && el_pt<= 25 && met> 0 && met<= 110) {eff = 0.0613195; errup = 0.00423119; errdown = 0.00398277;}
  else if (el_pt> 25 && el_pt<= 30 && met> 0 && met<= 110) {eff = 0.065407; errup = 0.00505952; errdown = 0.00473156;}
  else if (el_pt> 30 && el_pt<= 110 && met> 0 && met<= 110) {eff = 0.503611; errup = 0.0033815; errdown = 0.00338183;}
  else if (el_pt> 110 && el_pt<= 120 && met> 0 && met<= 110) {eff = 0.565506; errup = 0.0118904; errdown = 0.0119646;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 110) {eff = 0.833754; errup = 0.00543067; errdown = 0.00557368;}
  else if (el_pt> 20 && el_pt<= 25 && met> 110 && met<= 120) {eff = 0.294118; errup = 0.035356; errdown = 0.0331892;}
  else if (el_pt> 25 && el_pt<= 30 && met> 110 && met<= 120) {eff = 0.239766; errup = 0.0371352; errdown = 0.0338313;}
  else if (el_pt> 30 && el_pt<= 110 && met> 110 && met<= 120) {eff = 0.638531; errup = 0.0124216; errdown = 0.0126048;}
  else if (el_pt> 110 && el_pt<= 120 && met> 110 && met<= 120) {eff = 0.680851; errup = 0.0412595; errdown = 0.0440438;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 110 && met<= 120) {eff = 0.846939; errup = 0.0185027; errdown = 0.020394;}
  else if (el_pt> 20 && el_pt<= 25 && met> 120 && met<= 130) {eff = 0.34375; errup = 0.0376574; errdown = 0.0359119;}
  else if (el_pt> 25 && el_pt<= 30 && met> 120 && met<= 130) {eff = 0.320261; errup = 0.0421267; errdown = 0.0395842;}
  else if (el_pt> 30 && el_pt<= 110 && met> 120 && met<= 130) {eff = 0.6659; errup = 0.0133043; errdown = 0.0135659;}
  else if (el_pt> 110 && el_pt<= 120 && met> 120 && met<= 130) {eff = 0.744526; errup = 0.0388247; errdown = 0.0427257;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 120 && met<= 130) {eff = 0.870769; errup = 0.0189029; errdown = 0.0213686;}
  else if (el_pt> 20 && el_pt<= 25 && met> 130 && met<= 140) {eff = 0.493506; errup = 0.0434277; errdown = 0.0433369;}
  else if (el_pt> 25 && el_pt<= 30 && met> 130 && met<= 140) {eff = 0.433333; errup = 0.0497845; errdown = 0.0485776;}
  else if (el_pt> 30 && el_pt<= 110 && met> 130 && met<= 140) {eff = 0.72155; errup = 0.0129487; errdown = 0.0133175;}
  else if (el_pt> 110 && el_pt<= 120 && met> 130 && met<= 140) {eff = 0.724771; errup = 0.0448899; errdown = 0.0494292;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 130 && met<= 140) {eff = 0.847896; errup = 0.0208172; errdown = 0.0232425;}
  else if (el_pt> 20 && el_pt<= 25 && met> 140 && met<= 150) {eff = 0.529851; errup = 0.046418; errdown = 0.0468998;}
  else if (el_pt> 25 && el_pt<= 30 && met> 140 && met<= 150) {eff = 0.482759; errup = 0.0591637; errdown = 0.0587281;}
  else if (el_pt> 30 && el_pt<= 110 && met> 140 && met<= 150) {eff = 0.782133; errup = 0.0126483; errdown = 0.0131812;}
  else if (el_pt> 110 && el_pt<= 120 && met> 140 && met<= 150) {eff = 0.804348; errup = 0.0428773; errdown = 0.0503096;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 140 && met<= 150) {eff = 0.861199; errup = 0.0197526; errdown = 0.0222118;}
  else if (el_pt> 20 && el_pt<= 25 && met> 150 && met<= 160) {eff = 0.574074; errup = 0.0512024; errdown = 0.0526986;}
  else if (el_pt> 25 && el_pt<= 30 && met> 150 && met<= 160) {eff = 0.663043; errup = 0.0524528; errdown = 0.0563602;}
  else if (el_pt> 30 && el_pt<= 110 && met> 150 && met<= 160) {eff = 0.817901; errup = 0.0125485; errdown = 0.0132294;}
  else if (el_pt> 110 && el_pt<= 120 && met> 150 && met<= 160) {eff = 0.844828; errup = 0.0344943; errdown = 0.0411643;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 160) {eff = 0.919708; errup = 0.0165531; errdown = 0.0199497;}
  else if (el_pt> 20 && el_pt<= 25 && met> 160 && met<= 170) {eff = 0.654545; errup = 0.0481134; errdown = 0.0511868;}
  else if (el_pt> 25 && el_pt<= 30 && met> 160 && met<= 170) {eff = 0.792683; errup = 0.0465837; errdown = 0.0546294;}
  else if (el_pt> 30 && el_pt<= 110 && met> 160 && met<= 170) {eff = 0.852136; errup = 0.011889; errdown = 0.0126959;}
  else if (el_pt> 110 && el_pt<= 120 && met> 160 && met<= 170) {eff = 0.90099; errup = 0.0300397; errdown = 0.0392408;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 160 && met<= 170) {eff = 0.948339; errup = 0.0134475; errdown = 0.0172017;}
  else if (el_pt> 20 && el_pt<= 25 && met> 170 && met<= 180) {eff = 0.771084; errup = 0.0481984; errdown = 0.0555253;}
  else if (el_pt> 25 && el_pt<= 30 && met> 170 && met<= 180) {eff = 0.851351; errup = 0.0424361; errdown = 0.053383;}
  else if (el_pt> 30 && el_pt<= 110 && met> 170 && met<= 180) {eff = 0.897801; errup = 0.0109971; errdown = 0.0120867;}
  else if (el_pt> 110 && el_pt<= 120 && met> 170 && met<= 180) {eff = 0.885057; errup = 0.0347292; errdown = 0.0449884;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 170 && met<= 180) {eff = 0.953586; errup = 0.0136201; errdown = 0.0180267;}
  else if (el_pt> 20 && el_pt<= 25 && met> 180 && met<= 190) {eff = 0.804878; errup = 0.0454168; errdown = 0.0538207;}
  else if (el_pt> 25 && el_pt<= 30 && met> 180 && met<= 190) {eff = 0.848485; errup = 0.0453331; errdown = 0.0575786;}
  else if (el_pt> 30 && el_pt<= 110 && met> 180 && met<= 190) {eff = 0.919376; errup = 0.00988977; errdown = 0.0110519;}
  else if (el_pt> 110 && el_pt<= 120 && met> 180 && met<= 190) {eff = 0.927083; errup = 0.0265003; errdown = 0.0370657;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 180 && met<= 190) {eff = 0.979167; errup = 0.00896988; errdown = 0.0138476;}
  else if (el_pt> 20 && el_pt<= 25 && met> 190 && met<= 200) {eff = 0.783333; errup = 0.0556641; errdown = 0.0664662;}
  else if (el_pt> 25 && el_pt<= 30 && met> 190 && met<= 200) {eff = 0.888889; errup = 0.0400561; errdown = 0.0547111;}
  else if (el_pt> 30 && el_pt<= 110 && met> 190 && met<= 200) {eff = 0.931298; errup = 0.00994385; errdown = 0.0113613;}
  else if (el_pt> 110 && el_pt<= 120 && met> 190 && met<= 200) {eff = 1; errup = 0; errdown = 0.0287997;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 190 && met<= 200) {eff = 0.969565; errup = 0.011157; errdown = 0.0160084;}
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
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 110) {eff = 0.134868; errup = 0.00589481; errdown = 0.00568981;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 110) {eff = 0.491734; errup = 0.0099961; errdown = 0.00998962;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 110) {eff = 0.676738; errup = 0.00566981; errdown = 0.00572157;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 110) {eff = 0.932745; errup = 0.00254411; errdown = 0.00263477;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 110 && met<= 120) {eff = 0.392308; errup = 0.0325795; errdown = 0.0317007;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 110 && met<= 120) {eff = 0.738889; errup = 0.0339931; errdown = 0.0368625;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 110 && met<= 120) {eff = 0.83156; errup = 0.016017; errdown = 0.0172575;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 110 && met<= 120) {eff = 0.966785; errup = 0.00617036; errdown = 0.00738353;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 120 && met<= 130) {eff = 0.404545; errup = 0.0357413; errdown = 0.0348166;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 120 && met<= 130) {eff = 0.784884; errup = 0.0323551; errdown = 0.0359617;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 120 && met<= 130) {eff = 0.874187; errup = 0.015658; errdown = 0.0173931;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 120 && met<= 130) {eff = 0.965426; errup = 0.00666058; errdown = 0.00802104;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 130 && met<= 140) {eff = 0.588571; errup = 0.0393922; errdown = 0.0404783;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 130 && met<= 140) {eff = 0.806818; errup = 0.0306292; errdown = 0.0344349;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 130 && met<= 140) {eff = 0.90583; errup = 0.01397; errdown = 0.0159349;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 130 && met<= 140) {eff = 0.981162; errup = 0.00533576; errdown = 0.00706507;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 140 && met<= 150) {eff = 0.67; errup = 0.0347494; errdown = 0.0365714;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 140 && met<= 150) {eff = 0.856061; errup = 0.0312558; errdown = 0.0372898;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 140 && met<= 150) {eff = 0.936842; errup = 0.01254; errdown = 0.0150776;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 140 && met<= 150) {eff = 0.976705; errup = 0.0061118; errdown = 0.00791103;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 160) {eff = 0.732824; errup = 0.040383; errdown = 0.0442694;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 160) {eff = 0.813559; errup = 0.0370102; errdown = 0.0429202;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 160) {eff = 0.94864; errup = 0.0121479; errdown = 0.0151907;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 160) {eff = 0.989111; errup = 0.00431103; errdown = 0.00644732;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 160 && met<= 170) {eff = 0.753846; errup = 0.0393436; errdown = 0.0436253;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 160 && met<= 170) {eff = 0.931373; errup = 0.0249636; errdown = 0.0350079;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 160 && met<= 170) {eff = 0.934132; errup = 0.013645; errdown = 0.0165261;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 160 && met<= 170) {eff = 0.986275; errup = 0.0050482; errdown = 0.00731462;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 170 && met<= 180) {eff = 0.810526; errup = 0.0416155; errdown = 0.0489577;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 170 && met<= 180) {eff = 0.964706; errup = 0.0191475; errdown = 0.0331419;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 170 && met<= 180) {eff = 0.957529; errup = 0.0124757; errdown = 0.0165429;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 170 && met<= 180) {eff = 0.99115; errup = 0.00423058; errdown = 0.00694183;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 180 && met<= 190) {eff = 0.934066; errup = 0.0258435; errdown = 0.0373058;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 180 && met<= 190) {eff = 0.988372; errup = 0.00962116; errdown = 0.0262292;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 180 && met<= 190) {eff = 0.959459; errup = 0.011418; errdown = 0.0149723;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 180 && met<= 190) {eff = 0.995556; errup = 0.00287019; errdown = 0.00583174;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 190 && met<= 200) {eff = 0.894737; errup = 0.0355919; errdown = 0.0478087;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 190 && met<= 200) {eff = 0.956522; errup = 0.0235699; errdown = 0.0404889;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 190 && met<= 200) {eff = 0.995745; errup = 0.00352047; errdown = 0.0097167;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 190 && met<= 200) {eff = 0.99723; errup = 0.00229166; errdown = 0.00634082;}
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
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
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
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.965812; errup = 0.0118789; errdown = 0.0118789;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.935754; errup = 0.0129587; errdown = 0.0129587;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.984612; errup = 0.000361107; errdown = 0.000361107;}
  return eff;
});

}
