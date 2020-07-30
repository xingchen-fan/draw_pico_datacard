#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/apply_trigeffs2016.hpp"

namespace Higfuncs{

const NamedFunc get_0l_trigeff2016("get_0l_trigeff2016", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.421454; errup = 0.0122918; errdown = 0.0121974;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.511723; errup = 0.00553945; errdown = 0.0055423;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.524272; errup = 0.0229051; errdown = 0.0230034;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.502959; errup = 0.0412691; errdown = 0.0413067;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.575342; errup = 0.0630717; errdown = 0.0653595;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.478716; errup = 0.0135556; errdown = 0.0135252;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.58184; errup = 0.00573815; errdown = 0.00576013;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.523909; errup = 0.0237337; errdown = 0.0238375;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.564706; errup = 0.0404303; errdown = 0.0412473;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.5; errup = 0.0653624; errdown = 0.0653624;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.537445; errup = 0.0151993; errdown = 0.0152672;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.641288; errup = 0.00593655; errdown = 0.00597961;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.629712; errup = 0.023519; errdown = 0.0241214;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.580645; errup = 0.0474312; errdown = 0.0488429;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.682927; errup = 0.0786376; errdown = 0.0888868;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.556346; errup = 0.0176245; errdown = 0.0177627;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.699494; errup = 0.00600351; errdown = 0.00607174;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.699752; errup = 0.0235238; errdown = 0.0245674;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.59854; errup = 0.0445755; errdown = 0.0461321;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.680851; errup = 0.0733593; errdown = 0.0821322;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.638268; errup = 0.0184422; errdown = 0.0188431;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.749666; errup = 0.00603513; errdown = 0.00613206;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.733918; errup = 0.0246003; errdown = 0.0260493;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.72807; errup = 0.0436702; errdown = 0.0480675;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.574468; errup = 0.0799773; errdown = 0.083562;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.654701; errup = 0.02022; errdown = 0.0207713;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.784363; errup = 0.00597389; errdown = 0.00609421;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.75641; errup = 0.0249947; errdown = 0.0267432;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.734513; errup = 0.0434929; errdown = 0.0480594;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.634615; errup = 0.0726695; errdown = 0.0785167;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.741722; errup = 0.0210861; errdown = 0.022209;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.822244; errup = 0.00596603; errdown = 0.00612396;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.819672; errup = 0.0225045; errdown = 0.0247522;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.798319; errup = 0.0380713; errdown = 0.0436274;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.702703; errup = 0.0810427; errdown = 0.0937327;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.645533; errup = 0.0266391; errdown = 0.0275234;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.853995; errup = 0.00578957; errdown = 0.00598182;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.832753; errup = 0.0224984; errdown = 0.0249952;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.826531; errup = 0.0394452; errdown = 0.0469475;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.771429; errup = 0.0749545; errdown = 0.0932188;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.714286; errup = 0.0269016; errdown = 0.0284131;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.874115; errup = 0.00585518; errdown = 0.00609152;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.849817; errup = 0.022045; errdown = 0.0248186;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.818182; errup = 0.0425143; errdown = 0.0506834;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.880952; errup = 0.0504413; errdown = 0.0725148;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.75; errup = 0.0307629; errdown = 0.0333008;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.903775; errup = 0.00541855; errdown = 0.00569667;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.876068; errup = 0.0218991; errdown = 0.0254228;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.835443; errup = 0.042987; errdown = 0.0526851;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.870968; errup = 0.0607061; errdown = 0.0903257;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.794118; errup = 0.0236918; errdown = 0.0257459;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.926773; errup = 0.00363813; errdown = 0.00380788;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.915074; errup = 0.0129577; errdown = 0.0148636;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.9; errup = 0.0232988; errdown = 0.0285863;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.830769; errup = 0.0480266; errdown = 0.0597643;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.806122; errup = 0.0290347; errdown = 0.0324296;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.955567; errup = 0.00326322; errdown = 0.00349935;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.945026; errup = 0.0116963; errdown = 0.0142838;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.938356; errup = 0.0198781; errdown = 0.0268833;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.941176; errup = 0.0318434; errdown = 0.0539242;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.860656; errup = 0.0320534; errdown = 0.0387031;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.96675; errup = 0.00316544; errdown = 0.00346976;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.967262; errup = 0.00964044; errdown = 0.0128421;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.949153; errup = 0.0199859; errdown = 0.0291361;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.934426; errup = 0.0311194; errdown = 0.0488167;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.844444; errup = 0.0392546; errdown = 0.0479576;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.976041; errup = 0.00303115; errdown = 0.00342929;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.964968; errup = 0.0103099; errdown = 0.0137191;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.953488; errup = 0.0221284; errdown = 0.0352492;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.0449824;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.857143; errup = 0.0511428; errdown = 0.0684449;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.982009; errup = 0.00296575; errdown = 0.00348796;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.978102; errup = 0.00865233; errdown = 0.0128502;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.0180628;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.964286; errup = 0.0295635; errdown = 0.077387;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.873239; errup = 0.0402311; errdown = 0.0525007;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.987697; errup = 0.00186062; errdown = 0.00216111;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.986667; errup = 0.00490435; errdown = 0.0071078;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.995238; errup = 0.00393961; errdown = 0.0108643;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.0239329;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.870968; errup = 0.0607061; errdown = 0.0903257;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.992795; errup = 0.00183562; errdown = 0.00237029;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.99505; errup = 0.00319693; errdown = 0.00649192;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.984848; errup = 0.0097805; errdown = 0.0196341;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.983871; errup = 0.0133466; errdown = 0.0361115;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 9999) {eff = 0.846154; errup = 0.0721243; errdown = 0.105033;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.995964; errup = 0.00139558; errdown = 0.00198449;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.998106; errup = 0.00156681; errdown = 0.00434157;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.983784; errup = 0.00881271; errdown = 0.0155222;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 0.987179; errup = 0.0106082; errdown = 0.0288622;}
  return eff;
});

const NamedFunc get_0l_fakemet_trigeff2016("get_0l_fakemet_trigeff2016", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.521739; errup = 0.121562; errdown = 0.123779;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.32681; errup = 0.0153293; errdown = 0.0149793;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.223389; errup = 0.00947246; errdown = 0.00920154;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.31678; errup = 0.00327772; errdown = 0.00325977;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.267906; errup = 0.00247311; errdown = 0.00245881;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.230769; errup = 0.174724; errdown = 0.122762;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.353439; errup = 0.0162316; errdown = 0.0159111;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.235838; errup = 0.0106501; errdown = 0.0103362;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.359084; errup = 0.00368919; errdown = 0.00367273;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.307724; errup = 0.00283001; errdown = 0.00281572;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.384615; errup = 0.176146; errdown = 0.154567;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.346341; errup = 0.0174095; errdown = 0.0170212;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.274639; errup = 0.0119181; errdown = 0.0116136;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.406089; errup = 0.0041357; errdown = 0.00412254;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.346227; errup = 0.00319217; errdown = 0.00317848;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.222222; errup = 0.221429; errdown = 0.142118;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.434049; errup = 0.0202664; errdown = 0.0200563;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.301009; errup = 0.0133184; errdown = 0.0130004;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.42991; errup = 0.00453432; errdown = 0.0045227;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.380039; errup = 0.00354209; errdown = 0.00352946;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.125; errup = 0.23225; errdown = 0.103637;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.406305; errup = 0.0215811; errdown = 0.0212393;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.305254; errup = 0.0144889; errdown = 0.0141248;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.464126; errup = 0.00491894; errdown = 0.00491204;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.422041; errup = 0.00393632; errdown = 0.00392652;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.3; errup = 0.208262; errdown = 0.158328;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.478992; errup = 0.0239675; errdown = 0.0238753;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.372116; errup = 0.0159345; errdown = 0.0156697;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.505893; errup = 0.00532649; errdown = 0.00532781;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.453633; errup = 0.00425521; errdown = 0.00424851;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.444444; errup = 0.213463; errdown = 0.198267;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.514925; errup = 0.0260974; errdown = 0.0261752;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.363305; errup = 0.0175341; errdown = 0.0171901;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.53869; errup = 0.00572019; errdown = 0.00573029;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.491062; errup = 0.004678; errdown = 0.00467645;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.25; errup = 0.368402; errdown = 0.207731;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.533333; errup = 0.0281638; errdown = 0.0283671;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.370945; errup = 0.0190226; errdown = 0.0186447;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.57002; errup = 0.00616272; errdown = 0.00618423;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.527627; errup = 0.00500085; errdown = 0.00500635;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0; errup = 0.841345; errdown = 0;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.49226; errup = 0.0293427; errdown = 0.0292922;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.420118; errup = 0.0198312; errdown = 0.0195858;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.598532; errup = 0.00654851; errdown = 0.0065834;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.56653; errup = 0.00533824; errdown = 0.00535356;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.75; errup = 0.207731; errdown = 0.368402;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.552239; errup = 0.0319712; errdown = 0.032384;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.398614; errup = 0.0214137; errdown = 0.0210477;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.621147; errup = 0.00684845; errdown = 0.00689639;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.593498; errup = 0.0056997; errdown = 0.0057247;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.517544; errup = 0.0244284; errdown = 0.0245089;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.421053; errup = 0.016128; errdown = 0.0159663;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.662695; errup = 0.00512347; errdown = 0.00516153;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.629123; errup = 0.00436216; errdown = 0.00438313;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.5; errup = 0.417248; errdown = 0.417248;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.554545; errup = 0.0286551; errdown = 0.0290033;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.463415; errup = 0.0176084; errdown = 0.0175206;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.706727; errup = 0.00552328; errdown = 0.00558401;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.67454; errup = 0.00471951; errdown = 0.00475482;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.142857; errup = 0.257124; errdown = 0.11848;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.59387; errup = 0.031869; errdown = 0.0326317;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.461111; errup = 0.0193141; errdown = 0.0192021;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.727979; errup = 0.00610267; errdown = 0.00618832;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.70598; errup = 0.00519272; errdown = 0.00524613;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.841345;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.631579; errup = 0.0367831; errdown = 0.0382673;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.444079; errup = 0.0210507; errdown = 0.0208594;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.753563; errup = 0.00658886; errdown = 0.00670748;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.755499; errup = 0.0054641; errdown = 0.00554673;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.458642;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.60221; errup = 0.0384238; errdown = 0.0396347;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.471545; errup = 0.0235569; errdown = 0.0234362;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.788895; errup = 0.00690836; errdown = 0.00707453;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.775937; errup = 0.00589877; errdown = 0.00600944;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.841345;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.66436; errup = 0.0288543; errdown = 0.0300599;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.443362; errup = 0.018007; errdown = 0.0178642;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.798824; errup = 0.00501706; errdown = 0.00511097;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.816055; errup = 0.00401406; errdown = 0.00408216;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.8; errup = 0.166039; errdown = 0.324251;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.557692; errup = 0.0424246; errdown = 0.0432206;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.503289; errup = 0.0210759; errdown = 0.0210872;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.840655; errup = 0.00580347; errdown = 0.00597619;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.85274; errup = 0.00454817; errdown = 0.00466528;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 350) {eff = 9.53547e-322; errup = 7.16395e-322; errdown = 7.16395e-322;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.581006; errup = 0.0390692; errdown = 0.0400394;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.487421; errup = 0.0206099; errdown = 0.0205688;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.857458; errup = 0.00522276; errdown = 0.00538384;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.875; errup = 0.00402788; errdown = 0.00414019;}
  else if (ht> 0 && ht<= 200 && met> 350 && met<= 400) {eff = 1.0; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.560976; errup = 0.0596747; errdown = 0.0613142;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.482072; errup = 0.0335364; errdown = 0.0333849;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.884615; errup = 0.00709896; errdown = 0.00748629;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.909173; errup = 0.00502405; errdown = 0.00527927;}
  else if (ht> 0 && ht<= 200 && met> 400 && met<= 450) {eff = 1.0; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.657895; errup = 0.084041; errdown = 0.093599;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.584071; errup = 0.0497397; errdown = 0.0513603;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.904222; errup = 0.00952096; errdown = 0.0103984;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.935305; errup = 0.0061342; errdown = 0.00669751;}
  else if (ht> 0 && ht<= 200 && met> 450 && met<= 500) {eff = 1.0; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.5; errup = 0.195182; errdown = 0.195182;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.688889; errup = 0.0743826; errdown = 0.0839837;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.93763; errup = 0.0110824; errdown = 0.0130713;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.941648; errup = 0.00796507; errdown = 0.00904797;}
  else if (ht> 0 && ht<= 200 && met> 500 && met<= 9999) {eff = 1.0; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 500 && met<= 9999) {eff = 0.666667; errup = 0.211935; errdown = 0.282023;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.9; errup = 0.0643201; errdown = 0.116971;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.958042; errup = 0.0118128; errdown = 0.0154801;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.948919; errup = 0.00978575; errdown = 0.0117251;}
  return eff;
});

const NamedFunc get_1el_trigeff2016("get_1el_trigeff2016", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = b.el_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (el_pt> 20 && el_pt<= 25 && met> 150 && met<= 110) {eff = 0.0582445; errup = 0.00513861; errdown = 0.00476062;}
  else if (el_pt> 25 && el_pt<= 30 && met> 150 && met<= 110) {eff = 0.157254; errup = 0.00903781; errdown = 0.00864611;}
  else if (el_pt> 30 && el_pt<= 110 && met> 150 && met<= 110) {eff = 0.530553; errup = 0.00429568; errdown = 0.00430018;}
  else if (el_pt> 110 && el_pt<= 120 && met> 150 && met<= 110) {eff = 0.650194; errup = 0.0151726; errdown = 0.015473;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 110) {eff = 0.89841; errup = 0.0059834; errdown = 0.00630236;}
  else if (el_pt> 20 && el_pt<= 25 && met> 155 && met<= 120) {eff = 0.149123; errup = 0.0411162; errdown = 0.0341911;}
  else if (el_pt> 25 && el_pt<= 30 && met> 155 && met<= 120) {eff = 0.338028; errup = 0.0652582; errdown = 0.0601709;}
  else if (el_pt> 30 && el_pt<= 110 && met> 155 && met<= 120) {eff = 0.672389; errup = 0.0181985; errdown = 0.0187112;}
  else if (el_pt> 110 && el_pt<= 120 && met> 155 && met<= 120) {eff = 0.686275; errup = 0.06986; errdown = 0.078157;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 155 && met<= 120) {eff = 0.969231; errup = 0.0146686; errdown = 0.02366;}
  else if (el_pt> 20 && el_pt<= 25 && met> 160 && met<= 130) {eff = 0.2; errup = 0.049646; errdown = 0.0425693;}
  else if (el_pt> 25 && el_pt<= 30 && met> 160 && met<= 130) {eff = 0.314286; errup = 0.0651353; errdown = 0.0592025;}
  else if (el_pt> 30 && el_pt<= 110 && met> 160 && met<= 130) {eff = 0.684039; errup = 0.0192486; errdown = 0.0198735;}
  else if (el_pt> 110 && el_pt<= 120 && met> 160 && met<= 130) {eff = 0.886364; errup = 0.0481931; errdown = 0.0695663;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 160 && met<= 130) {eff = 0.95082; errup = 0.0193366; errdown = 0.02822;}
  else if (el_pt> 20 && el_pt<= 25 && met> 165 && met<= 140) {eff = 0.425287; errup = 0.0593282; errdown = 0.0574387;}
  else if (el_pt> 25 && el_pt<= 30 && met> 165 && met<= 140) {eff = 0.513514; errup = 0.0641717; errdown = 0.0645758;}
  else if (el_pt> 30 && el_pt<= 110 && met> 165 && met<= 140) {eff = 0.723894; errup = 0.0192634; errdown = 0.0200923;}
  else if (el_pt> 110 && el_pt<= 120 && met> 165 && met<= 140) {eff = 0.909091; errup = 0.0429965; errdown = 0.0660622;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 165 && met<= 140) {eff = 0.892157; errup = 0.0311385; errdown = 0.0399966;}
  else if (el_pt> 20 && el_pt<= 25 && met> 170 && met<= 150) {eff = 0.453333; errup = 0.0643662; errdown = 0.0629896;}
  else if (el_pt> 25 && el_pt<= 30 && met> 170 && met<= 150) {eff = 0.464286; errup = 0.0755319; errdown = 0.0741026;}
  else if (el_pt> 30 && el_pt<= 110 && met> 170 && met<= 150) {eff = 0.765531; errup = 0.0193923; errdown = 0.0205114;}
  else if (el_pt> 110 && el_pt<= 120 && met> 170 && met<= 150) {eff = 0.844444; errup = 0.0555301; errdown = 0.0737002;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 170 && met<= 150) {eff = 0.947917; errup = 0.0223143; errdown = 0.0336975;}
  else if (el_pt> 20 && el_pt<= 25 && met> 175 && met<= 160) {eff = 0.606557; errup = 0.0681557; errdown = 0.0720644;}
  else if (el_pt> 25 && el_pt<= 30 && met> 175 && met<= 160) {eff = 0.731707; errup = 0.0738471; errdown = 0.0869218;}
  else if (el_pt> 30 && el_pt<= 110 && met> 175 && met<= 160) {eff = 0.816934; errup = 0.0188601; errdown = 0.0203981;}
  else if (el_pt> 110 && el_pt<= 120 && met> 175 && met<= 160) {eff = 0.765957; errup = 0.0651312; errdown = 0.0782067;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 160) {eff = 0.987013; errup = 0.010746; errdown = 0.029229;}
  else if (el_pt> 20 && el_pt<= 25 && met> 180 && met<= 170) {eff = 0.555556; errup = 0.0827174; errdown = 0.0855148;}
  else if (el_pt> 25 && el_pt<= 30 && met> 180 && met<= 170) {eff = 0.75; errup = 0.0815095; errdown = 0.0999181;}
  else if (el_pt> 30 && el_pt<= 110 && met> 180 && met<= 170) {eff = 0.836364; errup = 0.0192059; errdown = 0.0210707;}
  else if (el_pt> 110 && el_pt<= 120 && met> 180 && met<= 170) {eff = 0.965517; errup = 0.0285434; errdown = 0.0748731;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 180 && met<= 170) {eff = 1; errup = 0; errdown = 0.0287997;}
  else if (el_pt> 20 && el_pt<= 25 && met> 185 && met<= 180) {eff = 0.625; errup = 0.0842352; errdown = 0.0913836;}
  else if (el_pt> 25 && el_pt<= 30 && met> 185 && met<= 180) {eff = 0.729167; errup = 0.0683073; errdown = 0.0792509;}
  else if (el_pt> 30 && el_pt<= 110 && met> 185 && met<= 180) {eff = 0.869333; errup = 0.0176707; errdown = 0.019788;}
  else if (el_pt> 110 && el_pt<= 120 && met> 185 && met<= 180) {eff = 0.92; errup = 0.051501; errdown = 0.0959187;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 185 && met<= 180) {eff = 0.947368; errup = 0.02502; errdown = 0.0396606;}
  else if (el_pt> 20 && el_pt<= 25 && met> 190 && met<= 190) {eff = 0.714286; errup = 0.0643308; errdown = 0.0730111;}
  else if (el_pt> 25 && el_pt<= 30 && met> 190 && met<= 190) {eff = 0.771429; errup = 0.0749545; errdown = 0.0932188;}
  else if (el_pt> 30 && el_pt<= 110 && met> 190 && met<= 190) {eff = 0.897143; errup = 0.0164322; errdown = 0.018897;}
  else if (el_pt> 110 && el_pt<= 120 && met> 190 && met<= 190) {eff = 1; errup = 0; errdown = 0.0738409;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 190 && met<= 190) {eff = 0.984848; errup = 0.0125375; errdown = 0.0339781;}
  else if (el_pt> 20 && el_pt<= 25 && met> 195 && met<= 200) {eff = 0.759259; errup = 0.0613529; errdown = 0.0723448;}
  else if (el_pt> 25 && el_pt<= 30 && met> 195 && met<= 200) {eff = 0.75; errup = 0.0767393; errdown = 0.0929856;}
  else if (el_pt> 30 && el_pt<= 110 && met> 195 && met<= 200) {eff = 0.907534; errup = 0.0171339; errdown = 0.0202008;}
  else if (el_pt> 110 && el_pt<= 120 && met> 195 && met<= 200) {eff = 0.851852; errup = 0.0695101; errdown = 0.101731;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 195 && met<= 200) {eff = 0.969231; errup = 0.0198493; errdown = 0.0391458;}
  else if (el_pt> 20 && el_pt<= 25 && met> 200 && met<= 210) {eff = 0.857143; errup = 0.0551652; errdown = 0.0755709;}
  else if (el_pt> 25 && el_pt<= 30 && met> 200 && met<= 210) {eff = 0.916667; errup = 0.0450072; errdown = 0.0744674;}
  else if (el_pt> 30 && el_pt<= 110 && met> 200 && met<= 210) {eff = 0.94386; errup = 0.0136633; errdown = 0.0171713;}
  else if (el_pt> 110 && el_pt<= 120 && met> 200 && met<= 210) {eff = 0.9375; errup = 0.051761; errdown = 0.129429;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 200 && met<= 210) {eff = 1; errup = 0; errdown = 0.0292574;}
  else if (el_pt> 20 && el_pt<= 25 && met> 210 && met<= 9999) {eff = 0.918919; errup = 0.026025; errdown = 0.034831;}
  else if (el_pt> 25 && el_pt<= 30 && met> 210 && met<= 9999) {eff = 0.865169; errup = 0.0369968; errdown = 0.0464051;}
  else if (el_pt> 30 && el_pt<= 110 && met> 210 && met<= 9999) {eff = 0.965271; errup = 0.00585514; errdown = 0.0068881;}
  else if (el_pt> 110 && el_pt<= 120 && met> 210 && met<= 9999) {eff = 0.989796; errup = 0.00844284; errdown = 0.0230719;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 210 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00822172;}
  return eff;
});

const NamedFunc get_1mu_trigeff2016("get_1mu_trigeff2016", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = b.mu_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 110) {eff = 0.115848; errup = 0.00702761; errdown = 0.00668609;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 110) {eff = 0.50953; errup = 0.0129078; errdown = 0.0129202;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 110) {eff = 0.646394; errup = 0.00726734; errdown = 0.00733457;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 110) {eff = 0.893396; errup = 0.00388525; errdown = 0.00401146;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 155 && met<= 120) {eff = 0.343284; errup = 0.0458476; errdown = 0.043308;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 155 && met<= 120) {eff = 0.709091; errup = 0.0455501; errdown = 0.049726;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 155 && met<= 120) {eff = 0.771331; errup = 0.0252171; errdown = 0.0271938;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 155 && met<= 120) {eff = 0.938636; errup = 0.0114961; errdown = 0.0136851;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 160 && met<= 130) {eff = 0.408333; errup = 0.0496489; errdown = 0.0479882;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 160 && met<= 130) {eff = 0.757895; errup = 0.0459589; errdown = 0.0519967;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 160 && met<= 130) {eff = 0.825758; errup = 0.0238784; errdown = 0.0265402;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 160 && met<= 130) {eff = 0.930946; errup = 0.0129022; errdown = 0.0153227;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 165 && met<= 140) {eff = 0.406977; errup = 0.0596105; errdown = 0.0572279;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 165 && met<= 140) {eff = 0.80303; errup = 0.0509296; errdown = 0.0614258;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 165 && met<= 140) {eff = 0.88; errup = 0.0220141; errdown = 0.0257281;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 165 && met<= 140) {eff = 0.950685; errup = 0.011346; errdown = 0.0141095;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 170 && met<= 150) {eff = 0.559524; errup = 0.0589468; errdown = 0.0605076;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 170 && met<= 150) {eff = 0.814286; errup = 0.0481926; errdown = 0.0584532;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 170 && met<= 150) {eff = 0.873786; errup = 0.0235445; errdown = 0.027542;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 170 && met<= 150) {eff = 0.951429; errup = 0.0114981; errdown = 0.0143934;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 160) {eff = 0.597015; errup = 0.0652245; errdown = 0.0684491;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 160) {eff = 0.90566; errup = 0.0401398; errdown = 0.0587823;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 160) {eff = 0.867299; errup = 0.0237886; errdown = 0.027612;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 160) {eff = 0.968254; errup = 0.00979478; errdown = 0.0132384;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 180 && met<= 170) {eff = 0.716418; errup = 0.0584326; errdown = 0.065696;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 180 && met<= 170) {eff = 0.93617; errup = 0.0345373; errdown = 0.0582116;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 180 && met<= 170) {eff = 0.955128; errup = 0.0164018; errdown = 0.0233331;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 180 && met<= 170) {eff = 0.9699; errup = 0.00977778; errdown = 0.0134446;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 185 && met<= 180) {eff = 0.784615; errup = 0.0532994; errdown = 0.0632744;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 185 && met<= 180) {eff = 0.872727; errup = 0.0457212; errdown = 0.061809;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 185 && met<= 180) {eff = 0.951515; errup = 0.016605; errdown = 0.0230425;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 185 && met<= 180) {eff = 0.977358; errup = 0.00894499; errdown = 0.0132786;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 190 && met<= 190) {eff = 0.733333; errup = 0.0702056; errdown = 0.082142;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 190 && met<= 190) {eff = 0.983051; errup = 0.0140254; errdown = 0.037896;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 190 && met<= 190) {eff = 0.94702; errup = 0.0181267; errdown = 0.0250908;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 190 && met<= 190) {eff = 0.987069; errup = 0.00702944; errdown = 0.0124182;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 195 && met<= 200) {eff = 0.863636; errup = 0.0527265; errdown = 0.0725604;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 195 && met<= 200) {eff = 1; errup = 0; errdown = 0.0449824;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 195 && met<= 200) {eff = 0.934211; errup = 0.0201296; errdown = 0.0267482;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 195 && met<= 200) {eff = 0.988142; errup = 0.00644659; errdown = 0.0113996;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 200 && met<= 210) {eff = 0.853659; errup = 0.0564709; errdown = 0.0771695;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 200 && met<= 210) {eff = 0.916667; errup = 0.0450072; errdown = 0.0744674;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 200 && met<= 210) {eff = 0.974359; errup = 0.0139228; errdown = 0.024313;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 200 && met<= 210) {eff = 0.979381; errup = 0.00984218; errdown = 0.0160023;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 210 && met<= 9999) {eff = 0.92562; errup = 0.0239126; errdown = 0.0321199;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 210 && met<= 9999) {eff = 0.955752; errup = 0.018981; errdown = 0.0288242;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 210 && met<= 9999) {eff = 0.976959; errup = 0.00712399; errdown = 0.00966982;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 210 && met<= 9999) {eff = 0.993455; errup = 0.00282399; errdown = 0.00440304;}
  return eff;
});

const NamedFunc get_2el_trigeff2016("get_2el_trigeff2016", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = b.el_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.819209; errup = 0.0289267; errdown = 0.0289267;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.887324; errup = 0.0216654; errdown = 0.0216654;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.865741; errup = 0.0231974; errdown = 0.0231974;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.880597; errup = 0.0198075; errdown = 0.0198075;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.906383; errup = 0.019002; errdown = 0.019002;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.851351; errup = 0.0238758; errdown = 0.0238758;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.885; errup = 0.0225583; errdown = 0.0225583;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.920455; errup = 0.0203964; errdown = 0.0203964;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.917197; errup = 0.021994; errdown = 0.021994;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.902439; errup = 0.0267544; errdown = 0.0267544;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.936508; errup = 0.0217235; errdown = 0.0217235;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.936508; errup = 0.0217235; errdown = 0.0217235;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.895238; errup = 0.0298866; errdown = 0.0298866;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.93617; errup = 0.025213; errdown = 0.025213;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.988395; errup = 0.00210652; errdown = 0.00210652;}
  return eff;
});

const NamedFunc get_2mu_trigeff2016("get_2mu_trigeff2016", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = b.mu_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.933594; errup = 0.0155619; errdown = 0.0155619;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.968288; errup = 0.00805725; errdown = 0.00805725;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.972572; errup = 0.000505438; errdown = 0.000505438;}
  return eff;
});

}
