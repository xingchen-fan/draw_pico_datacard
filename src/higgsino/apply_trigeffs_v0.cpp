#include <vector>
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/apply_trigeffs_v0.hpp"

namespace Higfuncs{

const NamedFunc get_0l_trigeff2016_v0("get_0l_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_0l_fakemet_trigeff2016_v0("get_0l_fakemet_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.34375; errup = 0.0551838; errdown = 0.0516032;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.329126; errup = 0.00748338; errdown = 0.00739907;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.340784; errup = 0.00565815; errdown = 0.00561368;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.318462; errup = 0.00311091; errdown = 0.00309492;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.267794; errup = 0.00247084; errdown = 0.00245656;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.30303; errup = 0.066999; errdown = 0.0602987;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.361364; errup = 0.00827805; errdown = 0.00819792;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.377637; errup = 0.00635196; errdown = 0.00631077;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.359295; errup = 0.00349402; errdown = 0.00347927;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.30763; errup = 0.00282814; errdown = 0.00281386;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.28; errup = 0.0777781; errdown = 0.067729;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.4; errup = 0.0093527; errdown = 0.00928171;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.418022; errup = 0.00701682; errdown = 0.00698425;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.407438; errup = 0.00392148; errdown = 0.00390982;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.346154; errup = 0.00319023; errdown = 0.00317655;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.186047; errup = 0.0788698; errdown = 0.0616867;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.458898; errup = 0.0104849; errdown = 0.0104493;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.456691; errup = 0.00778; errdown = 0.00775924;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.436683; errup = 0.00430337; errdown = 0.00429393;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.379798; errup = 0.00353981; errdown = 0.00352716;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.277778; errup = 0.0944452; errdown = 0.0800801;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.47762; errup = 0.0115283; errdown = 0.011505;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.496503; errup = 0.00850186; errdown = 0.00849987;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.470105; errup = 0.00468488; errdown = 0.00467967;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.421681; errup = 0.00393157; errdown = 0.00392175;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.324324; errup = 0.0946276; errdown = 0.0836685;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.508585; errup = 0.0124515; errdown = 0.0124619;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.532835; errup = 0.00914594; errdown = 0.00916767;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.511022; errup = 0.00506406; errdown = 0.0050663;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.45314; errup = 0.00425117; errdown = 0.0042444;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.423077; errup = 0.116981; errdown = 0.110076;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.535191; errup = 0.0138396; errdown = 0.0138925;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.562858; errup = 0.00986103; errdown = 0.00991;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.54331; errup = 0.0054092; errdown = 0.00541933;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.490859; errup = 0.00467494; errdown = 0.00467336;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.578947; errup = 0.1303; errdown = 0.140168;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.557904; errup = 0.0148052; errdown = 0.0149059;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.595641; errup = 0.0105275; errdown = 0.0106144;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.57839; errup = 0.00581006; errdown = 0.00583159;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.526404; errup = 0.00499295; errdown = 0.00499819;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.416667; errup = 0.183083; errdown = 0.166178;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.543564; errup = 0.0161155; errdown = 0.0162044;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.616804; errup = 0.0108963; errdown = 0.0110123;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.605413; errup = 0.0062108; errdown = 0.00624459;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.566151; errup = 0.00533269; errdown = 0.00534789;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.6; errup = 0.145412; errdown = 0.161475;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.579744; errup = 0.0173172; errdown = 0.0175091;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.632; errup = 0.0117334; errdown = 0.0118879;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.628879; errup = 0.00647922; errdown = 0.00652528;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.593062; errup = 0.0056954; errdown = 0.00572023;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.666667; errup = 0.135514; errdown = 0.162494;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.607942; errup = 0.0133945; errdown = 0.0135546;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.646354; errup = 0.00896477; errdown = 0.0090669;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.666143; errup = 0.00486205; errdown = 0.00489725;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.627819; errup = 0.00436025; errdown = 0.00438096;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.684211; errup = 0.117339; errdown = 0.140625;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.659864; errup = 0.0150872; errdown = 0.015408;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.66796; errup = 0.00991916; errdown = 0.0100672;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.711698; errup = 0.00522735; errdown = 0.00528364;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.673255; errup = 0.00471847; errdown = 0.00475342;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.375; errup = 0.156049; errdown = 0.137247;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.666667; errup = 0.0172876; errdown = 0.0177309;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.687773; errup = 0.0109918; errdown = 0.0112019;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.734524; errup = 0.00580289; errdown = 0.00588383;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.704943; errup = 0.00519146; errdown = 0.00524446;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.454545; errup = 0.189662; errdown = 0.179582;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.709898; errup = 0.0192106; errdown = 0.0199588;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.692726; errup = 0.0122306; errdown = 0.0124999;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.760484; errup = 0.00621321; errdown = 0.00632373;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.754891; errup = 0.00546306; errdown = 0.00554532;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.9; errup = 0.082873; errdown = 0.194135;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.696121; errup = 0.0219634; errdown = 0.0228502;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.71073; errup = 0.0135217; errdown = 0.0138949;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.79273; errup = 0.00655972; errdown = 0.00671362;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.774911; errup = 0.00590379; errdown = 0.00601387;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.777778; errup = 0.142118; errdown = 0.221429;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.71466; errup = 0.01669; errdown = 0.0172739;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.709596; errup = 0.0103435; errdown = 0.0105604;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.807188; errup = 0.0047099; errdown = 0.00479779;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.815172; errup = 0.00401727; errdown = 0.00408503;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0; errup = 0.308024; errdown = 0;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.694118; errup = 0.0230215; errdown = 0.0239817;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.722388; errup = 0.0124316; errdown = 0.0127735;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.845351; errup = 0.00547786; errdown = 0.00563772;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.851237; errup = 0.00456296; errdown = 0.00467933;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 350) {eff = 0.777778; errup = 0.142118; errdown = 0.221429;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.735802; errup = 0.0225042; errdown = 0.0237324;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.716895; errup = 0.0126343; errdown = 0.0129743;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.862974; errup = 0.00495149; errdown = 0.00510348;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.874064; errup = 0.00403698; errdown = 0.00414878;}
  else if (ht> 0 && ht<= 200 && met> 350 && met<= 400) {eff = 1; errup = 0; errdown = 0.458642;}
  else if (ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.709677; errup = 0.0380961; errdown = 0.0410288;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.748673; errup = 0.0186589; errdown = 0.0195808;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.889896; errup = 0.00670688; errdown = 0.00707213;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.908459; errup = 0.00503657; errdown = 0.00529077;}
  else if (ht> 0 && ht<= 200 && met> 400 && met<= 450) {eff = 1; errup = 0; errdown = 0.841345;}
  else if (ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.78481; errup = 0.0482156; errdown = 0.0563445;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.792531; errup = 0.0268436; errdown = 0.0294573;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.909605; errup = 0.00886479; errdown = 0.00967703;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.933579; errup = 0.00620482; errdown = 0.00676451;}
  else if (ht> 0 && ht<= 200 && met> 450 && met<= 500) {eff = 0.5; errup = 0.417248; errdown = 0.417248;}
  else if (ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.772727; errup = 0.0944237; errdown = 0.12452;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.87037; errup = 0.0329774; errdown = 0.0407707;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.94084; errup = 0.0103533; errdown = 0.0121897;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.941714; errup = 0.00795617; errdown = 0.009038;}
  else if (ht> 0 && ht<= 200 && met> 500 && met<= 9999) {eff = 9.53547e-322; errup = 8.74496e-322; errdown = 8.74496e-322;}
  else if (ht> 200 && ht<= 600 && met> 500 && met<= 9999) {eff = 0.714286; errup = 0.182129; errdown = 0.259938;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.971014; errup = 0.0187; errdown = 0.0369543;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.961165; errup = 0.0109425; errdown = 0.0143599;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.948919; errup = 0.00978575; errdown = 0.0117251;}
  return eff;
});

const NamedFunc get_1el_trigeff2016_v0("get_1el_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 20 && el_pt<= 25 && met> 0 && met<= 110) {eff = 0.0582445; errup = 0.00513861; errdown = 0.00476062;}
  else if (el_pt> 25 && el_pt<= 30 && met> 0 && met<= 110) {eff = 0.157254; errup = 0.00903781; errdown = 0.00864611;}
  else if (el_pt> 30 && el_pt<= 110 && met> 0 && met<= 110) {eff = 0.530553; errup = 0.00429568; errdown = 0.00430018;}
  else if (el_pt> 110 && el_pt<= 120 && met> 0 && met<= 110) {eff = 0.650194; errup = 0.0151726; errdown = 0.015473;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 110) {eff = 0.89841; errup = 0.0059834; errdown = 0.00630236;}
  else if (el_pt> 20 && el_pt<= 25 && met> 110 && met<= 120) {eff = 0.149123; errup = 0.0411162; errdown = 0.0341911;}
  else if (el_pt> 25 && el_pt<= 30 && met> 110 && met<= 120) {eff = 0.338028; errup = 0.0652582; errdown = 0.0601709;}
  else if (el_pt> 30 && el_pt<= 110 && met> 110 && met<= 120) {eff = 0.672389; errup = 0.0181985; errdown = 0.0187112;}
  else if (el_pt> 110 && el_pt<= 120 && met> 110 && met<= 120) {eff = 0.686275; errup = 0.06986; errdown = 0.078157;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 110 && met<= 120) {eff = 0.969231; errup = 0.0146686; errdown = 0.02366;}
  else if (el_pt> 20 && el_pt<= 25 && met> 120 && met<= 130) {eff = 0.2; errup = 0.049646; errdown = 0.0425693;}
  else if (el_pt> 25 && el_pt<= 30 && met> 120 && met<= 130) {eff = 0.314286; errup = 0.0651353; errdown = 0.0592025;}
  else if (el_pt> 30 && el_pt<= 110 && met> 120 && met<= 130) {eff = 0.684039; errup = 0.0192486; errdown = 0.0198735;}
  else if (el_pt> 110 && el_pt<= 120 && met> 120 && met<= 130) {eff = 0.886364; errup = 0.0481931; errdown = 0.0695663;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 120 && met<= 130) {eff = 0.95082; errup = 0.0193366; errdown = 0.02822;}
  else if (el_pt> 20 && el_pt<= 25 && met> 130 && met<= 140) {eff = 0.425287; errup = 0.0593282; errdown = 0.0574387;}
  else if (el_pt> 25 && el_pt<= 30 && met> 130 && met<= 140) {eff = 0.513514; errup = 0.0641717; errdown = 0.0645758;}
  else if (el_pt> 30 && el_pt<= 110 && met> 130 && met<= 140) {eff = 0.723894; errup = 0.0192634; errdown = 0.0200923;}
  else if (el_pt> 110 && el_pt<= 120 && met> 130 && met<= 140) {eff = 0.909091; errup = 0.0429965; errdown = 0.0660622;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 130 && met<= 140) {eff = 0.892157; errup = 0.0311385; errdown = 0.0399966;}
  else if (el_pt> 20 && el_pt<= 25 && met> 140 && met<= 150) {eff = 0.453333; errup = 0.0643662; errdown = 0.0629896;}
  else if (el_pt> 25 && el_pt<= 30 && met> 140 && met<= 150) {eff = 0.464286; errup = 0.0755319; errdown = 0.0741026;}
  else if (el_pt> 30 && el_pt<= 110 && met> 140 && met<= 150) {eff = 0.765531; errup = 0.0193923; errdown = 0.0205114;}
  else if (el_pt> 110 && el_pt<= 120 && met> 140 && met<= 150) {eff = 0.844444; errup = 0.0555301; errdown = 0.0737002;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 140 && met<= 150) {eff = 0.947917; errup = 0.0223143; errdown = 0.0336975;}
  else if (el_pt> 20 && el_pt<= 25 && met> 150 && met<= 160) {eff = 0.606557; errup = 0.0681557; errdown = 0.0720644;}
  else if (el_pt> 25 && el_pt<= 30 && met> 150 && met<= 160) {eff = 0.731707; errup = 0.0738471; errdown = 0.0869218;}
  else if (el_pt> 30 && el_pt<= 110 && met> 150 && met<= 160) {eff = 0.816934; errup = 0.0188601; errdown = 0.0203981;}
  else if (el_pt> 110 && el_pt<= 120 && met> 150 && met<= 160) {eff = 0.765957; errup = 0.0651312; errdown = 0.0782067;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 160) {eff = 0.987013; errup = 0.010746; errdown = 0.029229;}
  else if (el_pt> 20 && el_pt<= 25 && met> 160 && met<= 170) {eff = 0.555556; errup = 0.0827174; errdown = 0.0855148;}
  else if (el_pt> 25 && el_pt<= 30 && met> 160 && met<= 170) {eff = 0.75; errup = 0.0815095; errdown = 0.0999181;}
  else if (el_pt> 30 && el_pt<= 110 && met> 160 && met<= 170) {eff = 0.836364; errup = 0.0192059; errdown = 0.0210707;}
  else if (el_pt> 110 && el_pt<= 120 && met> 160 && met<= 170) {eff = 0.965517; errup = 0.0285434; errdown = 0.0748731;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 160 && met<= 170) {eff = 1; errup = 0; errdown = 0.0287997;}
  else if (el_pt> 20 && el_pt<= 25 && met> 170 && met<= 180) {eff = 0.625; errup = 0.0842352; errdown = 0.0913836;}
  else if (el_pt> 25 && el_pt<= 30 && met> 170 && met<= 180) {eff = 0.729167; errup = 0.0683073; errdown = 0.0792509;}
  else if (el_pt> 30 && el_pt<= 110 && met> 170 && met<= 180) {eff = 0.869333; errup = 0.0176707; errdown = 0.019788;}
  else if (el_pt> 110 && el_pt<= 120 && met> 170 && met<= 180) {eff = 0.92; errup = 0.051501; errdown = 0.0959187;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 170 && met<= 180) {eff = 0.947368; errup = 0.02502; errdown = 0.0396606;}
  else if (el_pt> 20 && el_pt<= 25 && met> 180 && met<= 190) {eff = 0.714286; errup = 0.0643308; errdown = 0.0730111;}
  else if (el_pt> 25 && el_pt<= 30 && met> 180 && met<= 190) {eff = 0.771429; errup = 0.0749545; errdown = 0.0932188;}
  else if (el_pt> 30 && el_pt<= 110 && met> 180 && met<= 190) {eff = 0.897143; errup = 0.0164322; errdown = 0.018897;}
  else if (el_pt> 110 && el_pt<= 120 && met> 180 && met<= 190) {eff = 1; errup = 0; errdown = 0.0738409;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 180 && met<= 190) {eff = 0.984848; errup = 0.0125375; errdown = 0.0339781;}
  else if (el_pt> 20 && el_pt<= 25 && met> 190 && met<= 200) {eff = 0.759259; errup = 0.0613529; errdown = 0.0723448;}
  else if (el_pt> 25 && el_pt<= 30 && met> 190 && met<= 200) {eff = 0.75; errup = 0.0767393; errdown = 0.0929856;}
  else if (el_pt> 30 && el_pt<= 110 && met> 190 && met<= 200) {eff = 0.907534; errup = 0.0171339; errdown = 0.0202008;}
  else if (el_pt> 110 && el_pt<= 120 && met> 190 && met<= 200) {eff = 0.851852; errup = 0.0695101; errdown = 0.101731;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 190 && met<= 200) {eff = 0.969231; errup = 0.0198493; errdown = 0.0391458;}
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

const NamedFunc get_1mu_trigeff2016_v0("get_1mu_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 110) {eff = 0.115848; errup = 0.00702761; errdown = 0.00668609;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 110) {eff = 0.50953; errup = 0.0129078; errdown = 0.0129202;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 110) {eff = 0.646394; errup = 0.00726734; errdown = 0.00733457;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 110) {eff = 0.893396; errup = 0.00388525; errdown = 0.00401146;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 110 && met<= 120) {eff = 0.343284; errup = 0.0458476; errdown = 0.043308;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 110 && met<= 120) {eff = 0.709091; errup = 0.0455501; errdown = 0.049726;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 110 && met<= 120) {eff = 0.771331; errup = 0.0252171; errdown = 0.0271938;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 110 && met<= 120) {eff = 0.938636; errup = 0.0114961; errdown = 0.0136851;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 120 && met<= 130) {eff = 0.408333; errup = 0.0496489; errdown = 0.0479882;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 120 && met<= 130) {eff = 0.757895; errup = 0.0459589; errdown = 0.0519967;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 120 && met<= 130) {eff = 0.825758; errup = 0.0238784; errdown = 0.0265402;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 120 && met<= 130) {eff = 0.930946; errup = 0.0129022; errdown = 0.0153227;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 130 && met<= 140) {eff = 0.406977; errup = 0.0596105; errdown = 0.0572279;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 130 && met<= 140) {eff = 0.80303; errup = 0.0509296; errdown = 0.0614258;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 130 && met<= 140) {eff = 0.88; errup = 0.0220141; errdown = 0.0257281;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 130 && met<= 140) {eff = 0.950685; errup = 0.011346; errdown = 0.0141095;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 140 && met<= 150) {eff = 0.559524; errup = 0.0589468; errdown = 0.0605076;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 140 && met<= 150) {eff = 0.814286; errup = 0.0481926; errdown = 0.0584532;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 140 && met<= 150) {eff = 0.873786; errup = 0.0235445; errdown = 0.027542;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 140 && met<= 150) {eff = 0.951429; errup = 0.0114981; errdown = 0.0143934;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 160) {eff = 0.597015; errup = 0.0652245; errdown = 0.0684491;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 160) {eff = 0.90566; errup = 0.0401398; errdown = 0.0587823;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 160) {eff = 0.867299; errup = 0.0237886; errdown = 0.027612;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 160) {eff = 0.968254; errup = 0.00979478; errdown = 0.0132384;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 160 && met<= 170) {eff = 0.716418; errup = 0.0584326; errdown = 0.065696;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 160 && met<= 170) {eff = 0.93617; errup = 0.0345373; errdown = 0.0582116;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 160 && met<= 170) {eff = 0.955128; errup = 0.0164018; errdown = 0.0233331;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 160 && met<= 170) {eff = 0.9699; errup = 0.00977778; errdown = 0.0134446;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 170 && met<= 180) {eff = 0.784615; errup = 0.0532994; errdown = 0.0632744;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 170 && met<= 180) {eff = 0.872727; errup = 0.0457212; errdown = 0.061809;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 170 && met<= 180) {eff = 0.951515; errup = 0.016605; errdown = 0.0230425;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 170 && met<= 180) {eff = 0.977358; errup = 0.00894499; errdown = 0.0132786;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 180 && met<= 190) {eff = 0.733333; errup = 0.0702056; errdown = 0.082142;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 180 && met<= 190) {eff = 0.983051; errup = 0.0140254; errdown = 0.037896;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 180 && met<= 190) {eff = 0.94702; errup = 0.0181267; errdown = 0.0250908;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 180 && met<= 190) {eff = 0.987069; errup = 0.00702944; errdown = 0.0124182;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 190 && met<= 200) {eff = 0.863636; errup = 0.0527265; errdown = 0.0725604;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 190 && met<= 200) {eff = 1; errup = 0; errdown = 0.0449824;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 190 && met<= 200) {eff = 0.934211; errup = 0.0201296; errdown = 0.0267482;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 190 && met<= 200) {eff = 0.988142; errup = 0.00644659; errdown = 0.0113996;}
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

const NamedFunc get_2el_trigeff2016_v0("get_2el_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
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

const NamedFunc get_2mu_trigeff2016_v0("get_2mu_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.933594; errup = 0.0155619; errdown = 0.0155619;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.968288; errup = 0.00805725; errdown = 0.00805725;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.972572; errup = 0.000505438; errdown = 0.000505438;}
  return eff;
});

const NamedFunc get_0l_trigeff2017_v0("get_0l_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.118202; errup = 0.00647114; errdown = 0.00618617;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.228561; errup = 0.00414684; errdown = 0.00409514;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.362438; errup = 0.0205543; errdown = 0.0200824;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.351759; errup = 0.0370682; errdown = 0.0354728;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.493333; errup = 0.0640194; errdown = 0.0638228;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.140786; errup = 0.00792831; errdown = 0.00758193;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.283364; errup = 0.0047502; errdown = 0.00470307;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.395731; errup = 0.0207976; errdown = 0.0204413;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.42328; errup = 0.038911; errdown = 0.0380424;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.493151; errup = 0.0649718; errdown = 0.0647641;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.166766; errup = 0.00959853; errdown = 0.00918874;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.324521; errup = 0.00523779; errdown = 0.00519473;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.458; errup = 0.0233451; errdown = 0.0231698;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.477157; errup = 0.038145; errdown = 0.0378973;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.38961; errup = 0.0632202; errdown = 0.0600452;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.205904; errup = 0.0115753; errdown = 0.0111268;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.390841; errup = 0.00582333; errdown = 0.00579277;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.5; errup = 0.0234834; errdown = 0.0234834;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.535032; errup = 0.042604; errdown = 0.043084;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.588235; errup = 0.0650517; errdown = 0.0679384;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.241224; errup = 0.0135216; errdown = 0.0130397;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.432834; errup = 0.00624755; errdown = 0.00622653;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.581395; errup = 0.0247251; errdown = 0.0251216;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.620915; errup = 0.0414935; errdown = 0.0431987;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.716981; errup = 0.0659706; errdown = 0.0752878;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.26533; errup = 0.0160304; errdown = 0.0154558;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.480098; errup = 0.00670629; errdown = 0.00669922;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.593301; errup = 0.0249575; errdown = 0.0254254;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.649007; errup = 0.0409427; errdown = 0.0430755;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.550725; errup = 0.0657743; errdown = 0.0674068;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.341679; errup = 0.0191662; errdown = 0.0186814;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.539748; errup = 0.00708134; errdown = 0.0070972;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.654731; errup = 0.0248756; errdown = 0.0257078;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.71875; errup = 0.0416165; errdown = 0.0453516;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.692308; errup = 0.0686825; errdown = 0.0770809;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.369732; errup = 0.0223259; errdown = 0.0218048;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.596124; errup = 0.00745613; errdown = 0.00750012;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.701333; errup = 0.0243702; errdown = 0.0255026;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.691176; errup = 0.0416075; errdown = 0.0446653;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.7; errup = 0.0780739; errdown = 0.0896007;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.383966; errup = 0.0236247; errdown = 0.0231127;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.633548; errup = 0.0078334; errdown = 0.00790347;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.667665; errup = 0.0267017; errdown = 0.027762;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.766129; errup = 0.0395363; errdown = 0.0442593;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.733333; errup = 0.0702056; errdown = 0.082142;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.430267; errup = 0.0286282; errdown = 0.0281924;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.681906; errup = 0.00783883; errdown = 0.00794142;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.775801; errup = 0.0255714; errdown = 0.0276701;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.704918; errup = 0.043372; errdown = 0.0470452;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.777778; errup = 0.0651643; errdown = 0.0795029;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.51153; errup = 0.0238841; errdown = 0.0239346;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.7266; errup = 0.00589238; errdown = 0.0059715;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.818681; errup = 0.0167789; errdown = 0.0180097;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.813084; errup = 0.027352; errdown = 0.0305252;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.818182; errup = 0.0454972; errdown = 0.0549007;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.545763; errup = 0.0304653; errdown = 0.030793;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.797219; errup = 0.00592227; errdown = 0.00605172;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.864198; errup = 0.0157643; errdown = 0.0173603;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.882716; errup = 0.0256997; errdown = 0.0309806;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.861111; errup = 0.041698; errdown = 0.0533274;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.585859; errup = 0.0369679; errdown = 0.0378947;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.841494; errup = 0.0061156; errdown = 0.0063088;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.893665; errup = 0.0148345; errdown = 0.0167507;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.948529; errup = 0.0187887; errdown = 0.0266233;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.892308; errup = 0.0388522; errdown = 0.0531815;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.740157; errup = 0.0406252; errdown = 0.0447687;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.885054; errup = 0.00587368; errdown = 0.00613912;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.929775; errup = 0.0136293; errdown = 0.0162897;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.920863; errup = 0.0230263; errdown = 0.0299998;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.942308; errup = 0.0312343; errdown = 0.0529489;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.7625; errup = 0.0498414; errdown = 0.0572049;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.921169; errup = 0.00569931; errdown = 0.00608683;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.959752; errup = 0.0109001; errdown = 0.0141436;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.972727; errup = 0.0148066; errdown = 0.0258178;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.953488; errup = 0.0299851; errdown = 0.0580755;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.82; errup = 0.0396683; errdown = 0.0468539;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.938905; errup = 0.00390786; errdown = 0.00414787;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.964561; errup = 0.00725305; errdown = 0.00883537;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.969828; errup = 0.0110614; errdown = 0.0158737;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.951456; errup = 0.0208096; errdown = 0.0315047;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.868421; errup = 0.0556306; errdown = 0.0792159;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.965547; errup = 0.00381497; errdown = 0.00424421;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.970954; errup = 0.00760867; errdown = 0.0098252;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.96988; errup = 0.0129497; errdown = 0.0198622;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.986486; errup = 0.0111817; errdown = 0.0303874;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 9999) {eff = 0.878788; errup = 0.0570897; errdown = 0.0855134;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.968519; errup = 0.00376154; errdown = 0.004222;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.968439; errup = 0.00710846; errdown = 0.0088446;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.981308; errup = 0.00892453; errdown = 0.0145322;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 0.928571; errup = 0.0259675; errdown = 0.0363535;}
  return eff;
});

const NamedFunc get_0l_fakemet_trigeff2017_v0("get_0l_fakemet_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.123596; errup = 0.0452794; errdown = 0.0355329;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.188815; errup = 0.00513021; errdown = 0.00502708;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.279651; errup = 0.00605781; errdown = 0.00597971;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.31749; errup = 0.00544657; errdown = 0.00539765;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.37227; errup = 0.00307408; errdown = 0.00306385;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.15; errup = 0.060943; errdown = 0.0473246;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.238992; errup = 0.00625925; errdown = 0.00615056;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.319574; errup = 0.00694655; errdown = 0.00686858;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.350866; errup = 0.00620533; errdown = 0.00615597;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.412794; errup = 0.00341814; errdown = 0.00340981;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.266667; errup = 0.082142; errdown = 0.0702056;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.266804; errup = 0.0072869; errdown = 0.00716463;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.361294; errup = 0.0079572; errdown = 0.00788304;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.378864; errup = 0.00686246; errdown = 0.00681498;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.443561; errup = 0.00374757; errdown = 0.00374121;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.104167; errup = 0.0643268; errdown = 0.044248;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.309732; errup = 0.00870177; errdown = 0.00857157;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.404604; errup = 0.00890297; errdown = 0.00884172;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.43023; errup = 0.00759906; errdown = 0.00756681;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.477556; errup = 0.00411404; errdown = 0.00411102;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.371429; errup = 0.0986578; errdown = 0.0901926;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.352278; errup = 0.0102298; errdown = 0.0100988;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.423592; errup = 0.00988964; errdown = 0.00982997;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.463069; errup = 0.00837741; errdown = 0.00835695;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.507831; errup = 0.00445136; errdown = 0.00445259;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.166667; errup = 0.112297; errdown = 0.0779889;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.38976; errup = 0.0114647; errdown = 0.0113469;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.486221; errup = 0.0110393; errdown = 0.0110261;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.49855; errup = 0.009136; errdown = 0.00913505;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.546107; errup = 0.00487945; errdown = 0.00488823;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.277778; errup = 0.144222; errdown = 0.114258;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.425186; errup = 0.0132496; errdown = 0.0131457;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.496301; errup = 0.0122113; errdown = 0.012207;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.540927; errup = 0.00996673; errdown = 0.00999895;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.574945; errup = 0.00522855; errdown = 0.00524519;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.235294; errup = 0.147312; errdown = 0.108959;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.493966; errup = 0.0145821; errdown = 0.0145721;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.561616; errup = 0.0131657; errdown = 0.0132508;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.575424; errup = 0.0103516; errdown = 0.0104168;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.623677; errup = 0.0055526; errdown = 0.0055849;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.272727; errup = 0.196072; errdown = 0.144396;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.507331; errup = 0.0161042; errdown = 0.016119;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.597656; errup = 0.0140129; errdown = 0.0141698;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.579246; errup = 0.0114304; errdown = 0.0115141;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.642391; errup = 0.0059573; errdown = 0.00600106;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.25; errup = 0.239567; errdown = 0.159659;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.539424; errup = 0.0181966; errdown = 0.0182987;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.653775; errup = 0.0147448; errdown = 0.0150369;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.623623; errup = 0.01221; errdown = 0.0123651;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.67077; errup = 0.00624379; errdown = 0.00630382;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.571429; errup = 0.222488; errdown = 0.247841;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.601674; errup = 0.0144865; errdown = 0.0146616;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.664822; errup = 0.0119066; errdown = 0.0121146;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.64353; errup = 0.00906688; errdown = 0.00916893;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.694319; errup = 0.00476164; errdown = 0.00480306;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.4; errup = 0.303366; errdown = 0.253348;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.661078; errup = 0.0167699; errdown = 0.0171696;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.753149; errup = 0.0126899; errdown = 0.013129;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.693149; errup = 0.0102007; errdown = 0.0103888;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.747028; errup = 0.0051753; errdown = 0.00524533;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.2; errup = 0.324251; errdown = 0.166039;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.721501; errup = 0.0174053; errdown = 0.018071;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.786091; errup = 0.0135196; errdown = 0.014146;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.708743; errup = 0.0110955; errdown = 0.0113436;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.778038; errup = 0.00561447; errdown = 0.00571617;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.5; errup = 0.417248; errdown = 0.417248;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.781053; errup = 0.019387; errdown = 0.0206344;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.792402; errup = 0.0151964; errdown = 0.0160253;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.739066; errup = 0.0121406; errdown = 0.0125059;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.806997; errup = 0.00609081; errdown = 0.00623776;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.75; errup = 0.207731; errdown = 0.368402;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.809249; errup = 0.0215917; errdown = 0.0234981;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.834813; errup = 0.0159028; errdown = 0.0171583;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.760676; errup = 0.0129268; errdown = 0.0134066;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.849689; errup = 0.00619514; errdown = 0.00640745;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.444444; errup = 0.213463; errdown = 0.198267;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.847122; errup = 0.0154938; errdown = 0.0168146;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.865031; errup = 0.011045; errdown = 0.0118259;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.772081; errup = 0.0095614; errdown = 0.00984496;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.854076; errup = 0.0045077; errdown = 0.00462405;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0; errup = 0.458642; errdown = 0;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.892857; errup = 0.0178531; errdown = 0.0206337;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.879159; errup = 0.0138072; errdown = 0.015219;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.799681; errup = 0.0114553; errdown = 0.0119503;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.892884; errup = 0.0049459; errdown = 0.00514997;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 350) {eff = 1; errup = 0; errdown = 0.841345;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.929515; errup = 0.0170808; errdown = 0.0213439;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.873702; errup = 0.0139923; errdown = 0.015365;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.773749; errup = 0.0114217; errdown = 0.0118313;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.911984; errup = 0.0042738; errdown = 0.00446466;}
  else if (ht> 0 && ht<= 200 && met> 350 && met<= 400) {eff = 1.80842e-308; errup = 9.35748e-307; errdown = 1.69349e-306;}
  else if (ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.923077; errup = 0.0279331; errdown = 0.038974;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.853061; errup = 0.0230683; errdown = 0.0262025;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.790164; errup = 0.0167972; errdown = 0.0177949;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.929591; errup = 0.00560601; errdown = 0.00603208;}
  else if (ht> 0 && ht<= 200 && met> 400 && met<= 450) {eff = 9.53547e-322; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.947368; errup = 0.0339218; errdown = 0.0652364;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.834862; errup = 0.0365832; errdown = 0.0434785;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.798535; errup = 0.0248984; errdown = 0.0272452;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.945255; errup = 0.00689888; errdown = 0.00776494;}
  else if (ht> 0 && ht<= 200 && met> 450 && met<= 500) {eff = 1; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.818182; errup = 0.116508; errdown = 0.191402;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.888889; errup = 0.0400561; errdown = 0.0547111;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.780488; errup = 0.0334179; errdown = 0.0371466;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.949275; errup = 0.00936678; errdown = 0.0111507;}
  else if (ht> 0 && ht<= 200 && met> 500 && met<= 9999) {eff = 1; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 500 && met<= 9999) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.807692; errup = 0.0804259; errdown = 0.109172;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.8625; errup = 0.0393773; errdown = 0.0498223;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.957377; errup = 0.0115359; errdown = 0.0149532;}
  return eff;
});

const NamedFunc get_1el_trigeff2017_v0("get_1el_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 20 && el_pt<= 25 && met> 0 && met<= 110) {eff = 0.0619935; errup = 0.00524257; errdown = 0.00487288;}
  else if (el_pt> 25 && el_pt<= 30 && met> 0 && met<= 110) {eff = 0.118907; errup = 0.00796992; errdown = 0.00754721;}
  else if (el_pt> 30 && el_pt<= 110 && met> 0 && met<= 110) {eff = 0.552149; errup = 0.00415803; errdown = 0.00416527;}
  else if (el_pt> 110 && el_pt<= 120 && met> 0 && met<= 110) {eff = 0.585569; errup = 0.0153569; errdown = 0.01552;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 110) {eff = 0.82829; errup = 0.0072596; errdown = 0.00750497;}
  else if (el_pt> 20 && el_pt<= 25 && met> 110 && met<= 120) {eff = 0.213793; errup = 0.0395552; errdown = 0.0352269;}
  else if (el_pt> 25 && el_pt<= 30 && met> 110 && met<= 120) {eff = 0.260504; errup = 0.0464748; errdown = 0.0420534;}
  else if (el_pt> 30 && el_pt<= 110 && met> 110 && met<= 120) {eff = 0.676737; errup = 0.0151551; errdown = 0.0155229;}
  else if (el_pt> 110 && el_pt<= 120 && met> 110 && met<= 120) {eff = 0.701493; errup = 0.0595131; errdown = 0.0662634;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 110 && met<= 120) {eff = 0.917647; errup = 0.0212488; errdown = 0.0268265;}
  else if (el_pt> 20 && el_pt<= 25 && met> 120 && met<= 130) {eff = 0.364407; errup = 0.0495787; errdown = 0.0470745;}
  else if (el_pt> 25 && el_pt<= 30 && met> 120 && met<= 130) {eff = 0.350515; errup = 0.0550056; errdown = 0.0516183;}
  else if (el_pt> 30 && el_pt<= 110 && met> 120 && met<= 130) {eff = 0.738876; errup = 0.0153181; errdown = 0.0158989;}
  else if (el_pt> 110 && el_pt<= 120 && met> 120 && met<= 130) {eff = 0.746269; errup = 0.056044; errdown = 0.0643443;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 120 && met<= 130) {eff = 0.913907; errup = 0.0230206; errdown = 0.0292709;}
  else if (el_pt> 20 && el_pt<= 25 && met> 130 && met<= 140) {eff = 0.359551; errup = 0.0578891; errdown = 0.0544101;}
  else if (el_pt> 25 && el_pt<= 30 && met> 130 && met<= 140) {eff = 0.464286; errup = 0.0604404; errdown = 0.0595044;}
  else if (el_pt> 30 && el_pt<= 110 && met> 130 && met<= 140) {eff = 0.734637; errup = 0.0168473; errdown = 0.0175299;}
  else if (el_pt> 110 && el_pt<= 120 && met> 130 && met<= 140) {eff = 0.710145; errup = 0.0579991; errdown = 0.0648324;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 130 && met<= 140) {eff = 0.93985; errup = 0.0205473; errdown = 0.0283269;}
  else if (el_pt> 20 && el_pt<= 25 && met> 140 && met<= 150) {eff = 0.639535; errup = 0.0554437; errdown = 0.0590257;}
  else if (el_pt> 25 && el_pt<= 30 && met> 140 && met<= 150) {eff = 0.6; errup = 0.0661652; errdown = 0.0695965;}
  else if (el_pt> 30 && el_pt<= 110 && met> 140 && met<= 150) {eff = 0.78254; errup = 0.0167488; errdown = 0.0176882;}
  else if (el_pt> 110 && el_pt<= 120 && met> 140 && met<= 150) {eff = 0.822222; errup = 0.0590681; errdown = 0.0759182;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 140 && met<= 150) {eff = 0.919643; errup = 0.0257972; errdown = 0.0345396;}
  else if (el_pt> 20 && el_pt<= 25 && met> 150 && met<= 160) {eff = 0.573529; errup = 0.0655803; errdown = 0.0679846;}
  else if (el_pt> 25 && el_pt<= 30 && met> 150 && met<= 160) {eff = 0.563636; errup = 0.0738937; errdown = 0.0764906;}
  else if (el_pt> 30 && el_pt<= 110 && met> 150 && met<= 160) {eff = 0.82459; errup = 0.0156505; errdown = 0.0167703;}
  else if (el_pt> 110 && el_pt<= 120 && met> 150 && met<= 160) {eff = 0.893617; errup = 0.0451725; errdown = 0.0655621;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 160) {eff = 0.963303; errup = 0.0174813; errdown = 0.028065;}
  else if (el_pt> 20 && el_pt<= 25 && met> 160 && met<= 170) {eff = 0.607143; errup = 0.0713147; errdown = 0.0756132;}
  else if (el_pt> 25 && el_pt<= 30 && met> 160 && met<= 170) {eff = 0.538462; errup = 0.0770645; errdown = 0.078728;}
  else if (el_pt> 30 && el_pt<= 110 && met> 160 && met<= 170) {eff = 0.873646; errup = 0.0142978; errdown = 0.0157315;}
  else if (el_pt> 110 && el_pt<= 120 && met> 160 && met<= 170) {eff = 0.827586; errup = 0.0723692; errdown = 0.0998191;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 160 && met<= 170) {eff = 0.906667; errup = 0.0337754; errdown = 0.04665;}
  else if (el_pt> 20 && el_pt<= 25 && met> 170 && met<= 180) {eff = 0.696429; errup = 0.0657537; errdown = 0.0736929;}
  else if (el_pt> 25 && el_pt<= 30 && met> 170 && met<= 180) {eff = 0.734694; errup = 0.0670433; errdown = 0.0780196;}
  else if (el_pt> 30 && el_pt<= 110 && met> 170 && met<= 180) {eff = 0.859922; errup = 0.0155312; errdown = 0.0170189;}
  else if (el_pt> 110 && el_pt<= 120 && met> 170 && met<= 180) {eff = 0.911765; errup = 0.0476324; errdown = 0.0784416;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 170 && met<= 180) {eff = 0.974026; errup = 0.0167592; errdown = 0.0332329;}
  else if (el_pt> 20 && el_pt<= 25 && met> 180 && met<= 190) {eff = 0.816327; errup = 0.0574355; errdown = 0.0725146;}
  else if (el_pt> 25 && el_pt<= 30 && met> 180 && met<= 190) {eff = 0.837209; errup = 0.0580173; errdown = 0.0766334;}
  else if (el_pt> 30 && el_pt<= 110 && met> 180 && met<= 190) {eff = 0.919725; errup = 0.0131195; errdown = 0.0152117;}
  else if (el_pt> 110 && el_pt<= 120 && met> 180 && met<= 190) {eff = 0.913043; errup = 0.0559625; errdown = 0.103371;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 180 && met<= 190) {eff = 0.985294; errup = 0.0121686; errdown = 0.0330032;}
  else if (el_pt> 20 && el_pt<= 25 && met> 190 && met<= 200) {eff = 0.767442; errup = 0.0679853; errdown = 0.0824369;}
  else if (el_pt> 25 && el_pt<= 30 && met> 190 && met<= 200) {eff = 0.791667; errup = 0.0613414; errdown = 0.0754508;}
  else if (el_pt> 30 && el_pt<= 110 && met> 190 && met<= 200) {eff = 0.936471; errup = 0.0118929; errdown = 0.0141481;}
  else if (el_pt> 110 && el_pt<= 120 && met> 190 && met<= 200) {eff = 1; errup = 0; errdown = 0.0576587;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 190 && met<= 200) {eff = 1; errup = 0; errdown = 0.0263287;}
  else if (el_pt> 20 && el_pt<= 25 && met> 200 && met<= 210) {eff = 0.878049; errup = 0.0516458; errdown = 0.0740833;}
  else if (el_pt> 25 && el_pt<= 30 && met> 200 && met<= 210) {eff = 0.810811; errup = 0.0670146; errdown = 0.0869581;}
  else if (el_pt> 30 && el_pt<= 110 && met> 200 && met<= 210) {eff = 0.960526; errup = 0.00996633; errdown = 0.0127051;}
  else if (el_pt> 110 && el_pt<= 120 && met> 200 && met<= 210) {eff = 0.931034; errup = 0.0444185; errdown = 0.0838105;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 200 && met<= 210) {eff = 0.958904; errup = 0.0222832; errdown = 0.0383633;}
  else if (el_pt> 20 && el_pt<= 25 && met> 210 && met<= 9999) {eff = 0.935714; errup = 0.0207171; errdown = 0.0279785;}
  else if (el_pt> 25 && el_pt<= 30 && met> 210 && met<= 9999) {eff = 0.938053; errup = 0.0225646; errdown = 0.0317723;}
  else if (el_pt> 30 && el_pt<= 110 && met> 210 && met<= 9999) {eff = 0.965116; errup = 0.00529352; errdown = 0.00612651;}
  else if (el_pt> 110 && el_pt<= 120 && met> 210 && met<= 9999) {eff = 0.976744; errup = 0.015007; errdown = 0.0298505;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 210 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00719369;}
  return eff;
});

const NamedFunc get_1mu_trigeff2017_v0("get_1mu_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 110) {eff = 0.140487; errup = 0.0069928; errdown = 0.00672041;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 110) {eff = 0.509969; errup = 0.0117043; errdown = 0.011715;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 110) {eff = 0.664176; errup = 0.00668724; errdown = 0.00675275;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 110) {eff = 0.915724; errup = 0.00327313; errdown = 0.00339022;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 110 && met<= 120) {eff = 0.397661; errup = 0.0408696; errdown = 0.039584;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 110 && met<= 120) {eff = 0.709924; errup = 0.0415572; errdown = 0.0450534;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 110 && met<= 120) {eff = 0.833773; errup = 0.0194869; errdown = 0.0213668;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 110 && met<= 120) {eff = 0.966499; errup = 0.00735202; errdown = 0.009091;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 120 && met<= 130) {eff = 0.447853; errup = 0.042223; errdown = 0.0415354;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 120 && met<= 130) {eff = 0.723881; errup = 0.0403806; errdown = 0.0440273;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 120 && met<= 130) {eff = 0.86921; errup = 0.0178717; errdown = 0.0200359;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 120 && met<= 130) {eff = 0.958084; errup = 0.00895556; errdown = 0.0109858;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 130 && met<= 140) {eff = 0.575; errup = 0.0483942; errdown = 0.0497523;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 130 && met<= 140) {eff = 0.762887; errup = 0.0451032; errdown = 0.0511298;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 130 && met<= 140) {eff = 0.852778; errup = 0.0189996; errdown = 0.0211017;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 130 && met<= 140) {eff = 0.969072; errup = 0.00782771; errdown = 0.0100131;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 140 && met<= 150) {eff = 0.554545; errup = 0.051138; errdown = 0.0522185;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 140 && met<= 150) {eff = 0.768116; errup = 0.0532749; errdown = 0.0620693;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 140 && met<= 150) {eff = 0.88755; errup = 0.0203114; errdown = 0.023728;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 140 && met<= 150) {eff = 0.973214; errup = 0.0075712; errdown = 0.00998961;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 160) {eff = 0.666667; errup = 0.0443388; errdown = 0.0472195;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 160) {eff = 0.866667; errup = 0.0366017; errdown = 0.0459433;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 160) {eff = 0.945607; errup = 0.0146741; errdown = 0.0189235;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 160) {eff = 0.99; errup = 0.00477986; errdown = 0.00783614;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 160 && met<= 170) {eff = 0.771429; errup = 0.0525767; errdown = 0.0613516;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 160 && met<= 170) {eff = 0.851852; errup = 0.0404926; errdown = 0.0504559;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 160 && met<= 170) {eff = 0.92766; errup = 0.0170011; errdown = 0.0210882;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 160 && met<= 170) {eff = 0.972637; errup = 0.00806857; errdown = 0.0107752;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 170 && met<= 180) {eff = 0.802469; errup = 0.0459391; errdown = 0.0543801;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 170 && met<= 180) {eff = 0.888889; errup = 0.0354029; errdown = 0.0466071;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 170 && met<= 180) {eff = 0.963303; errup = 0.0126003; errdown = 0.0176003;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 170 && met<= 180) {eff = 0.993671; errup = 0.00408698; errdown = 0.0082865;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 180 && met<= 190) {eff = 0.714286; errup = 0.0643308; errdown = 0.0730111;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 180 && met<= 190) {eff = 0.969697; errup = 0.0195489; errdown = 0.0385739;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 180 && met<= 190) {eff = 0.940887; errup = 0.0165669; errdown = 0.0215411;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 180 && met<= 190) {eff = 0.983165; errup = 0.00725289; errdown = 0.0112281;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 190 && met<= 200) {eff = 0.857143; errup = 0.0478753; errdown = 0.0628782;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 190 && met<= 200) {eff = 0.94; errup = 0.0324767; errdown = 0.054936;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 190 && met<= 200) {eff = 0.975309; errup = 0.0117802; errdown = 0.0190923;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 190 && met<= 200) {eff = 0.992647; errup = 0.00474792; errdown = 0.00961549;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 200 && met<= 210) {eff = 0.865385; errup = 0.0482807; errdown = 0.0649611;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 200 && met<= 210) {eff = 1; errup = 0; errdown = 0.0297298;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 200 && met<= 210) {eff = 0.992857; errup = 0.00590966; errdown = 0.0162324;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 200 && met<= 210) {eff = 0.985816; errup = 0.00677638; errdown = 0.0110731;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 210 && met<= 9999) {eff = 0.938547; errup = 0.0179637; errdown = 0.0236054;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 210 && met<= 9999) {eff = 0.97191; errup = 0.0120805; errdown = 0.0185555;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 210 && met<= 9999) {eff = 0.983264; errup = 0.00577123; errdown = 0.00815023;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 210 && met<= 9999) {eff = 0.994152; errup = 0.0025237; errdown = 0.00393672;}
  return eff;
});

const NamedFunc get_2el_trigeff2017_v0("get_2el_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.846154; errup = 0.0353795; errdown = 0.0353795;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.875; errup = 0.0268248; errdown = 0.0268248;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.888889; errup = 0.0240328; errdown = 0.0240328;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.911765; errup = 0.0217539; errdown = 0.0217539;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.872222; errup = 0.0248831; errdown = 0.0248831;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.862857; errup = 0.0260038; errdown = 0.0260038;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.897436; errup = 0.0242905; errdown = 0.0242905;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.926667; errup = 0.0212847; errdown = 0.0212847;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.881481; errup = 0.0278184; errdown = 0.0278184;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.900763; errup = 0.026122; errdown = 0.026122;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.925532; errup = 0.027078; errdown = 0.027078;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.897727; errup = 0.0323006; errdown = 0.0323006;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.93617; errup = 0.025213; errdown = 0.025213;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.880435; errup = 0.0338265; errdown = 0.0338265;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.975464; errup = 0.00289636; errdown = 0.00289636;}
  return eff;
});

const NamedFunc get_2mu_trigeff2017_v0("get_2mu_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.932039; errup = 0.0175353; errdown = 0.0175353;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.952381; errup = 0.011271; errdown = 0.011271;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.977828; errup = 0.000485423; errdown = 0.000485423;}
  return eff;
});

const NamedFunc get_0l_trigeff2018_v0("get_0l_trigeff2018_v0", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_0l_fakemet_trigeff2018_v0("get_0l_fakemet_trigeff2018_v0", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_1el_trigeff2018_v0("get_1el_trigeff2018_v0", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_1mu_trigeff2018_v0("get_1mu_trigeff2018_v0", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_2el_trigeff2018_v0("get_2el_trigeff2018_v0", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_2mu_trigeff2018_v0("get_2mu_trigeff2018_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.965812; errup = 0.0118789; errdown = 0.0118789;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.935754; errup = 0.0129587; errdown = 0.0129587;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.984612; errup = 0.000361107; errdown = 0.000361107;}
  return eff;
});

}
