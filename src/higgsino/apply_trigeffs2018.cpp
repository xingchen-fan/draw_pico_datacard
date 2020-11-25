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
  if (ht> 0 && ht<= 250 && met> 150 && met<= 155) {eff = 0.0967153; errup = 0.0150137; errdown = 0.0145224;}
  else if (ht> 250 && ht<= 350 && met> 150 && met<= 155) {eff = 0.167497; errup = 0.0211652; errdown = 0.0210723;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.241638; errup = 0.0301657; errdown = 0.0300795;}
  else if (ht> 450 && ht<= 600 && met> 150 && met<= 155) {eff = 0.296256; errup = 0.0366391; errdown = 0.03657;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.329746; errup = 0.0417888; errdown = 0.0416636;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.324786; errup = 0.0467917; errdown = 0.0461939;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.355705; errup = 0.0603952; errdown = 0.0589062;}
  else if (ht> 0 && ht<= 250 && met> 155 && met<= 160) {eff = 0.113527; errup = 0.0180531; errdown = 0.0174079;}
  else if (ht> 250 && ht<= 350 && met> 155 && met<= 160) {eff = 0.206503; errup = 0.0258924; errdown = 0.0258034;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.303987; errup = 0.0375402; errdown = 0.0374743;}
  else if (ht> 450 && ht<= 600 && met> 155 && met<= 160) {eff = 0.343447; errup = 0.0422467; errdown = 0.0421913;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.386041; errup = 0.0484809; errdown = 0.0483928;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.408163; errup = 0.0557802; errdown = 0.0554969;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.557692; errup = 0.0842001; errdown = 0.0849625;}
  else if (ht> 0 && ht<= 250 && met> 160 && met<= 165) {eff = 0.150068; errup = 0.0212125; errdown = 0.0205461;}
  else if (ht> 250 && ht<= 350 && met> 160 && met<= 165) {eff = 0.238877; errup = 0.0268564; errdown = 0.0267496;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.320253; errup = 0.035451; errdown = 0.0353715;}
  else if (ht> 450 && ht<= 600 && met> 160 && met<= 165) {eff = 0.378999; errup = 0.0415648; errdown = 0.0415101;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.450299; errup = 0.050112; errdown = 0.0500683;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.47076; errup = 0.0566384; errdown = 0.056548;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.442478; errup = 0.0690247; errdown = 0.0682026;}
  else if (ht> 0 && ht<= 250 && met> 165 && met<= 170) {eff = 0.167954; errup = 0.0251208; errdown = 0.0241645;}
  else if (ht> 250 && ht<= 350 && met> 165 && met<= 170) {eff = 0.305413; errup = 0.0337355; errdown = 0.0336594;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.400266; errup = 0.0436087; errdown = 0.0435682;}
  else if (ht> 450 && ht<= 600 && met> 165 && met<= 170) {eff = 0.453343; errup = 0.0491396; errdown = 0.0491197;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.50641; errup = 0.0558154; errdown = 0.055821;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.52459; errup = 0.0622959; errdown = 0.0623782;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.514019; errup = 0.074971; errdown = 0.0751716;}
  else if (ht> 0 && ht<= 250 && met> 170 && met<= 175) {eff = 0.233546; errup = 0.0321916; errdown = 0.0314222;}
  else if (ht> 250 && ht<= 350 && met> 170 && met<= 175) {eff = 0.354263; errup = 0.0389686; errdown = 0.0389045;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.447485; errup = 0.048566; errdown = 0.0485431;}
  else if (ht> 450 && ht<= 600 && met> 170 && met<= 175) {eff = 0.488976; errup = 0.0528465; errdown = 0.0528416;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.556686; errup = 0.0610976; errdown = 0.0611526;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.580645; errup = 0.0678513; errdown = 0.0681328;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.574074; errup = 0.0786405; errdown = 0.0796228;}
  else if (ht> 0 && ht<= 250 && met> 175 && met<= 180) {eff = 0.213115; errup = 0.0323576; errdown = 0.0311651;}
  else if (ht> 250 && ht<= 350 && met> 175 && met<= 180) {eff = 0.381985; errup = 0.0422952; errdown = 0.042226;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.4975; errup = 0.0538146; errdown = 0.0538134;}
  else if (ht> 450 && ht<= 600 && met> 175 && met<= 180) {eff = 0.56453; errup = 0.060499; errdown = 0.0605259;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.599411; errup = 0.0652644; errdown = 0.0653553;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.626506; errup = 0.0725999; errdown = 0.0730831;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.577982; errup = 0.078725; errdown = 0.0797419;}
  else if (ht> 0 && ht<= 250 && met> 180 && met<= 185) {eff = 0.242537; errup = 0.0365664; errdown = 0.0349584;}
  else if (ht> 250 && ht<= 350 && met> 180 && met<= 185) {eff = 0.468608; errup = 0.0456795; errdown = 0.0456595;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.54512; errup = 0.0523274; errdown = 0.0523529;}
  else if (ht> 450 && ht<= 600 && met> 180 && met<= 185) {eff = 0.628895; errup = 0.0596188; errdown = 0.0596832;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.678516; errup = 0.0649901; errdown = 0.0651589;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.708333; errup = 0.0724515; errdown = 0.0733919;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.72619; errup = 0.0840402; errdown = 0.0878291;}
  else if (ht> 0 && ht<= 250 && met> 185 && met<= 190) {eff = 0.256637; errup = 0.0400091; errdown = 0.038164;}
  else if (ht> 250 && ht<= 350 && met> 185 && met<= 190) {eff = 0.476334; errup = 0.0466341; errdown = 0.0466169;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.604627; errup = 0.0577631; errdown = 0.0578274;}
  else if (ht> 450 && ht<= 600 && met> 185 && met<= 190) {eff = 0.686888; errup = 0.0646849; errdown = 0.0647724;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.722523; errup = 0.0690384; errdown = 0.0692797;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.776699; errup = 0.077237; errdown = 0.0784073;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.8; errup = 0.0868232; errdown = 0.0916411;}
  else if (ht> 0 && ht<= 250 && met> 190 && met<= 195) {eff = 0.318182; errup = 0.0486661; errdown = 0.0469038;}
  else if (ht> 250 && ht<= 350 && met> 190 && met<= 195) {eff = 0.507812; errup = 0.0496426; errdown = 0.0496489;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.664678; errup = 0.0631778; errdown = 0.0632866;}
  else if (ht> 450 && ht<= 600 && met> 190 && met<= 195) {eff = 0.710407; errup = 0.0669595; errdown = 0.0670759;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.76234; errup = 0.0723177; errdown = 0.072583;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.693548; errup = 0.0727106; errdown = 0.0738229;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.785714; errup = 0.088376; errdown = 0.094047;}
  else if (ht> 0 && ht<= 250 && met> 195 && met<= 200) {eff = 0.392593; errup = 0.0587338; errdown = 0.0573817;}
  else if (ht> 250 && ht<= 350 && met> 195 && met<= 200) {eff = 0.575663; errup = 0.0564932; errdown = 0.0565812;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.717647; errup = 0.0678563; errdown = 0.0680036;}
  else if (ht> 450 && ht<= 600 && met> 195 && met<= 200) {eff = 0.768065; errup = 0.0719207; errdown = 0.072056;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.791075; errup = 0.0748891; errdown = 0.0752093;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.784211; errup = 0.0782049; errdown = 0.0795381;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.814815; errup = 0.0870468; errdown = 0.0918832;}
  else if (ht> 0 && ht<= 250 && met> 200 && met<= 210) {eff = 0.447115; errup = 0.0460682; errdown = 0.0456331;}
  else if (ht> 250 && ht<= 350 && met> 200 && met<= 210) {eff = 0.669557; errup = 0.0436146; errdown = 0.043731;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.73572; errup = 0.0467041; errdown = 0.0468034;}
  else if (ht> 450 && ht<= 600 && met> 200 && met<= 210) {eff = 0.796218; errup = 0.0499228; errdown = 0.0500172;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.844789; errup = 0.053137; errdown = 0.0533263;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.797251; errup = 0.0544585; errdown = 0.0554635;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.849315; errup = 0.0601831; errdown = 0.0630276;}
  else if (ht> 0 && ht<= 250 && met> 210 && met<= 220) {eff = 0.483051; errup = 0.0582314; errdown = 0.0579628;}
  else if (ht> 250 && ht<= 350 && met> 210 && met<= 220) {eff = 0.707895; errup = 0.0464989; errdown = 0.046708;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.823691; errup = 0.0517625; errdown = 0.0519056;}
  else if (ht> 450 && ht<= 600 && met> 210 && met<= 220) {eff = 0.853846; errup = 0.0531989; errdown = 0.053307;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.86014; errup = 0.0542692; errdown = 0.0545357;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.924528; errup = 0.0589093; errdown = 0.0599973;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.881818; errup = 0.0624045; errdown = 0.0667303;}
  else if (ht> 0 && ht<= 250 && met> 220 && met<= 230) {eff = 0.5875; errup = 0.0696658; errdown = 0.0717458;}
  else if (ht> 250 && ht<= 350 && met> 220 && met<= 230) {eff = 0.762238; errup = 0.0500787; errdown = 0.0504354;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.890011; errup = 0.055502; errdown = 0.0556851;}
  else if (ht> 450 && ht<= 600 && met> 220 && met<= 230) {eff = 0.905991; errup = 0.0561762; errdown = 0.0563066;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.897516; errup = 0.0562536; errdown = 0.0565501;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.877049; errup = 0.0577854; errdown = 0.0591199;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.921569; errup = 0.0624077; errdown = 0.0671418;}
  else if (ht> 0 && ht<= 250 && met> 230 && met<= 240) {eff = 0.759259; errup = 0.0821659; errdown = 0.0906687;}
  else if (ht> 250 && ht<= 350 && met> 230 && met<= 240) {eff = 0.809045; errup = 0.0616139; errdown = 0.062172;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.9; errup = 0.0657743; errdown = 0.0659915;}
  else if (ht> 450 && ht<= 600 && met> 230 && met<= 240) {eff = 0.915713; errup = 0.0665311; errdown = 0.0666614;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.965074; errup = 0.0349265; errdown = 0.0701522;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.948357; errup = 0.0516432; errdown = 0.07113;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.936364; errup = 0.0636364; errdown = 0.0748697;}
  else if (ht> 0 && ht<= 250 && met> 240 && met<= 250) {eff = 0.757576; errup = 0.0961604; errdown = 0.111793;}
  else if (ht> 250 && ht<= 350 && met> 240 && met<= 250) {eff = 0.888158; errup = 0.066507; errdown = 0.0673262;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.945899; errup = 0.0541012; errdown = 0.0690012;}
  else if (ht> 450 && ht<= 600 && met> 240 && met<= 250) {eff = 0.949398; errup = 0.0506024; errdown = 0.0689052;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.953722; errup = 0.0462777; errdown = 0.0696014;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.956098; errup = 0.0439024; errdown = 0.0715069;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.952381; errup = 0.047619; errdown = 0.0774568;}
  else if (ht> 0 && ht<= 250 && met> 250 && met<= 275) {eff = 0.75; errup = 0.0952577; errdown = 0.12006;}
  else if (ht> 250 && ht<= 350 && met> 250 && met<= 275) {eff = 0.923469; errup = 0.0210724; errdown = 0.0226616;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 275) {eff = 0.968597; errup = 0.0178036; errdown = 0.0181361;}
  else if (ht> 450 && ht<= 600 && met> 250 && met<= 275) {eff = 0.975101; errup = 0.0175303; errdown = 0.0177027;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.977588; errup = 0.0177678; errdown = 0.0181074;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.983092; errup = 0.0169082; errdown = 0.0194056;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.970588; errup = 0.0211717; errdown = 0.0257862;}
  else if (ht> 0 && ht<= 250 && met> 275 && met<= 300) {eff = 0.571429; errup = 0.222712; errdown = 0.248043;}
  else if (ht> 250 && ht<= 350 && met> 275 && met<= 300) {eff = 0.960265; errup = 0.0229593; errdown = 0.0284651;}
  else if (ht> 350 && ht<= 450 && met> 275 && met<= 300) {eff = 0.988679; errup = 0.0113208; errdown = 0.0185487;}
  else if (ht> 450 && ht<= 600 && met> 275 && met<= 300) {eff = 0.985045; errup = 0.0149551; errdown = 0.0179142;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.989262; errup = 0.0107383; errdown = 0.018086;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.993174; errup = 0.00682594; errdown = 0.0195361;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.981481; errup = 0.0185185; errdown = 0.0294305;}
  else if (ht> 0 && ht<= 250 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.369032;}
  else if (ht> 250 && ht<= 350 && met> 300 && met<= 9999) {eff = 0.97561; errup = 0.0188677; errdown = 0.0329513;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.994859; errup = 0.00514139; errdown = 0.0125718;}
  else if (ht> 450 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.996457; errup = 0.00354296; errdown = 0.0109897;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.990426; errup = 0.00957447; errdown = 0.0114221;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.997326; errup = 0.0026738; errdown = 0.0122737;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0144134;}
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
