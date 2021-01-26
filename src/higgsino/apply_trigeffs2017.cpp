#include <vector>
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/apply_trigeffs2017.hpp"

namespace Higfuncs{

const NamedFunc get_0l_trigeff2017("get_0l_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.123457; errup = 0.0511779; errdown = 0.0511779;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.151724; errup = 0.0624119; errdown = 0.0624119;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.228311; errup = 0.0857248; errdown = 0.0857248;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.228395; errup = 0.0873896; errdown = 0.0873896;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 180) {eff = 0.251366; errup = 0.0946626; errdown = 0.0946626;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 190) {eff = 0.285714; errup = 0.103601; errdown = 0.103601;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 200) {eff = 0.333333; errup = 0.127092; errdown = 0.127092;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 9999) {eff = 0.574468; errup = 0.130371; errdown = 0.130371;}
  else if (ht> 200 && ht<= 300 && met> 150 && met<= 155) {eff = 0.297837; errup = 0.116824; errdown = 0.116824;}
  else if (ht> 200 && ht<= 300 && met> 155 && met<= 160) {eff = 0.347426; errup = 0.136067; errdown = 0.136067;}
  else if (ht> 200 && ht<= 300 && met> 160 && met<= 165) {eff = 0.385093; errup = 0.138234; errdown = 0.138234;}
  else if (ht> 200 && ht<= 300 && met> 165 && met<= 170) {eff = 0.454955; errup = 0.162926; errdown = 0.162926;}
  else if (ht> 200 && ht<= 300 && met> 170 && met<= 175) {eff = 0.497354; errup = 0.178092; errdown = 0.178092;}
  else if (ht> 200 && ht<= 300 && met> 175 && met<= 180) {eff = 0.582192; errup = 0.208295; errdown = 0.208295;}
  else if (ht> 200 && ht<= 300 && met> 180 && met<= 185) {eff = 0.656934; errup = 0.223291; errdown = 0.223291;}
  else if (ht> 200 && ht<= 300 && met> 185 && met<= 190) {eff = 0.663551; errup = 0.225992; errdown = 0.225992;}
  else if (ht> 200 && ht<= 300 && met> 190 && met<= 195) {eff = 0.716495; errup = 0.243676; errdown = 0.243676;}
  else if (ht> 200 && ht<= 300 && met> 195 && met<= 200) {eff = 0.789189; errup = 0.210811; errdown = 0.267707;}
  else if (ht> 200 && ht<= 300 && met> 200 && met<= 210) {eff = 0.857143; errup = 0.142857; errdown = 0.180304;}
  else if (ht> 200 && ht<= 300 && met> 210 && met<= 220) {eff = 0.895028; errup = 0.104972; errdown = 0.188319;}
  else if (ht> 200 && ht<= 300 && met> 220 && met<= 230) {eff = 0.89726; errup = 0.10274; errdown = 0.189079;}
  else if (ht> 200 && ht<= 300 && met> 230 && met<= 240) {eff = 0.932039; errup = 0.0679612; errdown = 0.0964633;}
  else if (ht> 200 && ht<= 300 && met> 240 && met<= 250) {eff = 0.914286; errup = 0.0857143; errdown = 0.0973746;}
  else if (ht> 200 && ht<= 300 && met> 250 && met<= 275) {eff = 0.9375; errup = 0.0331599; errdown = 0.0331599;}
  else if (ht> 200 && ht<= 300 && met> 275 && met<= 9999) {eff = 0.977778; errup = 0.0222222; errdown = 0.0297024;}
  else if (ht> 300 && ht<= 400 && met> 150 && met<= 155) {eff = 0.419244; errup = 0.197077; errdown = 0.197077;}
  else if (ht> 300 && ht<= 400 && met> 155 && met<= 160) {eff = 0.440789; errup = 0.206929; errdown = 0.206929;}
  else if (ht> 300 && ht<= 400 && met> 160 && met<= 165) {eff = 0.565657; errup = 0.0858744; errdown = 0.0858744;}
  else if (ht> 300 && ht<= 400 && met> 165 && met<= 170) {eff = 0.5625; errup = 0.0866008; errdown = 0.0866008;}
  else if (ht> 300 && ht<= 400 && met> 170 && met<= 175) {eff = 0.609649; errup = 0.0929996; errdown = 0.0929996;}
  else if (ht> 300 && ht<= 400 && met> 175 && met<= 180) {eff = 0.696335; errup = 0.105018; errdown = 0.105018;}
  else if (ht> 300 && ht<= 400 && met> 180 && met<= 185) {eff = 0.73913; errup = 0.178528; errdown = 0.178528;}
  else if (ht> 300 && ht<= 400 && met> 185 && met<= 190) {eff = 0.783626; errup = 0.188336; errdown = 0.188336;}
  else if (ht> 300 && ht<= 400 && met> 190 && met<= 195) {eff = 0.804348; errup = 0.193564; errdown = 0.193564;}
  else if (ht> 300 && ht<= 400 && met> 195 && met<= 200) {eff = 0.845161; errup = 0.154839; errdown = 0.202364;}
  else if (ht> 300 && ht<= 400 && met> 200 && met<= 210) {eff = 0.88664; errup = 0.0943624; errdown = 0.0943624;}
  else if (ht> 300 && ht<= 400 && met> 210 && met<= 220) {eff = 0.925532; errup = 0.0744681; errdown = 0.098111;}
  else if (ht> 300 && ht<= 400 && met> 220 && met<= 230) {eff = 0.956522; errup = 0.0434783; errdown = 0.100737;}
  else if (ht> 300 && ht<= 400 && met> 230 && met<= 240) {eff = 0.978261; errup = 0.0217391; errdown = 0.0338717;}
  else if (ht> 300 && ht<= 400 && met> 240 && met<= 250) {eff = 0.964912; errup = 0.0350877; errdown = 0.0355423;}
  else if (ht> 300 && ht<= 400 && met> 250 && met<= 275) {eff = 0.970588; errup = 0.0227573; errdown = 0.0227573;}
  else if (ht> 300 && ht<= 400 && met> 275 && met<= 300) {eff = 0.988235; errup = 0.0117647; errdown = 0.0229914;}
  else if (ht> 300 && ht<= 400 && met> 300 && met<= 9999) {eff = 0.984848; errup = 0.0151515; errdown = 0.0912803;}
  else if (ht> 400 && ht<= 600 && met> 150 && met<= 155) {eff = 0.407002; errup = 0.0642462; errdown = 0.0642462;}
  else if (ht> 400 && ht<= 600 && met> 155 && met<= 160) {eff = 0.495192; errup = 0.0770015; errdown = 0.0770015;}
  else if (ht> 400 && ht<= 600 && met> 160 && met<= 165) {eff = 0.532374; errup = 0.0815655; errdown = 0.0815655;}
  else if (ht> 400 && ht<= 600 && met> 165 && met<= 170) {eff = 0.608808; errup = 0.0923941; errdown = 0.0923941;}
  else if (ht> 400 && ht<= 600 && met> 170 && met<= 175) {eff = 0.593548; errup = 0.0911365; errdown = 0.0911365;}
  else if (ht> 400 && ht<= 600 && met> 175 && met<= 180) {eff = 0.716981; errup = 0.107806; errdown = 0.107806;}
  else if (ht> 400 && ht<= 600 && met> 180 && met<= 185) {eff = 0.716475; errup = 0.0803522; errdown = 0.0803522;}
  else if (ht> 400 && ht<= 600 && met> 185 && met<= 190) {eff = 0.749049; errup = 0.0831922; errdown = 0.0831922;}
  else if (ht> 400 && ht<= 600 && met> 190 && met<= 195) {eff = 0.833333; errup = 0.0908854; errdown = 0.0908854;}
  else if (ht> 400 && ht<= 600 && met> 195 && met<= 200) {eff = 0.851282; errup = 0.0930868; errdown = 0.0930868;}
  else if (ht> 400 && ht<= 600 && met> 200 && met<= 210) {eff = 0.879795; errup = 0.089772; errdown = 0.089772;}
  else if (ht> 400 && ht<= 600 && met> 210 && met<= 220) {eff = 0.959375; errup = 0.040625; errdown = 0.0968661;}
  else if (ht> 400 && ht<= 600 && met> 220 && met<= 230) {eff = 0.954064; errup = 0.0459364; errdown = 0.0965082;}
  else if (ht> 400 && ht<= 600 && met> 230 && met<= 240) {eff = 0.953975; errup = 0.0213923; errdown = 0.0213923;}
  else if (ht> 400 && ht<= 600 && met> 240 && met<= 250) {eff = 0.969543; errup = 0.0208045; errdown = 0.0208045;}
  else if (ht> 400 && ht<= 600 && met> 250 && met<= 275) {eff = 0.978261; errup = 0.0217391; errdown = 0.0281273;}
  else if (ht> 400 && ht<= 600 && met> 275 && met<= 300) {eff = 0.98773; errup = 0.0122699; errdown = 0.0281309;}
  else if (ht> 400 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.998663; errup = 0.0013369; errdown = 0.0840732;}
  else if (ht> 600 && ht<= 950 && met> 150 && met<= 155) {eff = 0.415385; errup = 0.0737943; errdown = 0.0737943;}
  else if (ht> 600 && ht<= 950 && met> 155 && met<= 160) {eff = 0.518828; errup = 0.0899055; errdown = 0.0899055;}
  else if (ht> 600 && ht<= 950 && met> 160 && met<= 165) {eff = 0.555556; errup = 0.0702944; errdown = 0.0702944;}
  else if (ht> 600 && ht<= 950 && met> 165 && met<= 170) {eff = 0.580488; errup = 0.0726663; errdown = 0.0726663;}
  else if (ht> 600 && ht<= 950 && met> 170 && met<= 175) {eff = 0.615385; errup = 0.0768102; errdown = 0.0768102;}
  else if (ht> 600 && ht<= 950 && met> 175 && met<= 180) {eff = 0.716667; errup = 0.0858251; errdown = 0.0858251;}
  else if (ht> 600 && ht<= 950 && met> 180 && met<= 185) {eff = 0.675497; errup = 0.0734554; errdown = 0.0734554;}
  else if (ht> 600 && ht<= 950 && met> 185 && met<= 190) {eff = 0.795775; errup = 0.0813517; errdown = 0.0813517;}
  else if (ht> 600 && ht<= 950 && met> 190 && met<= 195) {eff = 0.742857; errup = 0.0783216; errdown = 0.0783216;}
  else if (ht> 600 && ht<= 950 && met> 195 && met<= 200) {eff = 0.843284; errup = 0.0844566; errdown = 0.0844566;}
  else if (ht> 600 && ht<= 950 && met> 200 && met<= 210) {eff = 0.885593; errup = 0.0605422; errdown = 0.0605422;}
  else if (ht> 600 && ht<= 950 && met> 210 && met<= 220) {eff = 0.921053; errup = 0.0623144; errdown = 0.0623144;}
  else if (ht> 600 && ht<= 950 && met> 220 && met<= 230) {eff = 0.940541; errup = 0.0594595; errdown = 0.0628678;}
  else if (ht> 600 && ht<= 950 && met> 230 && met<= 240) {eff = 0.947368; errup = 0.0261934; errdown = 0.0261934;}
  else if (ht> 600 && ht<= 950 && met> 240 && met<= 250) {eff = 0.983871; errup = 0.016129; errdown = 0.0215317;}
  else if (ht> 600 && ht<= 950 && met> 250 && met<= 275) {eff = 0.979167; errup = 0.0208333; errdown = 0.0244204;}
  else if (ht> 600 && ht<= 950 && met> 275 && met<= 300) {eff = 0.990148; errup = 0.00985222; errdown = 0.0241957;}
  else if (ht> 600 && ht<= 950 && met> 300 && met<= 9999) {eff = 0.997106; errup = 0.00289436; errdown = 0.0843903;}
  else if (ht> 950 && ht<= 9999 && met> 150 && met<= 160) {eff = 0.539216; errup = 0.0849513; errdown = 0.0849513;}
  else if (ht> 950 && ht<= 9999 && met> 160 && met<= 170) {eff = 0.567901; errup = 0.0880225; errdown = 0.0880225;}
  else if (ht> 950 && ht<= 9999 && met> 170 && met<= 180) {eff = 0.714286; errup = 0.101882; errdown = 0.101882;}
  else if (ht> 950 && ht<= 9999 && met> 180 && met<= 190) {eff = 0.825397; errup = 0.16131; errdown = 0.16131;}
  else if (ht> 950 && ht<= 9999 && met> 190 && met<= 200) {eff = 0.724138; errup = 0.147348; errdown = 0.147348;}
  else if (ht> 950 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.862745; errup = 0.119452; errdown = 0.119452;}
  else if (ht> 950 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.948718; errup = 0.0512821; errdown = 0.125276;}
  else if (ht> 950 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.875; errup = 0.122569; errdown = 0.122569;}
  else if (ht> 950 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.975; errup = 0.025; errdown = 0.032736;}
  else if (ht> 950 && ht<= 9999 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.0220517;}
  else if (ht> 950 && ht<= 9999 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.0198667;}
  else if (ht> 950 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.972222; errup = 0.0277778; errdown = 0.0335147;}
  else if (ht> 950 && ht<= 9999 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.089518;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_fakemet_trigeff2017("get_0l_fakemet_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 350 && met> 150 && met<= 160) {eff = 0.0974304; errup = 0.00970318; errdown = 0.00970318;}
  else if (ht> 0 && ht<= 350 && met> 160 && met<= 170) {eff = 0.14405; errup = 0.016044; errdown = 0.016044;}
  else if (ht> 0 && ht<= 350 && met> 170 && met<= 180) {eff = 0.184874; errup = 0.025163; errdown = 0.025163;}
  else if (ht> 0 && ht<= 350 && met> 180 && met<= 190) {eff = 0.321101; errup = 0.0447209; errdown = 0.0447209;}
  else if (ht> 0 && ht<= 350 && met> 190 && met<= 200) {eff = 0.348837; errup = 0.0513934; errdown = 0.0513934;}
  else if (ht> 0 && ht<= 350 && met> 200 && met<= 225) {eff = 0.545455; errup = 0.0567443; errdown = 0.0567443;}
  else if (ht> 0 && ht<= 350 && met> 225 && met<= 250) {eff = 0.62069; errup = 0.0901022; errdown = 0.0901022;}
  else if (ht> 0 && ht<= 350 && met> 250 && met<= 9999) {eff = 0.818182; errup = 0.116291; errdown = 0.116291;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.148239; errup = 0.0101691; errdown = 0.0101691;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.170517; errup = 0.012474; errdown = 0.012474;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.214286; errup = 0.0150635; errdown = 0.0150635;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.241993; errup = 0.0180663; errdown = 0.0180663;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.277533; errup = 0.0210154; errdown = 0.0210154;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.369444; errup = 0.0254381; errdown = 0.0254381;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.395833; errup = 0.0315667; errdown = 0.0315667;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.424731; errup = 0.036244; errdown = 0.036244;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.486188; errup = 0.0371505; errdown = 0.0371505;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.51049; errup = 0.0418029; errdown = 0.0418029;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.622642; errup = 0.0332911; errdown = 0.0332911;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.669173; errup = 0.0407985; errdown = 0.0407985;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.728261; errup = 0.0463795; errdown = 0.0463795;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.709677; errup = 0.0576468; errdown = 0.0576468;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.757576; errup = 0.0746009; errdown = 0.0746009;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 300) {eff = 0.916031; errup = 0.0242315; errdown = 0.0242315;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.954545; errup = 0.0314022; errdown = 0.0314022;}
  else if (ht> 450 && ht<= 550 && met> 150 && met<= 155) {eff = 0.189972; errup = 0.00841328; errdown = 0.00841328;}
  else if (ht> 450 && ht<= 550 && met> 155 && met<= 160) {eff = 0.228241; errup = 0.0103543; errdown = 0.0103543;}
  else if (ht> 450 && ht<= 550 && met> 160 && met<= 165) {eff = 0.243159; errup = 0.0119953; errdown = 0.0119953;}
  else if (ht> 450 && ht<= 550 && met> 165 && met<= 170) {eff = 0.323308; errup = 0.0153295; errdown = 0.0153295;}
  else if (ht> 450 && ht<= 550 && met> 170 && met<= 175) {eff = 0.377577; errup = 0.0174026; errdown = 0.0174026;}
  else if (ht> 450 && ht<= 550 && met> 175 && met<= 180) {eff = 0.407874; errup = 0.0195022; errdown = 0.0195022;}
  else if (ht> 450 && ht<= 550 && met> 180 && met<= 185) {eff = 0.418468; errup = 0.0218655; errdown = 0.0218655;}
  else if (ht> 450 && ht<= 550 && met> 185 && met<= 190) {eff = 0.486339; errup = 0.0261257; errdown = 0.0261257;}
  else if (ht> 450 && ht<= 550 && met> 190 && met<= 195) {eff = 0.454286; errup = 0.0266142; errdown = 0.0266142;}
  else if (ht> 450 && ht<= 550 && met> 195 && met<= 200) {eff = 0.582143; errup = 0.0294747; errdown = 0.0294747;}
  else if (ht> 450 && ht<= 550 && met> 200 && met<= 210) {eff = 0.590571; errup = 0.0244947; errdown = 0.0244947;}
  else if (ht> 450 && ht<= 550 && met> 210 && met<= 220) {eff = 0.671141; errup = 0.0272147; errdown = 0.0272147;}
  else if (ht> 450 && ht<= 550 && met> 220 && met<= 230) {eff = 0.795349; errup = 0.0275148; errdown = 0.0275148;}
  else if (ht> 450 && ht<= 550 && met> 230 && met<= 240) {eff = 0.823944; errup = 0.0319617; errdown = 0.0319617;}
  else if (ht> 450 && ht<= 550 && met> 240 && met<= 250) {eff = 0.796875; errup = 0.0355608; errdown = 0.0355608;}
  else if (ht> 450 && ht<= 550 && met> 250 && met<= 300) {eff = 0.911972; errup = 0.0168129; errdown = 0.0168129;}
  else if (ht> 450 && ht<= 550 && met> 300 && met<= 400) {eff = 0.975309; errup = 0.0121923; errdown = 0.0121923;}
  else if (ht> 450 && ht<= 550 && met> 400 && met<= 9999) {eff = 0.966667; errup = 0.0231741; errdown = 0.0231741;}
  else if (ht> 550 && ht<= 650 && met> 150 && met<= 155) {eff = 0.245086; errup = 0.00879666; errdown = 0.00879666;}
  else if (ht> 550 && ht<= 650 && met> 155 && met<= 160) {eff = 0.279679; errup = 0.0103794; errdown = 0.0103794;}
  else if (ht> 550 && ht<= 650 && met> 160 && met<= 165) {eff = 0.319121; errup = 0.0118475; errdown = 0.0118475;}
  else if (ht> 550 && ht<= 650 && met> 165 && met<= 170) {eff = 0.363946; errup = 0.0140301; errdown = 0.0140301;}
  else if (ht> 550 && ht<= 650 && met> 170 && met<= 175) {eff = 0.376166; errup = 0.0155941; errdown = 0.0155941;}
  else if (ht> 550 && ht<= 650 && met> 175 && met<= 180) {eff = 0.445813; errup = 0.0174432; errdown = 0.0174432;}
  else if (ht> 550 && ht<= 650 && met> 180 && met<= 185) {eff = 0.507184; errup = 0.0189505; errdown = 0.0189505;}
  else if (ht> 550 && ht<= 650 && met> 185 && met<= 190) {eff = 0.532957; errup = 0.021651; errdown = 0.021651;}
  else if (ht> 550 && ht<= 650 && met> 190 && met<= 195) {eff = 0.599109; errup = 0.0231283; errdown = 0.0231283;}
  else if (ht> 550 && ht<= 650 && met> 195 && met<= 200) {eff = 0.631043; errup = 0.02434; errdown = 0.02434;}
  else if (ht> 550 && ht<= 650 && met> 200 && met<= 210) {eff = 0.685155; errup = 0.0187591; errdown = 0.0187591;}
  else if (ht> 550 && ht<= 650 && met> 210 && met<= 220) {eff = 0.694382; errup = 0.0218378; errdown = 0.0218378;}
  else if (ht> 550 && ht<= 650 && met> 220 && met<= 230) {eff = 0.740413; errup = 0.0238111; errdown = 0.0238111;}
  else if (ht> 550 && ht<= 650 && met> 230 && met<= 240) {eff = 0.8; errup = 0.0245718; errdown = 0.0245718;}
  else if (ht> 550 && ht<= 650 && met> 240 && met<= 250) {eff = 0.856436; errup = 0.0246715; errdown = 0.0246715;}
  else if (ht> 550 && ht<= 650 && met> 250 && met<= 275) {eff = 0.880342; errup = 0.0173238; errdown = 0.0173238;}
  else if (ht> 550 && ht<= 650 && met> 275 && met<= 300) {eff = 0.903382; errup = 0.0205343; errdown = 0.0205343;}
  else if (ht> 550 && ht<= 650 && met> 300 && met<= 400) {eff = 0.933333; errup = 0.0117589; errdown = 0.0117589;}
  else if (ht> 550 && ht<= 650 && met> 400 && met<= 9999) {eff = 0.978102; errup = 0.00510454; errdown = 0.00510454;}
  else if (ht> 650 && ht<= 800 && met> 150 && met<= 155) {eff = 0.322707; errup = 0.00710557; errdown = 0.00710557;}
  else if (ht> 650 && ht<= 800 && met> 155 && met<= 160) {eff = 0.366185; errup = 0.00810742; errdown = 0.00810742;}
  else if (ht> 650 && ht<= 800 && met> 160 && met<= 165) {eff = 0.406324; errup = 0.00905657; errdown = 0.00905657;}
  else if (ht> 650 && ht<= 800 && met> 165 && met<= 170) {eff = 0.439822; errup = 0.00997531; errdown = 0.00997531;}
  else if (ht> 650 && ht<= 800 && met> 170 && met<= 175) {eff = 0.463355; errup = 0.0110594; errdown = 0.0110594;}
  else if (ht> 650 && ht<= 800 && met> 175 && met<= 180) {eff = 0.515652; errup = 0.0119226; errdown = 0.0119226;}
  else if (ht> 650 && ht<= 800 && met> 180 && met<= 185) {eff = 0.554072; errup = 0.0128428; errdown = 0.0128428;}
  else if (ht> 650 && ht<= 800 && met> 185 && met<= 190) {eff = 0.598348; errup = 0.0134323; errdown = 0.0134323;}
  else if (ht> 650 && ht<= 800 && met> 190 && met<= 195) {eff = 0.642009; errup = 0.0144877; errdown = 0.0144877;}
  else if (ht> 650 && ht<= 800 && met> 195 && met<= 200) {eff = 0.696538; errup = 0.0146713; errdown = 0.0146713;}
  else if (ht> 650 && ht<= 800 && met> 200 && met<= 210) {eff = 0.715486; errup = 0.0110539; errdown = 0.0110539;}
  else if (ht> 650 && ht<= 800 && met> 210 && met<= 220) {eff = 0.781983; errup = 0.0113604; errdown = 0.0113604;}
  else if (ht> 650 && ht<= 800 && met> 220 && met<= 230) {eff = 0.827618; errup = 0.0114987; errdown = 0.0114987;}
  else if (ht> 650 && ht<= 800 && met> 230 && met<= 240) {eff = 0.839763; errup = 0.0115368; errdown = 0.0115368;}
  else if (ht> 650 && ht<= 800 && met> 240 && met<= 250) {eff = 0.884467; errup = 0.0114531; errdown = 0.0114531;}
  else if (ht> 650 && ht<= 800 && met> 250 && met<= 275) {eff = 0.898712; errup = 0.00730038; errdown = 0.00730038;}
  else if (ht> 650 && ht<= 800 && met> 275 && met<= 300) {eff = 0.930855; errup = 0.00691768; errdown = 0.00691768;}
  else if (ht> 650 && ht<= 800 && met> 300 && met<= 350) {eff = 0.941468; errup = 0.00522821; errdown = 0.00522821;}
  else if (ht> 650 && ht<= 800 && met> 350 && met<= 400) {eff = 0.959596; errup = 0.00528902; errdown = 0.00528902;}
  else if (ht> 650 && ht<= 800 && met> 400 && met<= 450) {eff = 0.976963; errup = 0.00485453; errdown = 0.00485453;}
  else if (ht> 650 && ht<= 800 && met> 450 && met<= 500) {eff = 0.991124; errup = 0.00360739; errdown = 0.00360739;}
  else if (ht> 650 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.973333; errup = 0.00515956; errdown = 0.00515956;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.489496; errup = 0.00369257; errdown = 0.00369257;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.522675; errup = 0.00395583; errdown = 0.00395583;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.568861; errup = 0.00421634; errdown = 0.00421634;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.592914; errup = 0.00450707; errdown = 0.00450707;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.635472; errup = 0.00469609; errdown = 0.00469609;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.674262; errup = 0.00491711; errdown = 0.00491711;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.686402; errup = 0.00517974; errdown = 0.00517974;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.726592; errup = 0.00534945; errdown = 0.00534945;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.737255; errup = 0.00551234; errdown = 0.00551234;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.775837; errup = 0.00559582; errdown = 0.00559582;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.802235; errup = 0.00410939; errdown = 0.00410939;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.83659; errup = 0.00423424; errdown = 0.00423424;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.849854; errup = 0.00442152; errdown = 0.00442152;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.880527; errup = 0.00454664; errdown = 0.00454664;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.893832; errup = 0.00466479; errdown = 0.00466479;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.916677; errup = 0.00306379; errdown = 0.00306379;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.936391; errup = 0.00323515; errdown = 0.00323515;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.941998; errup = 0.0028653; errdown = 0.0028653;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.952878; errup = 0.00376834; errdown = 0.00376834;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.963603; errup = 0.00465149; errdown = 0.00465149;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.954495; errup = 0.00694311; errdown = 0.00694311;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.983751; errup = 0.00328979; errdown = 0.00328979;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.387888; errup = 0.001472; errdown = 0.001472;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.417927; errup = 0.00160609; errdown = 0.00160609;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.451373; errup = 0.0017459; errdown = 0.0017459;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.475537; errup = 0.00189276; errdown = 0.00189276;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.512738; errup = 0.00204099; errdown = 0.00204099;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.53824; errup = 0.00218953; errdown = 0.00218953;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.567544; errup = 0.0023476; errdown = 0.0023476;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.598143; errup = 0.00249669; errdown = 0.00249669;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.624385; errup = 0.00264351; errdown = 0.00264351;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.649303; errup = 0.00279531; errdown = 0.00279531;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.681851; errup = 0.002131; errdown = 0.002131;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.729521; errup = 0.00231462; errdown = 0.00231462;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.762535; errup = 0.00251095; errdown = 0.00251095;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.798438; errup = 0.00267273; errdown = 0.00267273;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.821241; errup = 0.00285038; errdown = 0.00285038;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.856744; errup = 0.00196232; errdown = 0.00196232;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.886954; errup = 0.00223653; errdown = 0.00223653;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.916497; errup = 0.00185908; errdown = 0.00185908;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.936176; errup = 0.00234582; errdown = 0.00234582;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.947549; errup = 0.00292349; errdown = 0.00292349;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.954474; errup = 0.0037065; errdown = 0.0037065;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.967807; errup = 0.00238727; errdown = 0.00238727;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1el_trigeff2017("get_1el_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.00606061; errup = 0.00297522; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.252727; errup = 0.0137849; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.47488; errup = 0.0101789; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.0580311; errup = 0.00853854; errdown = 0.00756;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.273902; errup = 0.0169699; errdown = 0.0163626;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.512061; errup = 0.0119679; errdown = 0.0119814;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.239788; errup = 0.0123203; errdown = 0.0119141;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.432836; errup = 0.0161873; errdown = 0.0160495;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.610241; errup = 0.00967965; errdown = 0.00976557;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.368855; errup = 0.0057901; errdown = 0.00575302;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.58404; errup = 0.00788592; errdown = 0.0079285;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.718103; errup = 0.00425673; errdown = 0.00429573;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.358079; errup = 0.0165383; errdown = 0.0162179;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.574176; errup = 0.0270377; errdown = 0.0274661;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.721596; errup = 0.01343; errdown = 0.0138268;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.727555; errup = 0.00887959; errdown = 0.00906034;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.776404; errup = 0.014199; errdown = 0.0148447;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.830705; errup = 0.00707675; errdown = 0.00731432;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.0128205; errup = 0.00857983; errdown = 0.00552679;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.251185; errup = 0.0228961; errdown = 0.0216526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.515088; errup = 0.0135744; errdown = 0.0135961;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.07; errup = 0.0180052; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.359756; errup = 0.0284325; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.58209; errup = 0.0154441; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.313953; errup = 0.0239721; errdown = 0.0230631;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.564612; errup = 0.0229432; errdown = 0.0232113;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.668443; errup = 0.0116455; errdown = 0.0118501;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.478528; errup = 0.0103138; errdown = 0.0102959;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.690422; errup = 0.0112097; errdown = 0.0114324;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.76283; errup = 0.00509721; errdown = 0.00517276;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.404682; errup = 0.0303426; errdown = 0.0296689;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.722973; errup = 0.0383999; errdown = 0.0416756;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.766862; errup = 0.0165056; errdown = 0.0173229;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.73375; errup = 0.0159396; errdown = 0.0165468;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.814136; errup = 0.0203184; errdown = 0.0220684;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.843819; errup = 0.00860726; errdown = 0.0089989;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.0526316; errup = 0.0271984; errdown = 0.019208;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.295082; errup = 0.0375578; errdown = 0.0351452;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.577465; errup = 0.0167048; errdown = 0.0168781;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.196429; errup = 0.0448753; errdown = 0.0388449;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.467105; errup = 0.0438353; errdown = 0.0433692;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.620134; errup = 0.0182735; errdown = 0.018608;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.495098; errup = 0.037384; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.59375; errup = 0.03452; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.738619; errup = 0.0125052; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.584366; errup = 0.0156025; errdown = 0.0157683;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.737811; errup = 0.0147447; errdown = 0.0152789;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.796101; errup = 0.00566351; errdown = 0.00578094;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.607843; errup = 0.0518085; errdown = 0.0541227;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.77551; errup = 0.0439117; errdown = 0.050173;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.801205; errup = 0.0182366; errdown = 0.0195126;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.767516; errup = 0.0244905; errdown = 0.0263047;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.791878; errup = 0.0297937; errdown = 0.0330066;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.863636; errup = 0.00987432; errdown = 0.0104889;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.0338983; errup = 0.0429671; errdown = 0.021865;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.357895; errup = 0.0557905; errdown = 0.0525014;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.63606; errup = 0.0202397; errdown = 0.0207128;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.230769; errup = 0.0733859; errdown = 0.0614831;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.532468; errup = 0.0624335; errdown = 0.0633651;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.701613; errup = 0.0211067; errdown = 0.0219585;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.645161; errup = 0.0455298; errdown = 0.0480775;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.700787; errup = 0.0426733; errdown = 0.0461241;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.807339; errup = 0.0135565; errdown = 0.0142912;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.716075; errup = 0.0211522; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.797227; errup = 0.0170527; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.84023; errup = 0.0061067; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.785714; errup = 0.0664487; errdown = 0.0823485;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.891892; errup = 0.0510106; errdown = 0.077259;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.827922; errup = 0.0219648; errdown = 0.0242503;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.841667; errup = 0.0342095; errdown = 0.0405813;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.893443; errup = 0.0283253; errdown = 0.0356815;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.905502; errup = 0.0102026; errdown = 0.0112294;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.2; errup = 0.130119; errdown = 0.0931225;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.604167; errup = 0.0776099; errdown = 0.0825214;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.801508; errup = 0.0204268; errdown = 0.0220347;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.222222; errup = 0.109068; errdown = 0.0843892;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.575; errup = 0.0872414; errdown = 0.0915181;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.788406; errup = 0.0225238; errdown = 0.0243027;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.745098; errup = 0.0646476; errdown = 0.0756564;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.813333; errup = 0.0466366; errdown = 0.0561452;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.85; errup = 0.0150446; errdown = 0.0163204;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.815686; errup = 0.0248726; errdown = 0.027542;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.876254; errup = 0.0193393; errdown = 0.0220715;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.895428; errup = 0.00633502; errdown = 0.00668102;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.666667; errup = 0.106107; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.892857; errup = 0.0577336; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.885965; errup = 0.0213678; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.982456; errup = 0.0145177; errdown = 0.0391869;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.970588; errup = 0.0189746; errdown = 0.0374789;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.93956; errup = 0.0102474; errdown = 0.0120003;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.461538; errup = 0.172004; errdown = 0.164847;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.612903; errup = 0.0974957; errdown = 0.105934;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.865546; errup = 0.0225156; errdown = 0.0258705;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.454545; errup = 0.189662; errdown = 0.179582;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.666667; errup = 0.106107; errdown = 0.122517;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.883178; errup = 0.0223037; errdown = 0.0262543;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.75; errup = 0.0943497; errdown = 0.119341;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.906977; errup = 0.0439838; errdown = 0.0674608;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.914286; errup = 0.0132412; errdown = 0.0152122;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.873418; errup = 0.026946; errdown = 0.0322169;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.929648; errup = 0.0182155; errdown = 0.0231163;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.929748; errup = 0.00601652; errdown = 0.00650973;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.809524; errup = 0.0888157; errdown = 0.125184;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0709947;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.907801; errup = 0.0246102; errdown = 0.0312058;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.942857; errup = 0.0368224; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.975; errup = 0.0206905; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.951351; errup = 0.011195; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.8; errup = 0.106751; errdown = 0.157061;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.916667; errup = 0.0536391; errdown = 0.0995072;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.955645; errup = 0.0130229; errdown = 0.017253;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.769231; errup = 0.122762; errdown = 0.174724;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.869565; errup = 0.0701228; errdown = 0.110814;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.920833; errup = 0.0175659; errdown = 0.0214819;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.944444; errup = 0.046004; errdown = 0.116415;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.97619; errup = 0.0197048; errdown = 0.0526298;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.961039; errup = 0.00899222; errdown = 0.0112252;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.950413; errup = 0.0194949; errdown = 0.0284436;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.978022; errup = 0.0104893; errdown = 0.0170363;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.969503; errup = 0.00381731; errdown = 0.00430889;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.916667; errup = 0.0690403; errdown = 0.16652;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0802771;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.964706; errup = 0.013912; errdown = 0.0204859;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.966667; errup = 0.0275914; errdown = 0.072517;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0384134;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.980815; errup = 0.00661202; errdown = 0.00932515;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.913043; errup = 0.0559625; errdown = 0.103371;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.981651; errup = 0.00725455; errdown = 0.0107985;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0769247;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 0.984674; errup = 0.00732055; errdown = 0.0119517;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0659133;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.042887;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 0.997938; errup = 0.00170573; errdown = 0.00472518;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.990099; errup = 0.00819202; errdown = 0.0223979;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00964279;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.992439; errup = 0.00171528; errdown = 0.00215231;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.115502;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00804214;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.108691;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0283562;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.998397; errup = 0.00132575; errdown = 0.0036754;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1mu_trigeff2017("get_1mu_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.0708592; errup = 0.00848378; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.386643; errup = 0.0149712; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.79144; errup = 0.00809149; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.271298; errup = 0.0170541; errdown = 0.0164307;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.550075; errup = 0.0198882; errdown = 0.0200436;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.843672; errup = 0.0091337; errdown = 0.00957456;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.431312; errup = 0.012669; errdown = 0.0125819;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.704082; errup = 0.0120989; errdown = 0.0123844;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.886586; errup = 0.00491944; errdown = 0.00510812;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.887264; errup = 0.00786673; errdown = 0.00835697;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.92204; errup = 0.00652681; errdown = 0.00704377;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.936299; errup = 0.00313807; errdown = 0.00328508;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.0614251; errup = 0.0143316; errdown = 0.0119574;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.414583; errup = 0.0236888; errdown = 0.0233169;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.814067; errup = 0.00972747; errdown = 0.0101237;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.325843; errup = 0.0311932; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.614198; errup = 0.0281699; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.857381; errup = 0.0102039; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.543956; errup = 0.0190702; errdown = 0.0191953;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.780458; errup = 0.0146207; errdown = 0.0153252;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.893736; errup = 0.0056122; errdown = 0.00587802;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.921429; errup = 0.0093482; errdown = 0.0104152;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.952077; errup = 0.00699194; errdown = 0.00802783;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.940981; errup = 0.0035208; errdown = 0.00372265;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.15; errup = 0.0337209; errdown = 0.0288809;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.45098; errup = 0.0332565; errdown = 0.0328488;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.842715; errup = 0.0105919; errdown = 0.0111815;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.520325; errup = 0.0487155; errdown = 0.0490739;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.647668; errup = 0.0360589; errdown = 0.037699;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.861354; errup = 0.0115486; errdown = 0.0123748;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.666667; errup = 0.025072; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.806922; errup = 0.0171572; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.903462; errup = 0.00599439; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.941414; errup = 0.0106018; errdown = 0.0125537;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.954397; errup = 0.00843586; errdown = 0.0100577;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.95625; errup = 0.00337964; errdown = 0.00363757;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.14; errup = 0.0672437; errdown = 0.0501519;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.556452; errup = 0.0479412; errdown = 0.0489288;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.868726; errup = 0.0122566; errdown = 0.0132553;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.619048; errup = 0.0664109; errdown = 0.0706357;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.816514; errup = 0.0382778; errdown = 0.044765;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.919672; errup = 0.0110898; errdown = 0.0125669;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.73262; errup = 0.0336118; errdown = 0.0362961;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.855596; errup = 0.0215151; errdown = 0.0242958;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.931906; errup = 0.00593101; errdown = 0.0064273;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.956364; errup = 0.0122799; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.972569; errup = 0.00808855; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.968242; errup = 0.00341417; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.423077; errup = 0.116981; errdown = 0.110076;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.666667; errup = 0.0654592; errdown = 0.0717057;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.912046; errup = 0.0124958; errdown = 0.0141913;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.724138; errup = 0.089299; errdown = 0.107524;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.850746; errup = 0.044684; errdown = 0.0568246;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.943262; errup = 0.0112896; errdown = 0.0136025;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.94186; errup = 0.0248848; errdown = 0.0374166;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.928571; errup = 0.0191788; errdown = 0.024547;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.956884; errup = 0.00560047; errdown = 0.0063382;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.980392; errup = 0.0106526; errdown = 0.0187052;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.973913; errup = 0.0102998; errdown = 0.0152564;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.981845; errup = 0.00295218; errdown = 0.0034642;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.533333; errup = 0.152935; errdown = 0.15827;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.882353; errup = 0.0554381; errdown = 0.0832913;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.929412; errup = 0.01398; errdown = 0.0167685;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.769231; errup = 0.122762; errdown = 0.174724;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.904762; errup = 0.0450174; errdown = 0.0689193;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.981308; errup = 0.0073897; errdown = 0.0109973;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0283562;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.94898; errup = 0.0218627; errdown = 0.0330405;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.973348; errup = 0.00524904; errdown = 0.00636007;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.988636; errup = 0.00940245; errdown = 0.0256444;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.993711; errup = 0.0052034; errdown = 0.0143129;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.985469; errup = 0.00306041; errdown = 0.00377044;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.142229;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.888889; errup = 0.0598486; errdown = 0.0963981;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.98374; errup = 0.00643122; errdown = 0.00958564;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.947368; errup = 0.0435805; errdown = 0.110836;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0738409;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.988473; errup = 0.00550887; errdown = 0.00902056;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0392319;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0182418;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.990955; errup = 0.00295229; errdown = 0.00410363;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0249041;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.993711; errup = 0.0052034; errdown = 0.0143129;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.995131; errup = 0.00168324; errdown = 0.00239245;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0709947;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 0.992788; errup = 0.00392226; errdown = 0.00696502;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.15411;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0972223;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 0.99505; errup = 0.00319693; errdown = 0.00649192;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0542609;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0224723;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 0.994958; errup = 0.00199818; errdown = 0.00299934;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0461088;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.012972;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.997174; errup = 0.00112047; errdown = 0.0016842;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2el_trigeff2017("get_2el_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.897959; errup = 0.0249664; errdown = 0.0249664;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.949309; errup = 0.0148916; errdown = 0.0148916;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.939502; errup = 0.0142222; errdown = 0.0142222;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.930233; errup = 0.0158603; errdown = 0.0158603;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.919014; errup = 0.0161885; errdown = 0.0161885;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.931818; errup = 0.0155131; errdown = 0.0155131;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.94382; errup = 0.0140922; errdown = 0.0140922;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.973913; errup = 0.0105101; errdown = 0.0105101;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.941704; errup = 0.01569; errdown = 0.01569;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.960591; errup = 0.0136558; errdown = 0.0136558;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.919075; errup = 0.0207345; errdown = 0.0207345;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.941176; errup = 0.0180462; errdown = 0.0180462;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.972222; errup = 0.0122488; errdown = 0.0122488;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.934211; errup = 0.0201085; errdown = 0.0201085;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.983325; errup = 0.00197631; errdown = 0.00197631;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2mu_trigeff2017("get_2mu_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.948413; errup = 0.0139338; errdown = 0.0139338;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.963134; errup = 0.0090451; errdown = 0.0090451;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.983369; errup = 0.000433117; errdown = 0.000433117;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

}
