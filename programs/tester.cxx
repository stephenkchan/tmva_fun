#include "bdt_trainer.h"
#include "bdt_testing.h"
#include "bdt_validate.h"
#include "bdt_ranker.h"

int main(){
  string in="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20160510/eff_77/",out="/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/results/";
//   string in="/afs/cern.ch/work/s/stchan/vhbb/eos/atlas/user/s/stchan/tmva_fun/20160510/eff_77/",out="/afs/cern.ch/work/s/stchan/vhbb/tmva_fun/results/";
  bdt_ranker surge(in,out,{"j0_j1","angle_j0_j1"},{"l0_l1","angle_inv_lb0","j0_l0"});
  TCanvas*c=new TCanvas("derp","derp",500,300);
//   surge.finish_ranking("blarg",false,2);  surge.print_training_plot(c,"blarg",".",2);
//   return 0;
  vector<TString>vars={"mBB","dRBB","dPhiVBB","MV2c20B1","pTB2","MET","dEtaVBB","MV2c20B2","mLL","dEtaBB","pTB1","pTV"};
  bdt_base pikachu(in,out,"blarg");
//   TCanvas*c=new TCanvas("derp","derp",500,300);
//   cout<<"TRAIN "<<pikachu.significance(true).second<<endl;
  bdt_trainer pichu(vars,true,pikachu,false,1);//,string("cut"),vars);
  cout<<"VAL "<<pichu.print_bdt_plots(c,out).second<<endl;
  bdt_trainer raichu(vars,true,pikachu,false,2);//,string("cut"),vars);
  cout<<"VAL "<<raichu.print_bdt_plots(c,out).second<<endl;
  bdt_trainer choochoo(vars,true,pikachu,true);//,string("cut"),vars);
  cout<<"VAL "<<choochoo.print_bdt_plots(c,out).second<<endl;
//   pichu.print_validation_bdts(c,".",-1);
//   for(auto i:{2,4,5,6})bdt_testing(vars,false,pikachu).print_test_bdts(c,".",i,false);
  return 0;
}
