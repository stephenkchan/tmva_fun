#include "rfli_plot.h"
// #include "TMVA/correlations.h"
const bool lxplus=true;

//Make the ranking plots and store the names of the TMVA .root files so we can access them later to make correlations, ROC plots, input plots, etc.
int main(int argc,char**argv){
  string top="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20170412/",pdir="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/plots/tmva_fun/20170412/";
  int tF=2;
  if(argc>1)top  =argv[1];
  if(argc>2)pdir =argv[2];
  if(argc>3)tF   =atoi(argv[3]);
  rfli_plot navis(top,top+"/output/");vector<TString>tmva_files=navis.all_plots(pdir,tF);
  for(const auto&f:tmva_files)cout<<f<<endl;
  return 0;
}
