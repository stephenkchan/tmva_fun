#include "rfli_plot.h"
// #include "TMVA/correlations.h"
const bool lxplus=true;

//Make the ranking plots and store the names of the TMVA .root files so we can access them later to make correlations, ROC plots, input plots, etc.
int main(int argc,char**argv){
  string top="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20170412/",pdir="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/plots/tmva_fun/20170412/";
  int tF=2,whic=1,nj=2,ptv=2;
  if(!(argc==7||argc==8))std::cerr<<"You're probably using this wrong.  Defaults only exist for easy testing; if you're trying to use this for production, provide 6 extra arguments after the executable name.";
  if(argc>1)top  =argv[1];
  if(argc>2)pdir =argv[2];
  if(argc>3)tF   =atoi(argv[3]);
  if(argc>4)whic =atoi(argv[4]);
  if(argc>5)nj   =atoi(argv[5]);
  if(argc>6)ptv  =atoi(argv[6]);
  rfli_plot navis(top,top+"output/",ptv,nj,false);
  //have the extra be an override to force plot making
  if(argc>7) navis.single_plots(whic,nj,ptv,pdir,tF,true);
  else navis.single_plots_log_check(whic,nj,ptv,pdir,tF,true);
  return 0;
}
