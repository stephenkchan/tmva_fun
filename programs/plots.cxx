#include "sig_plotter.h"
#include "correlations.C"
#include <dirent.h>
#include <cstdlib>

int main(int argc,char**argv){
  string ddir="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20160226_srli/output/",pdir=ddir+"plots/";
  if(argc>1)ddir=argv[1];
  if(argc>2)pdir=argv[2];
  //ensure slashes are at the end
  if(ddir.substr(ddir.length()-1)!="/")ddir+="/";
  if(pdir.substr(pdir.length()-1)!="/")pdir+="/";
  if(opendir(ddir.c_str())==NULL)return 0;//if no data, don't do anything
  if(opendir(pdir.c_str())==NULL)system(("mkdir "+pdir).c_str());
  sig_plotter puerh(false);
  //significance plots
  puerh.do_plots(ddir,pdir);
  //correlations
  for(int iv=puerh.n_std_min();iv<=puerh.n_li_cases();iv++){
    for(int bg=-1;bg<puerh.n_bg_cases();bg++){
      string ftag=string(puerh.li_tag(iv).Data())+string(puerh.bg_tag(bg).Data());
      TString fn=ddir+"TMVA"+ftag+".root";
      cout<<"Correlations for "<<fn<<endl;
      if(TFile::Open(fn)==NULL)continue;
      correlations(fn);
      for(auto letter:vector<string>{"S","B"}){
	system(("epstopdf plots/CorrelationMatrix"+letter+".eps").c_str());
	system(("mv plots/CorrelationMatrix"+letter+".pdf "+pdir+"CorrelationMatrix"+letter+ftag+".pdf").c_str());
      }
    }
  }
  system("rm -rf plots/");
  return 0;
}
