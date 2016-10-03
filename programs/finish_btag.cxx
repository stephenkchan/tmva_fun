#include "btag_ranking.h"

int main(int argc,char**argv){
  string tmva_top="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20160307_srli/";
  int cut=-99,ratio=20,nj=2,tF=2;bool cts=true,lptv=true;bool overwrite=false;
  if(argc>1) tmva_top=argv[1];
  if(argc>2)      cut=atoi(argv[2]);
  if(argc>3)    ratio=atoi(argv[3]);
  if(argc>4)      cts=atoi(argv[4]);
  if(argc>5)       nj=atoi(argv[5]);
  if(argc>6)     lptv=atoi(argv[6]);
  if(argc>7)overwrite=atoi(argv[7]);
  if(argc>8)       tF=atoi(argv[8]);
  // cout<<"initialize with..."<<tmva_top<<","<<cut<<","<<ratio<<","<<cts<<","<<nj<<","<<lptv<<endl;
  btag_ranking odin(tmva_top,cut,ratio,cts,nj,lptv);
  if(cut==-99)odin.all_ranks(overwrite,tF);
  else odin.finish_rank(overwrite,tF);
  return 0;
}
