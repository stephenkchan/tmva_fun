#include "std_ranking.h"

int main(int argc,char**argv){
  string tmva_in="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20160307_srli/",tmva_out=tmva_in+"output/";
  int nt=-99,nj=2,tF=2;bool lptv=true;bool overwrite=false;
  if(argc>1)  tmva_in=argv[1];
  if(argc>2) tmva_out=argv[2];
  if(argc>3)       nt=atoi(argv[3]);
  if(argc>4)       nj=atoi(argv[4]);
  if(argc>5)     lptv=atoi(argv[5]);
  if(argc>6)overwrite=atoi(argv[6]);
  if(argc>7)       tF=atoi(argv[7]);
  std_ranking odin(tmva_in,tmva_out,nt,nj,lptv);
  if(nt==-99)odin.all_ranks(overwrite,tF);
  else odin.finish_rank(overwrite,tF);
  return 0;
}
