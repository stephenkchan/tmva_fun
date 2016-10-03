#include "li_ranking.h"

int main(int argc,char**argv){
  string tmva_in="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20160307_srli/",tmva_out=tmva_in+"output/";
  int which=-99,nj=2,tF=2;bool lptv=true,overwrite=false,btv=true;
  if(argc>1)  tmva_in=argv[1];
  if(argc>2) tmva_out=argv[2];
  if(argc>3)    which=atoi(argv[3]);
  if(argc>4)       nj=atoi(argv[4]);
  if(argc>5)     lptv=atoi(argv[5]);
  if(argc>6)      btv=atoi(argv[6]);
  if(argc>7)overwrite=atoi(argv[7]);
  if(argc>8)       tF=atoi(argv[8]);
  li_ranking odin(tmva_in,tmva_out,lptv,nj,btv);
  if(which==-99)odin.all_ranks(overwrite,tF);
  else odin.finish_rank(which,nj,lptv,overwrite,tF);
  return 0;
}
