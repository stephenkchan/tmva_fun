#include "rfli_rank.h"

int main(int argc,char**argv){
  string tmva_in="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/201704/",tmva_out=tmva_in+"output/";
  int which=-99,nj=2,tF=2,lptv=0;bool overwrite=false,btv=true;
  if(argc>1)  tmva_in=argv[1];
  if(argc>2) tmva_out=argv[2];
  if(argc>3)    which=atoi(argv[3]);
  if(argc>4)       nj=atoi(argv[4]);
  if(argc>5)     lptv=atoi(argv[5]);
  if(argc>6)      btv=atoi(argv[6]);
  if(argc>7)overwrite=atoi(argv[7]);
  if(argc>8)       tF=atoi(argv[8]);
  rfli_rank odin(tmva_in,tmva_out,lptv,nj,btv,which);
  if(which==-99)odin.all_ranks(overwrite,tF);
  else if(argc>9){
    if(std::string(argv[9])=="front-matter")odin.ranking_front_matter(odin.tag(),overwrite,tF);
    else{
      odin.set_tftag(tF);
      std::vector<TString>done;std::vector<double>daisy_chain;
      for(int i=9;i<argc;i++){
	std::string split=argv[i];
	done.push_back(TString(split.substr(0,split.find(","))));
	daisy_chain.push_back(std::stof(split.substr(split.find(",")+1)));
      }
      ifstream fin(odin.log().c_str());
      std::string cable;bool print=true;
      while(fin.good()){
	fin>>cable;
	if(cable=="DAISY_CHAIN:"){
	  std::cout<<"Results seem to exist in your ranking log "<<odin.log()<<"; won't double print"<<std::endl;
	  print=false;
	}
      }
      ofstream fout(odin.log().c_str(),std::ofstream::app);
      odin.set_tag(odin.tag());//setting tf tag adds the transformation * tag to the identifier; remove it once again so the kfold files don't have this information
      //for now, don't transform the kfold xml's...things are rebinned in the WSMaker
      odin.ranking_back_matter(done,daisy_chain,fout,overwrite,0,print);
    }
  }
  else{
    cout<<"Humans shouldn't be using this program directly.  You haven't specified the right kind of callable stuff.  You had used the command:"<<endl;
    for(int i=0;i<argc;i++)cout<<" "<<argv[i];
    cout<<endl;
  }
  return 0;
}
