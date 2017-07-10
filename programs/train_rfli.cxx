#include "rfli_rank.h"

int main(int argc,char**argv){
  //these are the lines we want to turn into an executable for the python
  /*
    if(buffer.find(var)!=buffer.end())continue;
    vector<TString>frodo(done);frodo.push_back(var);
  */
  // so we want to take an rfli_rank object and cast down to a bdt_validate
  // a lot of these arguments will be inefficient to use this way
  // but construction isn't so expensive
  if(argc<10){
    std::cout<<"You're doing this wrong.  Don't use this executable unless you know what you're doing."<<std::endl;
    return 0;
  }
  // construct the objects and do the training/significance calculation
  string tmva_in=argv[1],tmva_out=argv[2];
  int which=atoi(argv[3]),nj=atoi(argv[4]),lptv=atoi(argv[5]);
  bool btv=atoi(argv[6]),overwrite=atoi(argv[7]);
  int tF=atoi(argv[8]);std::vector<TString>frodo;
  for(int i=9;i<argc;i++){
    frodo.push_back(TString(argv[i]));
//     std::cout<<argv[i]<<std::endl;
  }
  rfli_rank odin(tmva_in,tmva_out,lptv,nj,btv,which);
  odin.set_tftag(tF);
  double sig=bdt_validate(frodo,overwrite,(bdt_base)odin).val_sig(-1,tF).second;

  // print the result to console but also save to a text file (this is a little dumb, but whatever)
  std::string logdir=tmva_out+"rank-"+odin.get_identifier()+"/",name=logdir; odin.check_mkdir(logdir);
  for(const auto&s:frodo)name+=(s+"-");
  name=name.substr(0,name.length()-1);
  name+=".txt";
  ofstream fout(name.c_str());
  std::cout<<frodo.back()<<","<<sig<<" ";
  fout<<sig<<std::endl;
  fout.close();
  return 0;
}
