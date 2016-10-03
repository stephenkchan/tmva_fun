#include "std_plotter.h"

std_plotter::std_plotter(const string&in,const string&out):std_ranking(in,out){
  TStyle *style=set_style();style->cd();
  c=new TCanvas(region_str(),region_str(),500,300);
}

// std_plotter::~std_plotter(){delete c;}

vector<TString> std_plotter::all_plots(const TString&pdir,int tDF){
  vector<TString>tmva_files,xls,yls;double totsqr=0.;vector<vector<double> >sigmtrx(ptvs().size(),vector<double>(jets().size(),0.));
  for(int nj=0;nj<(int)jets().size();nj++)yls.push_back(Form("%i jet",jets()[nj]));
  for(auto ptv:ptvs()){
    xls.push_back(ptv?"p_{T}^{V}<120 GeV":"p_{T}^{V}>120 GeV");
    for(int nj=0;nj<(int)jets().size();nj++){
      pair<TString,double> castle=single_plots(2,jets()[nj],ptv,pdir,tDF);
      tmva_files.push_back(castle.first);
      totsqr+=castle.second*castle.second;sigmtrx[!ptv][nj]=castle.second;
    }
  }
  totsqr=sqrt(totsqr);ostringstream nelson;nelson<<"S/#sqrt{S+B}_{tot}="<<setprecision(3)<<totsqr;
  print_info_plot_2d(sigmtrx,xls,yls,"std_sig_mtrx",pdir,TString(nelson.str()));
  return tmva_files;
}

pair<TString,double> std_plotter::single_plots(int nt,int nj,bool ptv,const TString&pdir,int tDF){
  set_njet(nj);set_ntag(nt);set_ptv(ptv);set_tag(string(region_str())+tdf_str(tDF));
  TString dir=check_mkdir(pdir);
  if(opendir(dir.Data())==NULL)system(("mkdir -p "+string(dir.Data())).c_str());
  cout<<"looking for "<<log(region_str().Data())<<endl;
  if(!file_exists(log(region_str().Data())))single_rank(true,tDF);
  vector<TString>vars=print_training_plot(c,region_str(),dir,tDF);
  vector<int>bgs{1,2,-1};double sig=0.;
  for(auto i:bgs)sig=bdt_testing(vars,false,(bdt_base)*this).print_test_bdts(c,dir,i,tDF).second;
  return {filename(vars),sig};
}

