#include "li_plotter.h"

li_plotter::li_plotter(const string&in,const string&out,bool lptv,int nj,bool btvs):li_ranking(in,out,lptv,nj,btvs){}

vector<TString> li_plotter::all_plots(const TString&pdir,int tDF){
  vector<TString>tmva_files,xls,yls={"combined"};vector<vector<double> >sigmtrx(cases().size(),vector<double>(jets().size()*ptvs().size()+1,0.));
  for(auto whic:cases())xls.push_back(var_tag(whic));
  for(auto j:jets()){
    for(auto lpt:ptvs()) yls.push_back(Form("%ijet,%s-p^{V}_{T}",j,(lpt?string("lo").c_str():string("hi").c_str())));
  }
  for(int i=0;i<(int)jets().size();i++){
    for(int k=0;k<(int)cases().size();k++){
      for(int j=0;j<(int)ptvs().size();j++){
	pair<TString,double> wompa=single_plots(cases()[k],jets()[i],ptvs()[j],pdir,tDF);
	tmva_files.push_back(wompa.first);sigmtrx[k][i*ptvs().size()+j+1]=wompa.second;sigmtrx[k][0]+=wompa.second*wompa.second;
	cout<<jets()[i]<<" jet, "<<(ptvs()[j]?"LO":"HI")<<" pTV, "<<var_tag(cases()[k])<<", sigmtrx["<<k<<"]["<<i*ptvs().size()+j+1<<"]="<<wompa.second<<"; sigmtrx["<<k<<"][0]="<<sigmtrx[k][0]<<endl;
      }
    }
  }
  for(int k=0;k<(int)cases().size();k++)sigmtrx[k][0]=sqrt(sigmtrx[k][0]);
  print_info_plot_2d(sigmtrx,xls,yls,"li_sig_mtrx"+tdf_str(tDF),pdir,"LI S/#sqrt{S+B} Summary");
  return tmva_files;
}

pair<TString,double> li_plotter::single_plots(int whic,int nj,bool lpt,const TString&pdir,int tDF,bool split){
  TStyle *style=set_style();style->cd();set_pars(whic,nj,lpt);set_tag(tag()+tdf_str(tDF));
  TCanvas *c=new TCanvas(pdir+tag(),pdir+tag(),500,300);
  if(opendir(pdir.Data())==NULL)system(("mkdir -p "+string(pdir.Data())).c_str());
  if(!file_exists(var_log(tDF)))single_rank(whic,nj,lpt);
  vector<TString>vars=print_training_plot(c,tag(),pdir,tDF);
  double sig=0.;vector<int>bgs{-1};//{1,2,-1};
  for(auto i:bgs)sig=bdt_testing(vars,false,(bdt_base)*this).print_test_bdts(c,pdir,i,tDF,split).second;
  delete c;
  return {filename(vars),sig};
}

