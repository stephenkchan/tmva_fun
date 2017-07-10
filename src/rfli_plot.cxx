#include "rfli_plot.h"

rfli_plot::rfli_plot(const string&in,const string&out,int ptv,int nj,bool btvs):rfli_rank(in,out,ptv,nj,btvs){}

vector<TString> rfli_plot::all_plots(const TString&pdir,int tDF){
  //to do: code in a check to see if some log file exists with the output of a single_plots() call
  vector<TString>tmva_files,xls,yls={"combined"};vector<vector<double> >sigmtrx(cases().size(),vector<double>(jets().size()*ptvs().size()+1,0.));
  for(auto whic:cases())xls.push_back(var_tag(whic));
  for(auto j:jets()){
    for(auto lpt:ptvs()) yls.push_back(Form("%ijet,%s",j,ptv_label(lpt).Data()));
  }
  for(int i=0;i<(int)jets().size();i++){
    for(int k=0;k<(int)cases().size();k++){
      for(int j=0;j<(int)ptvs().size();j++){
	pair<TString,double> wompa=single_plots_log_check(cases()[k],jets()[i],ptvs()[j],pdir,tDF);
	tmva_files.push_back(wompa.first);sigmtrx[k][i*ptvs().size()+j+1]=wompa.second;sigmtrx[k][0]+=wompa.second*wompa.second;
	cout<<jets()[i]<<" jet, "<<ptv_label(ptvs()[j])<<", "<<var_tag(cases()[k])<<", sigmtrx["<<k<<"]["<<i*ptvs().size()+j+1<<"]="<<wompa.second<<"; sigmtrx["<<k<<"][0]="<<sigmtrx[k][0]<<endl;
      }
    }
  }
  for(int k=0;k<(int)cases().size();k++)sigmtrx[k][0]=sqrt(sigmtrx[k][0]);
  print_info_plot_2d(sigmtrx,xls,yls,"rfli_sig_mtrx"+tdf_str(tDF),pdir,"RFLI S/#sqrt{S+B} Summary");
  ofstream fout(pdir+"rfli_tmva_files"+tdf_str(tDF)+".txt");
  for(const auto&f:tmva_files)fout<<f<<endl;
  fout.close();
  return tmva_files;
}

pair<TString,double> rfli_plot::single_plots_log_check(int whic,int nj,int ptv,const TString&pdir,int tDF,bool split){
  set_pars(whic,nj,ptv);set_tftag(tDF);
  std::string logname=std::string(pdir.Data())+"test_sig_"+get_identifier()+".txt";
  if(file_exists(logname)){
    TString f;double s;ifstream fin(logname.c_str());
    fin>>f>>s;
    return {f,s};
  }
  else std::cout<<"Cannot find log "<<logname<<"...making plots..."<<std::endl;
  return single_plots(whic,nj,ptv,pdir,tDF,split);
}

pair<TString,double> rfli_plot::single_plots(int whic,int nj,int ptv,const TString&pdir,int tDF,bool split){
  TStyle *style=set_style();style->cd();set_pars(whic,nj,ptv);set_tag(tag()+tdf_str(tDF));
  TCanvas *c=new TCanvas(pdir+tag(),pdir+tag(),500,300);
  if(opendir(pdir.Data())==NULL)system(("mkdir -p "+string(pdir.Data())).c_str());
  if(!file_exists(var_log(tDF)))single_rank(whic,nj,ptv);
  vector<TString>vars=print_training_plot(c,tag(),pdir,tDF);
  double sig=0.;vector<int>bgs{-1};//{1,2,-1};
  bdt_base test_shell(*this);test_shell.set_tag(var_tag(which));
//   for(auto i:bgs)bdt_testing(tmva_in_dir,tmva_out_dir,var_tag(which),vars,false,ntag,njet,beff_cut,charm_ratio,cts_mv2,ptvbin,use_btagvs,widemll,metcut).print_test_bdts(c,pdir,i,tDF,split).second;
  for(auto i:bgs)bdt_testing(vars,false,test_shell).print_test_bdts(c,pdir,i,tDF,split).second;
  delete c;
  std::string logname=std::string(pdir.Data())+"test_sig_"+get_identifier()+".txt";
  ofstream fout(logname.c_str());
  fout<<filename(vars)<<" "<<sig<<std::endl;
  fout.close();
  return {filename(vars),sig};
}

