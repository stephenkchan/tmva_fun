#include "rfli_plot.h"

rfli_plot::rfli_plot(const string&in,const string&out,int ptv,int nj,bool btvs,int whic):rfli_rank(in,out,ptv,nj,btvs,whic){debug=true;}

vector<TString> rfli_plot::all_plots(const TString&pdir,int tDF){
  //to do: code in a check to see if some log file exists with the output of a single_plots() call
  vector<TString>xls,yls={"combined"};vector<vector<double> >sigmtrx(cases().size(),vector<double>(jets().size()*ptvs().size()+1,0.)),errmtrx(cases().size(),vector<double>(jets().size()*ptvs().size()+1,0.));
  for(auto whic:cases())xls.push_back(var_tag(whic));
  for(auto j:jets()){
    for(auto lpt:ptvs()) yls.push_back(Form("%i%sjet,%s p_{T}^{V}",j,j==3?"+":"",lpt>1?"hi":"lo"));
  }
  for(int i=0;i<(int)jets().size();i++){
    for(int k=0;k<(int)cases().size();k++){
      for(int j=0;j<(int)ptvs().size();j++){
	pair<double,double> wompa=single_plots_log_check(cases()[k],jets()[i],ptvs()[j],pdir,tDF);
	sigmtrx[k][i*ptvs().size()+j+1]=wompa.first; //sigmtrx[k][0]+=wompa.first*wompa.first;
	errmtrx[k][i*ptvs().size()+j+1]=wompa.second;//errmtrx[k][0]+=wompa.second*wompa.second;
// 	cout<<jets()[i]<<" jet, "<<ptv_label(ptvs()[j])<<", "<<var_tag(cases()[k])<<", sigmtrx["<<k<<"]["<<i*ptvs().size()+j+1<<"]="<<wompa.second<<"; sigmtrx["<<k<<"][0]="<<sigmtrx[k][0]<<endl;
      }
    }
  }
  for(int k=0;k<(int)cases().size();k++)quad_sum_error(sigmtrx[k],errmtrx[k]);
  print_info_plot_2d(sigmtrx,xls,yls,"rfli_sig_mtrx"+tdf_str(tDF),pdir,"RFLI S/#sqrt{S+B} Summary",errmtrx);
  /*
  ofstream fout(pdir+"rfli_tmva_files"+tdf_str(tDF)+".txt");
  for(const auto&f:tmva_files)fout<<f<<endl;
  fout.close();
  */
  return vector<TString>(1,"don't use this value");
}

void rfli_plot::quad_sum_error(vector<double>&sigs,vector<double>&errors)const{//assume the format is the one above
  //|d/dx sqrt(x^2+blah)| = x/sqrt(x^2+blah)
  int nent=sigs.size();
  if(sigs.size()!=errors.size())return;
  double sigtot=0.,errtot=0.;
  for(int i=1;i<nent;i++)sigtot += sigs[i]*sigs[i];
  sigs[0]=sqrt(sigtot);
  for(int i=1;i<nent;i++)errtot += sigs[i]/sigs[0]*errors[i]*errors[i];
  errors[0]=sqrt(errtot);
}

pair<double,double> rfli_plot::single_plots_log_check(int whic,int nj,int ptv,const TString&pdir,int tDF,bool split){
  set_pars(whic,nj,ptv);set_tftag(tDF);
  std::string logname=std::string(pdir.Data())+"test_sig_"+get_identifier()+".txt";
  if(file_exists(logname)){
    TString f;double s,e;ifstream fin(logname.c_str());
    fin>>f>>s>>e;
    return {s,e};
  }
  else std::cout<<"Cannot find log "<<logname<<"...making plots..."<<std::endl;
  return single_plots(whic,nj,ptv,pdir,tDF,split);
}

pair<double,double> rfli_plot::single_plots(int whic,int nj,int ptv,const TString&pdir,int tDF,bool split){
  TStyle *style=set_style();style->cd();set_pars(whic,nj,ptv);set_tag(tag()+tdf_str(tDF));
  TCanvas *c=new TCanvas(pdir+tag(),pdir+tag(),500,300);
  if(opendir(pdir.Data())==NULL)system(("mkdir -p "+string(pdir.Data())).c_str());
  if(!file_exists(var_log(tDF)))single_rank(whic,nj,ptv);
  vector<TString>vars=print_training_plot(c,tag(),pdir,tDF);
  double sig=0.,err=0.;vector<int>bgs{-1};//{1,2,-1};
  bdt_base test_shell(*this);test_shell.set_tag(var_tag(which));
//   for(auto i:bgs)bdt_testing(tmva_in_dir,tmva_out_dir,var_tag(which),vars,false,ntag,njet,beff_cut,charm_ratio,cts_mv2,ptvbin,use_btagvs,widemll,metcut).print_test_bdts(c,pdir,i,tDF,split).second;
  for(auto i:bgs){
    pair<double,double>es=bdt_testing(vars,false,test_shell).print_test_bdts(c,pdir,i,tDF,split);
    err=es.first;sig=es.second;
  }
  delete c;
  std::string logname=std::string(pdir.Data())+"test_sig_"+get_identifier()+".txt";
  ofstream fout(logname.c_str());
  fout<<filename(vars)<<" "<<sig<<" "<<err<<std::endl;
  fout.close();
  return {sig,err};
}

