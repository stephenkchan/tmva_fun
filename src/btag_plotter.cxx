#include "btag_plotter.h"

btag_plotter::btag_plotter(const string&top):btag_ranking(top){
  TStyle *style=set_style();style->cd();
  c=new TCanvas(region_str(),region_str(),500,300);
}

// btag_plotter::~btag_plotter(){delete c;}

vector<TString> btag_plotter::all_plots(const TString&pdir,int tDF){
  //large matrix is variables (pTV(2) x jet(2)) x btagging cuts (4)
  //one for each (cts/binned) x (charm ratio) examined; 
  //summary plot with combined significances 
  vector<TString>tmva_files,xls,yls={"combined"},summary_labels;      vector<vector<double> >sigmtrxquad(ratios().size()*2,vector<double>(cuts().size(),0.));
  //efficiency cut labels
  for(int jc=0;jc<(int)cuts().size();jc++){xls.push_back(Form("%i%%",cuts()[jc]));}
  //phase space labels
  for(int mnj=0;mnj<(int)jets().size();mnj++){
    for(int lpt=0;lpt<(int)ptvs().size();lpt++)yls.push_back(Form("%ijet,%s-p^{V}_{T}",jets()[mnj],(ptvs()[lpt]?string("lo").c_str():string("hi").c_str())));
  }

  vector<vector<double> >big_sigmtrx(cuts().size(),vector<double>(cts().size()*ratios().size(),0.));
  for(int kct=0;kct<(int)cts().size();kct++){
    for(int ir=0;ir<(int)ratios().size();ir++){
      vector<vector<double> >sigmtrx(ptvs().size()*jets().size()+1,vector<double>(cuts().size(),0.));
      for(int jc=0;jc<(int)cuts().size();jc++){
	for(int mnj=0;mnj<(int)jets().size();mnj++){
	  for(int lpt=0;lpt<(int)ptvs().size();lpt++){
	    pair<TString,double> castle=single_plots(cuts()[jc],ratios()[ir],cts()[kct],jets()[mnj],ptvs()[lpt],pdir,tDF);
	    if(ratios()[ir]==20)tmva_files.push_back(castle.first);
	    sigmtrx[mnj*ptvs().size()+lpt+1][jc]=castle.second;sigmtrx[0][jc]+=castle.second*castle.second;
	    cout<<region_label()<<"; sigmtrx["<<jc<<"]["<<mnj*ptvs().size()+lpt+1<<"]="<<castle.second<<endl;
	  }
	}
	sigmtrx[0][jc]=sqrt(sigmtrx[0][jc]);//once all the significances have been squared and summed, take the sqrt for the total bin
	big_sigmtrx[jc][kct*ratios().size()+ir]=sigmtrx[0][jc];
	cout<<"big_sigmtrx["<<kct*ratios().size()+ir<<"]["<<jc<<"]="<<sigmtrx[0][jc]<<endl;
      }
      summary_labels.push_back(btag_var());
      print_info_plot_2d(sigmtrx,yls,xls,TString("btag_sig_mtrx-")+btag_var()+tdf_tstr(tDF),pdir,"S/#sqrt{S+B} summary: "+btag_var());
    }
  }
  if(big_sigmtrx.front().size()>1)print_info_plot_2d(big_sigmtrx,xls,summary_labels,TString("btag_sig_mtrx_tot")+tdf_tstr(tDF),pdir,"total btag S/#sqrt{S+B} summary");
  return tmva_files;
}

pair<TString,double> btag_plotter::single_plots(int cut,int ratio,bool cts,int nj,bool lptv,const TString&pdir,int tDF){
  set_btag(cut,ratio,cts);set_njet(nj);set_ptv(lptv);set_tag(string(btag_str()+region_str())+tdf_str(tDF));
  cout<<"Processing "<<btag_label()<<" "<<region_label()<<endl;
//   cout<<cut<<","<<ratio<<","<<cts<<" "<<these_cuts()<<endl;
  TString dir=check_mkdir(pdir);
  if(debug)cout<<"does this exist? "<<log(string(btag_str()+region_str()))<<endl;
  if(!file_exists(log(string(btag_str()+region_str()))))single_rank(true,tDF);
  vector<TString>vars=print_training_plot(c,btag_str()+region_str(),dir,tDF);
  double sig=0.;vector<int>bgs{-1};//{1,2,
  for(auto i:bgs)    sig=bdt_testing(vars,false,(bdt_base)*this).print_test_bdts(c,dir,i,tDF).second;
  return {filename(vars),sig};
}

