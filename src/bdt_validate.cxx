#include "bdt_validate.h"
bdt_validate::bdt_validate(const string&in,const string&out,const string&tag,const vector<TString>&vs,bool overwrite,int nt,int nj,int eff_cut,int ratio,bool cts,bool lptv,bool wmll,bool met):
  bdt_trainer(in,out,tag,vs,overwrite,nt,nj,eff_cut,ratio,cts,lptv,wmll,met){}

bdt_validate::bdt_validate(const vector<TString>&vs,bool overwrite,const bdt_base&base):
  bdt_trainer(vs,overwrite,base){}

bdt_validate::bdt_validate(const bdt_trainer&train): bdt_trainer(train){}

pair<TH1D*,TH1D*> bdt_validate::get_validation_bdt_dists(int bg,int trans_DF){
  TFile*f=TFile::Open(root_file());
  TDirectoryFile*df=(TDirectoryFile*)f->Get("Method_BDT/BDT");
  if(df==NULL)return {NULL,NULL};
  TH1D*sig=(TH1D*)df->Get("MVA_BDT_S"),*big=(TH1D*)df->Get("MVA_BDT_B");
  pair<double,double>sb(evts_S_B(bg));
  if(debug)cerr<<"sb.first="<<sb.first<<", sb.second="<<sb.second<<endl;
  if(trans_DF<=0){
    sig->Scale(sb.first/sig->Integral());big->Scale(sb.second/big->Integral());
    return {sig,big};
  }
  pair<TH1D*,TH1D*>valid=transformation_DF(start.first,start.second,trans_DF>1);
  valid.first->Reset();valid.second->Reset();
  for(int ibin=1;ibin<sig->GetNbinsX();ibin++){
    valid.first->Fill(sig->GetBinCenter(ibin),sig->GetBinContent(ibin));
    valid.second->Fill(big->GetBinCenter(ibin),big->GetBinContent(ibin));
  }
  valid.first->Scale(sb.first/valid.first->Integral());
  valid.second->Scale(sb.second/valid.second->Integral());
  opfiles.push_back(f);
  return valid;
}

