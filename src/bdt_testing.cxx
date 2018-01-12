#include "bdt_testing.h"
bdt_testing::bdt_testing(const string&in,const string&out,const string&tag,const vector<TString>&vs,bool overwrite,int nt,int nj,int eff_cut,int ratio,bool cts,int lptv,bool btv,bool wmll,bool met):
  bdt_trainer(in,out,tag,vs,overwrite,nt,nj,eff_cut,ratio,cts,lptv,btv,wmll,met){}

bdt_testing::bdt_testing(const vector<TString>&vs,bool overwrite,const bdt_base&base):
  bdt_trainer(vs,overwrite,base){}

bdt_testing::bdt_testing(const bdt_trainer&train):bdt_trainer(train){}

pair<TH1D*,TH1D*> bdt_testing::make_testing_bdt_dists(int bg,int trans_DF,bool split)const{
  pair<double,double>sb(evts_S_B(bg));
  if(debug){
    cout<<filename(vars)<<" S="<<sb.first<<", B="<<sb.second<<" "<<start.first<<" "<<start.second<<endl;
    cout<<" "<<start.first->GetEntries()<<" "<<start.second->GetEntries()<<endl;
    cout<<(trans_DF>0)<<", and start's first has name "<<start.first->GetName()<<endl;
  }
  pair<TH1D*,TH1D*>valid=(trans_DF>0?transformation_DF(start,trans_DF>1):pair<TH1D*,TH1D*>((TH1D*)start.first->Clone(),(TH1D*)start.second->Clone()));
  valid.first->Reset();valid.second->Reset();
  TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );    
  map<TString,pair<float,TBranch*>>lisbon;
  for(auto nom:vars)lisbon[nom]=pair<float,TBranch*>(0.,NULL);
  for(auto nom:vars)reader->AddVariable(nom,&lisbon.at(nom).first);
  vector<TString>extra_vars={"EventWeight","mLL",btag_var(charm_ratio,cts_mv2)+"B1",btag_var(charm_ratio,cts_mv2)+"B2","pTV"};
  for(auto var:extra_vars){
    if(lisbon.find(var)==lisbon.end())lisbon[var]=pair<float,TBranch*>(0.,NULL);
  }
  TString xml=tmva_out_dir+xml_file(split);
  //if(debug)
  cout<<"Evaluation XML: "<<xml<<endl;
  reader->BookMVA("BDT",xml);TBranch*bev,*bnj;unsigned long long ev;int nj;

  TSystemDirectory dir(tmva_in_dir.c_str(),tmva_in_dir.c_str());  std::vector<TFile*> openfiles;
  TList* files = dir.GetListOfFiles(); TFile* inputFile;
  pair<double,double>sbg;TString cut=these_cuts();
  for (TObject* obj : *files) {
    TString fileName = obj->GetName();
    int sam=which_sam(fileName);
    if (!process_sample(bg,sam)) continue;
    if(debug)cout<<"Cracking open file "<<fileName<<endl;
    fileName.Prepend(dir.GetName());
    inputFile = TFile::Open(fileName);
    openfiles.push_back(inputFile);
    TTree* tree = (TTree*) inputFile->Get("Nominal");
    if (tree==NULL)continue;
    if (tree->GetEntries() == 0) continue;
    //set up branches
    for(auto&it:lisbon)tree->SetBranchAddress(it.first,&(it.second.first),&(it.second.second));
    tree->SetBranchAddress("EventNumber",&ev,&bev);tree->SetBranchAddress("nSigJet",&nj,&bnj);
    for(int iev=0;iev<tree->GetEntries();iev++){
//       cout<<"At "<<iev+1<<" of "<<tree->GetEntries()<<endl;
      for(auto&it:lisbon)it.second.second->GetEntry(iev);
      bev->GetEntry(iev);bnj->GetEntry(iev);
      /*
      if(passes_cut(lisbon,ev,nj)!=passes_cut(tree,ev,cut)){
	cout<<fileName<<"---Incompatible cut evaluation (branch="<<passes_cut(lisbon,ev,nj)<<",string="<<passes_cut(tree,ev,cut)<<"); the string cut was "<<TString(TString::Format("%s  && EventNumber==%llu",cut.Data(),ev)).Data()<<endl;
	cout<<"There would have been: "<<tree->Draw("",TString::Format("%s  && EventNumber==%llu",cut.Data(),ev))<<endl;
      }
      */
      if(!passes_cut(lisbon,ev,nj))continue;
//       if(!passes_cut(tree,ev,cut))continue;//this is actually very slow
      double val=reader->EvaluateMVA("BDT"),weight=lisbon.at("EventWeight").first;
      if(sam) valid.second->Fill(val,weight);
      else{
	valid.first->Fill(val,weight);
	/*
 	if(val>0.5||iev==92499){
	  cout<<iev<<" Fill "<<val<<" with weight "<<weight<<endl;
	  for(auto&it:lisbon)cout<<setw(10)<<it.first<<setprecision(3)<<setw(6)<<it.second.first;
	  cout<<endl;
	  }*/
      }
    }
  }
  for (TFile* f : openfiles) f->Close();
  delete reader;
  //test on only a third of events; scale up by a factor of three for proper estimates
  //actually, let's give it the "right" normalization; best apples to apples comparison
//   cout<<"Scale: "<<sb.first<<"/"<<valid.first->Integral()<<" "<<bgev<<"/"<<valid.second->Integral()<<endl;
  valid.first->Scale(sb.first/valid.first->Integral());
  valid.second->Scale(sb.second/valid.second->Integral());
  return valid;
}

bool bdt_testing::passes_cut(const map<TString,pair<float,TBranch*> >& vars,int ev,int nj)const{
  if(ev%3!=0)return false;
  if(nj!=njet&&njet>1)return false;
  if(vars.find("pTV")==vars.end()||vars.find(btag_var(charm_ratio,cts_mv2)+"B2")==vars.end()||vars.find(btag_var(charm_ratio,cts_mv2)+"B1")==vars.end()||vars.find("mLL")==vars.end())return false;
  int tags=2;
  if(vars.at(btag_var(charm_ratio,cts_mv2)+"B1").first<=mv2c_cut(beff_cut,charm_ratio))tags--;
  if(vars.at(btag_var(charm_ratio,cts_mv2)+"B2").first<=mv2c_cut(beff_cut,charm_ratio))tags--;
  if(tags<ntag)return false;
//   cout<<vars.at("pTV").first<<" "<<vars.at("MV2c20B2").first<<" "<<vars.at("MV2c20B1").first<<" "<<vars.at("mLL").first<<" "<<ev<<" "<<nj<<endl;
  double ptv=0.001*vars.at("pTV").first;
  if(ptv>75. && ptvbin==0)return false;//(0,75] GeV
  if((ptv<=75.||ptv>150.) && ptvbin==1)return false;//(75,150] GeV
  if(ptv<=150. && ptvbin==2)return false;//(150,\infty) GeV
  if(vars.at("mLL").first/1000.<=81 || vars.at("mLL").first/1000.>=101.)return false;
  return true;
}

bool bdt_testing::passes_cut(TTree*t,int ev,const TString& cut)const{
  if(debug)cout<<TString(TString::Format("%s  && EventNumber==%i",cut.Data(),ev)).Data()<<endl;
  if(ev%3!=0)return false;
  return t->Draw("",TString::Format("%s  && EventNumber==%i",cut.Data(),ev))>0;
}
