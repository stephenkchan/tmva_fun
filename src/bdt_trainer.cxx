#include "bdt_trainer.h"

bdt_trainer::bdt_trainer(const string&in,const string&out,const string&tag,const vector<TString>&vs,bool overwrite,int nt,int nj,int eff_cut,int ratio,bool cts,int lptv,bool btv,bool wmll,bool met,bool do_validation,int seed)
  :bdt_base(in,out,tag,nt,nj,eff_cut,ratio,cts,lptv,btv,wmll,met),vars(vs),validate(do_validation),kseed(seed){initialize(overwrite,do_validation,seed);}

bdt_trainer::bdt_trainer(const vector<TString>&vs,bool overwrite,const bdt_base&base,bool do_validation,int seed)
  :bdt_base(base),vars(vs),validate(do_validation),kseed(seed){initialize(overwrite,do_validation,seed);}

void bdt_trainer::initialize(bool overwrite,bool do_validation,int seed){
  if(opendir(tmva_in_dir.c_str())==NULL){
    cout<<"TMVA data dir "<<tmva_in_dir<<" does not exist. Exiting..."<<endl;
    exit(2);
  }
  if(debug)cout<<"file: "<<root_file()<<"; cuts: "<<these_cuts()<<endl;
  opfiles=vector<TFile*>();
  bdt_dists(overwrite,do_validation,seed);
}
bdt_trainer::~bdt_trainer(){
  for(auto f:opfiles)f->Close();
}
pair<double,double> bdt_trainer::significance(){
  if(debug){
    pair<TH1D*,TH1D*>powerup=transformation_DF(start,false),powedup=transformation_DF(start,true);
    cerr<<"train std: "<<optimal_sig(start).second<<" t_F: "<<optimal_sig(powerup).second<<" t_D: "<<optimal_sig(powedup).second<<endl;
  }
  return optimal_sig(start);
}

void bdt_trainer::bdt_dists(bool overwrite,bool do_validation,int seed,int bg){
  DIR *zoom=opendir(tmva_out_dir.c_str());struct dirent *file;
  if(vars.empty()){cout<<"derp"<<endl;return;}
  bool file_found=false;string target(root_file(bg).Data()),wudder=target.substr(target.find_last_of("/")+1);
  while((file=readdir(zoom))!=NULL){
    if(string(file->d_name)==wudder)file_found=true;
  }
  if(!file_found||overwrite)train(bg,do_validation,seed);
  TFile*f=TFile::Open(root_file(bg));opfiles.push_back(f);
  TDirectoryFile*df=(TDirectoryFile*)f->Get("Method_BDT/BDT");
  if(df==NULL)return;
  TH1D*sig=(TH1D*)df->Get("MVA_BDT_Train_S"),*big=(TH1D*)df->Get("MVA_BDT_Train_B");pair<double,double>sb(evts_S_B(-99));
  double snorm=sb.first/sig->Integral(),bnorm=sb.second/big->Integral();
  sig->Scale(snorm);big->Scale(bnorm);opfiles.push_back(f);
//   if(debug)cerr<<"S events: "<<sig->Integral()<<" B events: "<<big->Integral()<<endl;
  if(debug)cout<<"sig big="<<sig<<" " <<big<<endl;
  start={sig,big};
}

pair<double,double> bdt_trainer::evts_S_B(int bg)const{
  TFile*f=TFile::Open(root_file(-1));
  TH1D*events=(TH1D*)f->Get("events");
  double sigev=events->GetBinContent(1),bigev;
  if(bg==1)bigev=events->GetBinContent(4+1,-1);//Z+jets; back entries starting at 4
  else if(is_single_bg(bg))bigev=events->GetBinContent(bg+1);//single backgrounds
  else bigev=events->Integral(2,-1);//all background
  delete events;f->Close();
  return {sigev,bigev};
}

void bdt_trainer::train(int bg,bool do_validation,int seed)const{
  TString outfileName=root_file(bg,do_validation,seed);
  TString mycuts = these_cuts();
  if(do_validation) mycuts+=" && EventNumber%3!=0";//the last to set aside a third of events for validation
  if(debug)cout<<"Begin training to be saved in file "<<outfileName<<", cuts: "<<these_cuts()<<endl;
  // This loads the library
  TMVA::Tools::Instance();

  // --- Here the preparation phase begins
  std::vector<TFile*> openfiles;

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( xml_tag(do_validation,seed), outputFile,"!V:Silent:Color:DrawProgressBar:AnalysisType=Classification" );
  for(auto var:vars)factory->AddVariable(var,'F');
  factory->AddSpectator("EventNumber",'I');

  TH1D* events=new TH1D("events","events",n_sam(),-0.5,n_sam()-0.5);
  // global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0, backgroundWeight = 1.0;
  TSystemDirectory dir(tmva_in_dir.c_str(),tmva_in_dir.c_str());
  TList* files = dir.GetListOfFiles(); TFile* inputFile;
  pair<double,double>sbg;
  //test cut: if our seed is even, we want the test_cut to be EventNumber%2==0, if odd EventNumber%2==1 and vice versa for the train_cut
  int test_int=(seed%2)!=0,train_int=(seed%2)==0;
  TString cut=these_cuts(),test_cut=Form("%s && EventNumber%2==%i",mycuts.Data(),test_int),train_cut=Form("%s && EventNumber%2==%i",mycuts.Data(),train_int);bool is_manual=kseed>0;
  for (TObject* obj : *files) {
    TString fileName = obj->GetName();
    int sam=which_sam(fileName);
    //if the sample number is less than zero, skip
    if(!process_sample(bg,sam))continue;
    fileName.Prepend(dir.GetName());
    inputFile = TFile::Open(fileName);
    openfiles.push_back(inputFile);
    TTree* tree = (TTree*) inputFile->Get("Nominal");
    if (tree==NULL)continue;
    if (tree->GetEntries() == 0) continue;
    double event=tree->Draw("EventWeight>>hist",cut);
    TH1F *hist = (TH1F*)gDirectory->Get("hist");  event*=hist->GetMean();
    events->Fill(sam,event);
    if(is_manual){
      factory->AddTree(tree,(sam?"Background":"Signal"),(sam?backgroundWeight:signalWeight),TCut(train_cut),"train");
      factory->AddTree(tree,(sam?"Background":"Signal"),(sam?backgroundWeight:signalWeight), TCut(test_cut), "test");
    }
    else{
      if(sam)factory->AddBackgroundTree(tree, backgroundWeight);
      else factory->AddSignalTree(tree, signalWeight);
    }
  }
  factory->SetWeightExpression("EventWeight");
  TString test_train_config="SplitMode=Alternate:NormMode=EqualNumEvents";
  if(seed!=0)test_train_config=Form("SplitMode=Alternate:NormMode=EqualNumEvents:SplitSeed=%i",seed);
  factory->PrepareTrainingAndTestTree(TCut(mycuts),test_train_config);

  // ---- Book MVA methods
  factory->BookMethod( TMVA::Types::kBDT,"BDT","!H:!V:NTrees=200:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.15:SeparationType=GiniIndex:nCuts=100:PruneMethod=NoPruning:MinNodeSize=1.5%" );
  
  // ---- Now you can tell the factory to train, test, and evaluate the MVAs
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  // Save the output
  outputFile->cd();
  events->Write();
  outputFile->Close();

  std::cout << "==> TMVAClassification is done!" << std::endl;
  for (TFile* f : openfiles) f->Close();
  delete factory;
  string weightsfile="./weights/";weightsfile+=xml_file(do_validation,seed);
  string cmd="mv "+weightsfile+" "+tmva_out_dir;
  system(cmd.c_str());
}

pair<double,double> bdt_trainer::print_bdt_plots(TCanvas *c,pair<TH1D*,TH1D*>iris,const TString&pdir,const TString&tag,const TString&bglabel,bool funky_style,double scut){
  TStyle *style=set_style();style->cd();  c->cd();
  pair<TH1D*,TH1D*> bleuet={(TH1D*)start.first->Clone(),(TH1D*)start.second->Clone()};
  bool printcut=(scut>=iris.first->GetXaxis()->GetXmin()&&scut<=iris.first->GetXaxis()->GetXmax());
  pair<double,double>hyacinth=(printcut?cut_sig(iris,scut):cumulat_sig(iris));
  int bone=12,sone=8,btwo=0,stwo=20;
  bleuet.first->Scale(1./bleuet.first->Integral());
  bleuet.second->Scale(1./bleuet.second->Integral());
  iris.first->Scale(1./iris.first->Integral()*iris.first->GetNbinsX()/bleuet.first->GetNbinsX());
  iris.second->Scale(1./iris.second->Integral()*iris.first->GetNbinsX()/bleuet.first->GetNbinsX());
  iris.first->GetXaxis()->SetLabelSize(0.03);
  iris.first->GetXaxis()->SetTitleSize(0.04);
  iris.first->GetXaxis()->SetTitleOffset(0.85);
  iris.first->GetXaxis()->SetTitle("BDT score");
  iris.first->GetYaxis()->SetLabelSize(0.03);
  iris.first->GetYaxis()->SetTitleSize(0.04);
  iris.first->GetYaxis()->SetTitleOffset(0.85);
  iris.first->GetYaxis()->SetTitle("a.u.");
//   vector<TH1D*>bouquet={iris.first,iris.second,bleuet.first,bleuet.second};
//   auto_limits(bouquet);
  bleuet.first->SetMarkerColor(color(stwo));bleuet.first->SetLineColor(color(stwo));bleuet.first->Draw("E");
  bleuet.second->SetMarkerColor(color(btwo));bleuet.second->SetLineColor(color(btwo));bleuet.second->Draw("E same");
  if(funky_style){iris.first->SetFillStyle(3001);iris.second->SetFillStyle(3001);}
  if(debug){
    cout<<"The hists are bleuet: "<<bleuet.first<<","<<bleuet.second<<" | iris: "<<iris.first<<","<<iris.second<<endl;
    cout<<"The command in question is ("<<iris.first<<")->SetFillColorAlpha("<<color(sone)<<","<<"0.75)"<<endl;
  }
  iris.first->SetFillColorAlpha(color(sone),0.75);
//   iris.first->SetFillColor(color(sone));
  iris.first->SetLineColor(kBlack);iris.first->Draw("hist same");
  iris.second->SetFillColorAlpha(color(bone),0.75);
  iris.second->SetLineColor(kBlack);iris.second->Draw("hist same");
  double lo=style->GetPadLeftMargin(),hi=1.-style->GetPadRightMargin();
  double xmin=iris.first->GetXaxis()->GetXmin(),xmax=iris.first->GetXaxis()->GetXmax(),cut_val=printcut?hyacinth.first:0.35,line_x=lo+(hi-lo)*(cut_val-xmin)/(xmax-xmin);
  if(printcut){
    TLine*line=new TLine();  line->SetLineColor(kSpring-2);  line->SetLineWidth(3);
    line->DrawLineNDC(line_x,style->GetPadBottomMargin(),line_x,1.-style->GetPadTopMargin());
  }
  TLatex l=plot_latex();l.SetTextSize(0.05);double flip_flop=0.67,cut_x=(line_x<flip_flop?line_x+0.01:line_x-0.35),bg_x=(line_x<flip_flop?line_x+0.06:line_x-0.15),bg_x_long=(line_x<flip_flop?line_x+0.01:line_x-0.30);
  if(bglabel.Length()<1.05*tag.Length() && bglabel.Length()>0.5*tag.Length())bg_x=cut_x;
  else if(bglabel.Length()>1.05*tag.Length())bg_x=bg_x_long;
  ostringstream cut;
  if(printcut)cut<<"Cut="<<setprecision(3)<<cut_val<<", sig="<<hyacinth.second<<(funky_style?" (VAL)":" (TEST)");
  else cut<<"Sig="<<hyacinth.second<<(funky_style?" (VAL)":" (TEST)");
  l.DrawLatex(cut_x,0.8,cut.str().c_str());
  l.DrawLatex(bg_x,0.74,bglabel);
  if(iris.first->GetNbinsX()==bleuet.first->GetNbinsX()){
    double kolS=iris.first->KolmogorovTest(bleuet.first),kolB=iris.second->KolmogorovTest(bleuet.second);
    TString probatext = Form( "Kolmogorov-Smirnov test: signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
    l.DrawLatex(0.10,1.-style->GetPadTopMargin()*0.9,probatext);
  }
  c->Print(check_dir(pdir)+"bdt_dist-"+tag+".pdf");delete bleuet.first;delete bleuet.second;
  return hyacinth;
}

