#include "bdt_base.h"

bdt_base::bdt_base(const string&in,const string&out,const string&tag,int nt,int nj,int eff_cut,int ratio,bool cts,int lptv,bool use_btvs,bool wmll,bool met)
  :debug(false),tmva_in_dir(check_dir(in)),tmva_out_dir(check_mkdir(out)),ntag(nt),njet(nj),beff_cut(eff_cut),charm_ratio(ratio),cts_mv2(cts),loptv(lptv),use_btagvs(use_btvs),widemll(wmll),metcut(met){
  //akt4EM FixedCutBEff_90,85,77,70,60,50,30 (numbers are nominal tag efficiencies)
  //60,70,77,85 are supported by the tool right now
  btag_ratio_cuts={ {20,{{90,-0.9185}, {85,-0.7887}, {80,-0.5911}, {77,-0.4434}, {70,-0.0436}, {60,0.4496}, {50,0.7535}, {30,0.9540}}} };
  def_cratio=btag_ratio_cuts.begin()->first;
  if(btag_ratio_cuts.begin()->second.find(70)==btag_ratio_cuts.begin()->second.end()){
    int dist=9999;
    for(auto&cutpair:btag_ratio_cuts.begin()->second){
      if(abs(70-cutpair.first)<dist)def_beff_cut=cutpair.first;
    }
  }
  else def_beff_cut=70;
  set_style()->cd();
}

bool bdt_base::process_sample(int bg,int sam)const{
  if(sam==0)return true;//always run on signal
  if(sam<0)return false;
  if(bg==1)return sam==3||sam==4||sam==5||sam==6;
  else if(is_single_bg(bg))return bg==sam;
  return true;
}

int bdt_base::which_sam(TString sam)const{
  bool zh=sam.Contains("125.root")&&(sam.Contains("ZH")||sam.Contains("WH")),tt=sam.Contains("ttbar.root"),db=(sam.Contains("ZZ")||sam.Contains("WW"))&&sam.Contains(".root");
  bool zlep=(sam.Contains("Zee")||sam.Contains("Zmumu")||sam.Contains("Ztautau"));bool zl=zlep&&sam.Contains("L.root"),zc=zlep&&sam.Contains("C.root"),zb=zlep&&sam.Contains("B.root");
  if(zh)return 0;
//   if(zl||zc||zb)return 1;
  if(tt)return 2;
  if(db)return 3;
  if(zl)return 4;
  if(zc)return 5;
  if(zb)return 6;
  return -99;
}
TString bdt_base::sam_tag(int sam)const{
  if(sam==0)return "-zh125";
  if(sam==1)return "-zjets";
  if(sam==2)return "-ttbar";
  if(sam==3)return "_zz";
  if(sam==4)return "-zl";
  if(sam==5)return "-zc";
  if(sam==6)return "-zb";
  return "-all";
}
TString bdt_base::sam_label(int sam)const{
  if(sam==0)return "ZH125";
  if(sam==1)return "Z+jets";
  if(sam==2)return "t#bar{t}";
  if(sam==3)return "ZZ";
  if(sam==4)return "Zl";
  if(sam==5)return "Zc";
  if(sam==6)return "Zb";
  return "Z+jets, t#bar{t}";
}
pair<double,double> bdt_base::cumulat_sig(TH1D*sig,TH1D*big)const{
  double ssb=0.;int nbins=sig->GetNbinsX();
  if(nbins!=big->GetNbinsX())return pair<double,double>(-99.,0.);
  for(int i=1;i<=nbins;i++){
    double S=sig->GetBinContent(i),B=big->GetBinContent(i);
    if(B<=0)B=1.e-6;
    ssb+=S/sqrt(S+B);
  }
  return {-999999.,ssb};
}
pair<double,double> bdt_base::optimal_sig(TH1D*sig,TH1D*big,bool simple)const{
  vector<double>sigs;int nbins=sig->GetNbinsX();
  if(nbins!=big->GetNbinsX())return pair<double,double>(-99.,0.);
  for(int i=1;i<=nbins;i++){
    double S=sig->Integral(i,-1),B=big->Integral(i,-1),ssb=S/sqrt(S+B);
    sigs.push_back(isfinite(ssb)?ssb:-99.);
  }
  int index=distance(sigs.begin(),max_element(sigs.begin(),sigs.end()));
  double significance=*max_element(sigs.begin(),sigs.end()),cut=sig->GetBinLowEdge(index+1);
  if(simple) return pair<double,double>(cut,significance);
  double peak=sig->GetBinCenter(sig->GetMaximumBin()),spread=sig->GetRMS();
  TString n=Form("%s_fit",sig->GetName());TF1*kyrie=new TF1(n,"gaus",peak-spread,peak+spread);TFitResultPtr eleison=sig->Fit(n,"RQS");
  if((Int_t)eleison!=0)return {cut,significance};//TMinuit status; 0 is normal--if peak finding fails, return the simple result
  cut=kyrie->GetParameter(1);
  delete kyrie;
  return cut_sig(sig,big,cut);
}
pair<double,double> bdt_base::cut_sig(TH1D*sig,TH1D*big,double cut)const{
  int nbins=sig->GetNbinsX(),bin=sig->GetXaxis()->FindBin(cut);
  if(nbins!=big->GetNbinsX())return {-99.,-99.};
  double S=sig->Integral(bin,-1),B=big->Integral(bin,-1),ssb=S/sqrt(S+B),finite_cut=cut;
  while(bin>0&&!isfinite(ssb)){
    bin--;
    S=sig->Integral(bin,-1);B=big->Integral(bin,-1);ssb=S/sqrt(S+B);finite_cut=sig->GetBinLowEdge(bin);
  }
  return {finite_cut,(isfinite(ssb)?ssb:0)};
}
pair<TH1D*,TH1D*> bdt_base::transformation_DF(TH1D*sig,TH1D*big,bool husk,bool d,double zs,double zb)const{
  if(d&&zs==4.5&&zb==4.5){
    double flep=0.9;//default is nothing, except we always assume two leptons
    double ftag=(beff_cut<=70?1.4:0.6);//if LL 1.4, MM/TT 0.6
    double fjet=(njet==2?0.8:1.2);//0.8 for 2 jet, 1.2 for 3 jet
    double fptv=(loptv?0.6:1.4),prefactor=flep*ftag*fjet*fptv;
    zs=prefactor*6.;zb=prefactor*4.;
    zs=5.;zb=5.;//hack for harmonization with Liv/Bir trans_D
  }
  TH1D *flip=(TH1D*)sig->Clone(),*flop=(TH1D*)big->Clone();//this is probably unecessary, but something funny is going on, so clone just in case some weird is going on with histograms in the file...
  vector<int>dragon=HistoTransform().getRebinBins(flop,flip,(d?6:12));  TString tag=(d?"_td":"_tf"),sn=Form("%s_%s%s",identifier.c_str(),flip->GetName(),tag.Data()),bn=Form("%s_%s%s",identifier.c_str(),flop->GetName(),tag.Data());
  flip->SetNameTitle(sn,sn);HistoTransform().rebinHisto(flip,&dragon);  flop->SetNameTitle(bn,bn);HistoTransform().rebinHisto(flop,&dragon);
  if(husk){flip->Reset();flop->Reset();}
  return {flip,flop};
  if(debug)cerr<<"Transformation D for signal("<<flip->GetName()<<"), background("<<flop->GetName()<<")..."<<endl;
  int hibin=flip->GetNbinsX();double min=flip->GetXaxis()->GetXmin(),max=flip->GetXaxis()->GetXmax();
  if(hibin!=flop->GetNbinsX())return pair<TH1D*,TH1D*>(NULL,NULL);
  if(min!=flop->GetXaxis()->GetXmin())return pair<TH1D*,TH1D*>(NULL,NULL);
  if(max!=flop->GetXaxis()->GetXmax())return pair<TH1D*,TH1D*>(NULL,NULL);
  vector<double>newbin_edges={min,max};vector<int>bin_partition={0,hibin};double llrs=0,NB=flop->Integral(),NS=flip->Integral(),nb=0,ns=0;
  while(hibin>0){
    if(hibin==1)break;
    double s=flip->GetBinContent(hibin),b=flop->GetBinContent(hibin);
    if(b<=0.){
      NB-=b;
      double avg=0.5*(flop->GetBinContent(hibin-1)+flop->GetBinContent(hibin+1));
      if(avg<=0.)b=1.e-6;
      else b=avg;
      NB+=b;
    }
    if(s<=0.){
      NS-=s;
      double avg=0.5*(flop->GetBinContent(hibin-1)+flop->GetBinContent(hibin+1));
      if(avg<=0.)s=1.e-6;
      else s=avg;
      NS+=s;
    }
    llrs+=s*log(1+s/b);nb+=b;ns+=s;double Z=(d?zs*ns/NS+zb*nb/NB:sqrt(zb*nb/NB)+sqrt(zs*llrs));
    if(debug)cerr<<"At bin "<<hibin<<", s="<<s<<","<<"b="<<b<<", nb up to "<<nb<<", llrs up to "<<llrs<<", for Z of "<<Z<<endl;
    if(Z>1.){
      llrs=0.;nb=0.;ns=0.;
      bin_partition.insert(bin_partition.begin()+1,hibin-1);
      newbin_edges.insert(newbin_edges.begin()+1,flip->GetBinLowEdge(hibin));
    }
    hibin--;
  }
  if(debug){
    cerr<<"Bin partitioning: "<<endl;
    for(auto d:bin_partition)cerr<<setw(8)<<d;
    cerr<<endl;
    for(auto d:newbin_edges)cerr<<setprecision(3)<<setw(8)<<d;
    cerr<<endl;
  }
  int nbins=bin_partition.size()-1;
  TH1D *patrick=new TH1D(sn,sn,nbins,&newbin_edges[0]),*jane=new TH1D(bn,bn,nbins,&newbin_edges[0]);
  if(husk)return {patrick,jane};
  for(int i=0;i<nbins;i++){
    if(debug)cerr<<"Fill (pos,s,b)=("<<0.5*(newbin_edges[i]+newbin_edges[i+1])<<","<<flip->Integral(bin_partition[i]+1,bin_partition[i+1])<<","<<flop->Integral(bin_partition[i]+1,bin_partition[i+1])<<")"<<endl;
    patrick->Fill(0.5*(newbin_edges[i]+newbin_edges[i+1]),flip->Integral(bin_partition[i]+1,bin_partition[i+1]));
    jane->Fill(0.5*(newbin_edges[i]+newbin_edges[i+1]),flop->Integral(bin_partition[i]+1,bin_partition[i+1]));
  }
  delete flip;delete flop;
  return pair<TH1D*,TH1D*>(patrick,jane);
}

void bdt_base::print_info_plot_2d(const vector<vector<double> >&info,const vector<TString>&xlabels,const vector<TString>&ylabels,const TString&name,const TString&pdir,const TString&title)const{
  if(info.size()!=xlabels.size())return;
  if(ylabels.empty())return;
  if(info.front().size()!=ylabels.size())return;
  TStyle*style=set_style();  style->SetPadRightMargin(0.075);style->cd();TString doy=check_mkdir(pdir)+name+".pdf";
  TCanvas*c=new TCanvas(doy,doy,100*xlabels.size(),100*ylabels.size());c->cd();
  TH2D*compendium=new TH2D(name,name,xlabels.size(),-0.5,xlabels.size()-0.5,ylabels.size(),-0.5,ylabels.size()-0.5);
  compendium->GetXaxis()->SetLabelSize(0.035);  compendium->GetYaxis()->SetLabelSize(0.035);//  compendium->GetZaxis()->SetLabelSize(0.035);
  for(int x=0;x<(int)xlabels.size();x++)compendium->GetXaxis()->SetBinLabel(x+1,xlabels[x]);
  for(int y=0;y<(int)ylabels.size();y++)compendium->GetYaxis()->SetBinLabel(y+1,ylabels[y]);
  for(int x=0;x<(int)xlabels.size();x++){
    for(int y=0;y<(int)ylabels.size();y++)compendium->Fill(x,y,info[x][y]);
  }
  compendium->Draw("coltextz");
  TLatex l=plot_latex();l.SetTextSize(0.03);
  l.DrawLatex(0.08,0.95,title);
  c->Print(doy);
  delete compendium;delete c;
}

TString bdt_base::btag_var(int ratio,bool cts)const{
  TString var="MV2c";
  if(ratio==0)var+="00";
  else if(ratio==10)var+="10";
  else var+="20";
  if(!cts)var+="bin";
  return var;
}

TString bdt_base::cut_str(int cut,int ratio)const{
  if(btag_ratio_cuts.find(ratio)==btag_ratio_cuts.end())ratio=def_cratio;
  if(btag_ratio_cuts.at(ratio).find(cut)==btag_ratio_cuts.at(ratio).end()) cut=70;
  TString beff="-beff";
  beff+=cut;
  return beff;
}

TString bdt_base::cut_label(int cut,int ratio)const{
  if(btag_ratio_cuts.find(ratio)==btag_ratio_cuts.end())ratio=def_cratio;
  if(btag_ratio_cuts.at(ratio).find(cut)==btag_ratio_cuts.at(ratio).end()) cut=70;
  TString beff;beff+=cut;beff+="% b-tag eff";
  return beff;
}

TString bdt_base::region_str(int nt,int nj,bool lptv)const{
  TString str="-";
  //tags
  if(nt<0)str+=0;
  else if(nt>2)str+=2;
  else str+=nt;
  str+="tag";

  //jets
  if(nj<=2)str+=2;
  else if(nj==-99)str+="2p";
  else str+=3;//3p
  str+="jet-";

  if(lptv==1)str+="0_150ptv";
  else if(lptv==0) str+="150ptv";
  else str+="allptv";

  return str;
}

TString bdt_base::region_label(int nt,int nj,bool lptv)const{
  TString str;
  //tags
  if(nt<0)str+=0;
  else if(nt>2)str+=2;
  else str+=nt;
  str+="tag";

  //jets
  if(nj<=2)str+=2;
  else if(nj==-99)str+="2p";
  else str+=3;
  str+="jet, ";

  if(lptv==1)str+="p_{T}^{V}<150 GeV";
  else if(lptv==0) str+="p_{T}^{V}>150 GeV";

  return str;
}

TString bdt_base::tdf_tstr(int tDF)const{
  if(tDF==1)return "_tF";
  if(tDF==2)return "_tD";
  return "";
}

//short for the FixedCutBEff_x weight cuts
double bdt_base::mv2c_cut(int cut,int ratio)const{
  if(btag_ratio_cuts.find(ratio)==btag_ratio_cuts.end())ratio=def_cratio;//if the charm ratio isn't supported, revert to default working point cuts
  if(btag_ratio_cuts.at(ratio).find(cut)!=btag_ratio_cuts.at(ratio).end()) return btag_ratio_cuts.at(ratio).at(cut);

  return -0.0436;//MV2c20 70% efficiency working point default
}

TString bdt_base::cuts(int nt,int nj,int cut,int ratio,bool cts,bool lptv,bool wmll,bool met)const{
  double bcut=mv2c_cut(cut,ratio);TString var=btag_var(ratio,cts),b1=var+"B1",b2=var+"B2";
  ostringstream cuts;
  //common dR cut
  cuts<<"(pTV/1000.>=200 || (pTV/1000.<200&&dRBB>0.7))";
  //pTV cut
  if(lptv==1)cuts<<"&& pTV / 1000. < 150";
  else if(lptv==0) cuts<<"&& pTV / 1000. >= 150";
  //default is inclusive

  //mLL cut
  if(wmll)cuts<<" && mLL / 1000. > 71 && mLL / 1000. < 121.";
  else cuts<<" && mLL / 1000. > 83 && mLL / 1000. < 99.";

  //njet cut
  if(nj>=3)cuts<<" && nSigJet>=3";
  else if(nj==-99)cuts<<" && nSigJet>=2";
  else cuts<<" && nSigJet==2";

  //btag cut
  if(use_btagvs||always_tag){
    if(nt<=0)cuts<<"&& "<<b1<<"<"<<bcut<<" && "<<b2<<"<"<<bcut;
    else if(nt==1)cuts<<"&& (("<<b1<<">="<<bcut<<" && "<<b2<<"<"<<bcut<<") || ("<<b1<<"<"<<bcut<<" && "<<b2<<">="<<bcut<<"))";
    else cuts<<"&& "<<b1<<">="<<bcut<<" && "<<b2<<">="<<bcut;
  }
  //MET cut (60 GeV cut)
  if(met)cuts<<"&& MET/1000.<60.";
  return TString(cuts.str());
}

string bdt_base::check_dir(const string&dir)const{
  string gooddir(dir);
  if(gooddir.substr(gooddir.length()-1)!="/")gooddir+="/";
  return gooddir;
}

string bdt_base::check_isdir(const string&dir)const{
  string gooddir=check_dir(dir);
  if(!file_exists(gooddir.substr(0,gooddir.length()-1))){
    cout<<"Directory "<<gooddir<<" does not exist; check failed.  Aborting."<<endl;
    exit(9);
  }
  return gooddir;
}

string bdt_base::check_mkdir(const string&dir)const{
  string gooddir(dir);
  if(gooddir.substr(gooddir.length()-1)!="/")gooddir+="/";
  if(!file_exists(gooddir.substr(0,gooddir.length()-1)))system(("mkdir -p "+gooddir).c_str());
  return gooddir;
}

bool bdt_base::file_exists(const string&name)const{
  if(debug)cout<<"checking to see if "<<name<<" exists."<<endl;
  string path="./",file=name;
  if(file.find("/")!=string::npos){
    path=name.substr(0,name.find_last_of("/")+1);file=name.substr(name.find_last_of("/")+1);
  }
  DIR*red=opendir(path.c_str());struct dirent*andy;
  if(red==NULL)return false;
  while((andy=readdir(red))!=NULL){
    if(string(andy->d_name)==file)return true;
  }
  return false;
}

int bdt_base::color(int which)const{
  vector<int> colors={(int)kPink+7,(int)kRed,(int)kOrange+7,(int)kOrange-2,(int)kSpring-5,(int)kGreen+2,(int)kTeal-5,(int)kCyan+2,(int)kAzure+1,(int)kBlue+2,(int)kViolet-2,(int)kMagenta+2,(int)kPink+6,(int)kRed+2,(int)kOrange-3,(int)kYellow-3,(int)kSpring-2,(int)kGreen+1,(int)kTeal+2,(int)kCyan-2,(int)kAzure-4,(int)kBlue-3,(int)kViolet+1,(int)kMagenta-7,(int)kWhite,(int)kGray,(int)kGray+1,(int)kGray+2,(int)kGray+3,(int)kBlack};
  return colors[(which+1)%colors.size()];
}

TStyle* bdt_base::set_style() const{
  TStyle *style = AtlasStyle();//new TStyle("Plain","");
  style->SetNumberContours(100); 
  style->SetPadRightMargin(0.075);
  style->SetPadLeftMargin(0.10);
  style->SetPadTopMargin(0.075);
  style->SetPadBottomMargin(0.11);
  int n_pts=11;
//makes standard color palette slightly pastel (purple to red) with: crimson at red, lighter purple and 
  Double_t R[11]={1.00,0.60,0.20,0.20,0.40,0.20,0.20,0.60,1.00,1.00,0.80};
  Double_t G[11]={0.60,0.20,0.20,0.60,1.00,1.00,1.00,1.00,1.00,0.60,0.00};
  Double_t B[11]={1.00,1.00,1.00,1.00,1.00,0.60,0.20,0.20,0.20,0.20,0.00},length[11];
  for(int i=0;i<n_pts;i++)length[i]=0.1*i;
  TColor::CreateGradientColorTable(n_pts,length,R,G,B,100);
//   style->SetPalette(1,0);
  style->cd();
  return style;
}

TLatex bdt_base::plot_latex() const{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextSize(0.06);
//   l.SetTextFont(132);
  //label text color: black
  l.SetTextColor(kBlack);
  return l;
}

void bdt_base::auto_limits(vector<TH1D*> loki,double zero)const{
//   double fat_man=buff_hi*(*max_element(smoke.begin(),smoke.end())),lil_boy=buff_lo*(*min_element(mirrors.begin(),mirrors.end()));
//   if(lil_boy<0)lil_boy/=(buff_lo*buff_lo);
  vector<double> smoke,mirrors,shade;
  for(unsigned int jhist=0; jhist<loki.size();jhist++){
    for(int ibin=1;ibin<=loki[jhist]->GetNbinsX();ibin++){
      double val=loki[jhist]->GetBinContent(ibin),err=loki[jhist]->GetBinError(ibin);
      smoke.push_back(val+err);mirrors.push_back(val-err);shade.push_back(val);
    }
  }
  double fat_man=(*max_element(smoke.begin(),smoke.end())),lil_boy=(*min_element(mirrors.begin(),mirrors.end())),interval=fat_man-lil_boy,hval=fat_man;//(*max_element(shade.begin(),shade.end()));
  double bumper=(0.1*interval>(hval-zero)*0.01?0.12*interval:(hval-zero)*0.01);
//   cout<<"in hist "<<loki->GetName()<<", interval/10="<<0.1*interval<<", hval/100="<<hval/100.<<", bumper="<<bumper<<endl;
  bumper=0.15*interval;
  fat_man+=4.*bumper;lil_boy-=0.75*bumper;
  for(unsigned int i=0;i<loki.size();i++)loki[i]->SetAxisRange(lil_boy,fat_man,"Y");
}

TString bdt_base::ftag(const vector<TString>&vars)const{
  TString tag("-");tag+=identifier;tag+=btag_str();tag+=region_str();
  for(auto v:vars)tag+=("-"+v);
  return tag;
}
