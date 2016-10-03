#include "sig_plotter.h"

sig_plotter::sig_plotter(bool ks):keep_spectators(ks){}

void sig_plotter::do_plots(const string&sigdir,const string&pdir)const{
  set_style()->cd();
  TString ext=".pdf";
  TCanvas *c=new TCanvas(sigdir.c_str(),sigdir.c_str(),500,300);
  for(int iv=n_std_min();iv<=n_li_cases();iv++){
    std::vector<TString>vars=li_vars(iv).first;TString ltg=li_tag(iv);
    for(int bg=-1;bg<n_bg_cases();bg++){
      TString btg=bg_tag(bg);
      //plotting code is the first; variable type preceded by numbers to get order correct; each will be a separate entry in the map
      map<pair<TString,int>,pair<vector<double>,vector<double> > >grinfo;
      print_bdts(c,sigdir,pdir,iv,bg);
      for(int i=0;i<(int)vars.size();i++){
	pair<TString,int>key=var_type(vars[i]);
	if(grinfo.find(key)==grinfo.end()){
	  grinfo[key]=pair<vector<double>,vector<double> >(vector<double>(1,i),vector<double>(1,sig_from_files(sigdir,iv,bg,i).second));
	}
	else{
	  grinfo[key].first.push_back(i);grinfo[key].second.push_back(sig_from_files(sigdir,iv,bg,i).second);
	}
      }
      TString lol="sig_v_var_plot"+ltg+btg;
      TMultiGraph *telemachus=new TMultiGraph(lol,lol);
      for(auto&it:grinfo){
	TGraph*jabba=new TGraph(it.second.first.size(),&it.second.first[0],&it.second.second[0]);
	for(int i=0;i<(int)it.second.second.size();i++){
	  ostringstream num;num<<setprecision(3)<<it.second.second[i];
	  TLatex *latex = new TLatex(jabba->GetX()[i],jabba->GetY()[i],num.str().c_str());
	  latex->SetTextSize(0.04);
	  jabba->GetListOfFunctions()->Add(latex);
	}
	TString hutt(lol);hutt+=it.first.first;jabba->SetNameTitle(hutt,hutt);
	jabba->SetLineColor(it.first.second);jabba->SetMarkerColor(it.first.second);jabba->SetMarkerStyle(kFullTriangleDown);jabba->SetLineWidth(3);
	telemachus->Add(jabba);
      }
      telemachus->Draw("ACP");
      telemachus->GetYaxis()->SetTitle("S/#sqrt{S+B}");
      TAxis *xax=telemachus->GetXaxis();
      for(int i=0;i<(int)vars.size();i++) xax->SetBinLabel(xax->FindBin(i*1.),vars[i]);
      double totsig=grinfo.rbegin()->second.second.back(); TLatex l=plot_latex();
      ostringstream num;num<<"Total Sig = "<<setprecision(3)<<totsig<<", bg: "<<bg_label(bg);
      l.DrawLatex(0.2,0.8,num.str().c_str());
      c->Print(pdir+lol+ext);
    }
  }
}

TH1D* sig_plotter::bdt_dist(const string&file,int sbg)const{
  TFile*f=TFile::Open(file.c_str());
  if(f==NULL)return new TH1D();
  TDirectoryFile*df=(TDirectoryFile*)f->Get("Method_BDT/BDT");
  if(df==NULL)return new TH1D();
  TString hname="MVA_BDT_";
  hname+=(sbg?"S":"B");
  TH1D*hiro=(TH1D*)df->Get(hname);
  hiro->SetFillColorAlpha((sbg?kRed+2:kBlue+2),0.75);
  hiro->SetLineColor(kBlack);
  hiro->GetXaxis()->SetTitle("BDT score");
  hiro->GetYaxis()->SetTitle("a.u.");
//   df->Close();
//   f->Close();
//   delete df;  delete f;
  return hiro;
}

void sig_plotter::print_bdts(TCanvas*c,const string&ddir,const string&odir,int vt,int bg,int iv)const{
  ostringstream tag;tag<<li_tag(vt)<<bg_tag(bg);
  string file=ddir+"TMVA"+tag.str()+".root";
  TStyle *style=set_style();
  TH1D*s=bdt_dist(file,1),*b=bdt_dist(file,0);
  vector<TH1D*>yup;yup.push_back(s);yup.push_back(b);
  pair<double,double> sig=sig_from_files(ddir,vt,bg,iv);
  double lo=style->GetPadLeftMargin(),hi=1.-style->GetPadRightMargin();
  double xmin=s->GetXaxis()->GetXmin(),xmax=s->GetXaxis()->GetXmax(),cut_val=sig.first,line_x=lo+(hi-lo)*(cut_val-xmin)/(xmax-xmin);
//   cout<<"cut: "<<sig.first<<", sig: "<<sig.second<<", lo="<<lo<<", hi="<<hi<<", xmin="<<xmin<<", xmax="<<xmax<<", line_x="<<line_x<<", s "<<s<<", b "<<b<<endl;
  c->cd();
  s->Draw("hist");b->Draw("hist same");
  TLine*line=new TLine();
  line->SetLineColor(kSpring-2);
  line->SetLineWidth(3);
  line->DrawLineNDC(line_x,0.125,line_x,0.95);
  TLatex l=plot_latex();l.SetTextSize(0.04);
  ostringstream cut;cut<<"Cut="<<setprecision(3)<<cut_val<<", sig="<<sig.second;
  l.DrawLatex(line_x+0.01,0.8,cut.str().c_str());
  c->Print((odir+"bdt_dists"+tag.str()+".pdf").c_str());
  delete s;delete b;
}

pair<double,double> sig_plotter::sig_from_files(const string&ddir,int vt,int bg,int iv)const{
  //construct the TMVA output and log files
  ostringstream log,tmva;  log<<ddir<<"tmva_out"<<li_tag(vt)<<bg_tag(bg)<<ith_var_tag(vt,iv)<<".txt";  tmva<<ddir<<"TMVA"<<li_tag(vt)<<bg_tag(bg)<<ith_var_tag(vt,iv)<<".root";

  //read the log file to get the normalizations for signal and background
  ifstream fin(log.str());
  double sig_evts=0,bg_evts=0,bdtcut=0;string line;
  getline(fin,line);
  while(fin.good()){
    if(line.find("number of events passed:")!=string::npos){
      if(line.find("Signal")!=string::npos)sig_evts=atof(line.substr(line.find("weights:")+8).c_str());
      else if(line.find("Background")!=string::npos)bg_evts=atof(line.substr(line.find("weights:")+8).c_str());
    }
    else if(line.find(" BDT:")!=string::npos) bdtcut=atof(line.substr(line.find("BDT:")+4).c_str());
    getline(fin,line);
  }

  TFile*f=TFile::Open(tmva.str().c_str());
  if(f==NULL)return pair<int,double>(-3,0.);
  TDirectoryFile*df=(TDirectoryFile*)f->Get("Method_BDT/BDT");
  if(df==NULL)return pair<int,double>(-3,0.);
  TH1D*sig=(TH1D*)df->Get("MVA_BDT_S"),*big=(TH1D*)df->Get("MVA_BDT_B");
  pair<double,double> fun=optimal_sig(sig,big);
  double snorm=sig_evts/sig->Integral(),bnorm=bg_evts/big->Integral();
  int lobin=sig->GetXaxis()->FindBin(bdtcut);  int nbins=sig->GetNbinsX();
  vector<double>sigs;
  double St=sig->Integral(lobin,-1)*snorm,Bt=big->Integral(lobin,-1)*bnorm,filesig=St/sqrt(St+Bt);
  delete sig;delete big;f->Close();
  if(filesig>fun.second)return pair<double,double>(bdtcut,filesig);
//   df->Close();//   f->Close();//   delete df;//   delete f; 
  return fun;
}

pair<double,double> sig_plotter::optimal_sig(TH1D*sig,TH1D*big)const{
  vector<double>sigs;
  for(int i=1;i<=nbins;i++){
    double S=sig->Integral(i,-1)*snorm,B=big->Integral(i,-1)*bnorm;
    sigs.push_back(S/sqrt(S+B));
  }
  int index=distance(sigs.begin(),max_element(sigs.begin(),sigs.end()));double significance=sigs[0],cut=bdtcut;
  significance=*max_element(sigs.begin(),sigs.end());
  double xmin=sig->GetXaxis()->GetXmin(),xmax=sig->GetXaxis()->GetXmax();cut=xmin+(xmax-xmin)*(index+0.5)/nbins;
  return pair<double,double>(cut,significance);
}
double sig_plotter::sig_from_log(const string&fname)const{
  ifstream fin(fname);string word;double sig;
  fin>>word;
  while(fin.good()){
    if(word=="Significance")fin>>sig;
    fin>>word;
  }
  return sig;
}

bool sig_plotter::is_bkg(TString bg,int which)const{
  bool tt=bg.Contains("ttbar.root"),diboson=bg.Contains("ZZ.root"),zj=(bg.Contains("Zee")||bg.Contains("Zmumu")||bg.Contains("Ztautau"))&&(bg.Contains("L.root")||bg.Contains("C.root")||bg.Contains("B.root"));
  if(which==0) return zj;
  else if(which==1)return tt; 
  else if(which==2)return diboson;
  return zj||tt||diboson;
}

TString sig_plotter::ith_var_tag(int vt,int iv)const{
  if(iv>=0&&iv<(int)li_vars(vt).first.size()-1){
    TString tag="_";tag+=(iv+1);
    return tag;
  }
  return "";
}

TString sig_plotter::li_tag(int which)const{
  if(which==0)return "_li_met";
  else if(which==1) return "_li_ll";
  else if(which==2) return "_li_bb";
  else if(which==3) return "_li_bb_ll";
  else if(which==4) return "_li_only";
  else if(which==5) return "_li_lb0";
  else if(which==-1)return "_std";
  else return "_cut";
}
TString sig_plotter::bg_tag(int bg)const{
  if(bg==0)return "_zjets";
  else if(bg==1)return "_ttbar";
  else if(bg==2)return "_zz";
  return "_all";
}
TString sig_plotter::bg_label(int bg)const{
  if(bg==0)return "Z+jets";
  else if(bg==1)return "t#bar{t}";
  else if(bg==2)return "ZZ";
  return "Z+jets&t#bar{t}";
}
//figure out how to sort variables in the right order for orderly and consistent training
//order should look like:
// LI inner products(ips)/mBB,dRBB, LI angles/other std variables, btag, LI inv ips/MET, LI inv angles

std::pair<TString,int> sig_plotter::var_type(const TString& var)const{
  int iWhite =0,   iBlack =1,   iGray=920,  iRed   =632, iGreen =416, iBlue=600, iYellow=400, iMagenta=616, iCyan=432,  iOrange=800, iSpring=820, iTeal=840, iAzure =860, iViolet =880, iPink=900;
  int colors[]={iPink+7,iRed,iOrange+7,iOrange-2,iSpring-5,iGreen+2,iTeal-5,iCyan+2,iAzure+1,iBlue+2,iViolet-2,iMagenta+2,iPink+6,iRed+2,iOrange-3,iYellow-3,iSpring-2,iGreen+1,iTeal+2,iCyan-2,iAzure-4,iBlue-3,iViolet+1,iMagenta-7,iWhite,iGray,iGray+1,iGray+2,iGray+3,iBlack};
  int one_cycle=12,skip=1,skip2=one_cycle+skip,step=2;
  if(var=="mBB"||var=="dRBB") return std::pair<TString,int>("00_cut",colors[skip]);
  else if((var.Contains("0_")||var.Contains("1_"))&&(var.Contains("_j")||var.Contains("_l"))&&!var.Contains("inv")&&!var.Contains("angle"))return std::pair<TString,int>("01_li_ip",colors[skip2]);

//   else if(var=="dEtaBB"||var=="pTV"||var=="pTB1"||var=="pTB2"||var=="dPhiVBB"||var=="dEtaVBB"||var=="mLL")return std::pair<TString,int>("02_std",iOrange-2);
  else if(var.Contains("angle")&&!var.Contains("inv"))return std::pair<TString,int>("03_li_ang",colors[1*step+skip2]);

  else if(var.Contains("MV2"))return std::pair<TString,int>("04_btag",colors[2*step+skip]);

  else if(var=="MET")return std::pair<TString,int>("05_met",colors[3*step+skip]);
  else if(!var.Contains("angle")&&var.Contains("inv"))return std::pair<TString,int>("06_li_inv_ip",colors[3*step+skip2]);

  else if(var.Contains("angle")&&var.Contains("inv"))return std::pair<TString,int>("07_li_inv_ang",colors[4*step+skip2]);

  return std::pair<TString,int>("02_std",colors[1*step+skip]);
}

std::vector<TString> sig_plotter::sorted_vars(const std::vector<TString>&vars)const{
  std::map<TString,std::vector<TString> >sorter;
  for(auto var:vars){
    std::pair<TString,int>spater(var_type(var));
    if(sorter.find(spater.first)==sorter.end())sorter[spater.first]=std::vector<TString>(1,var);
    else sorter[spater.first].push_back(var);
  }
  std::vector<TString> shiny;
  for(auto& it:sorter){
    for(auto var:it.second)shiny.push_back(var);
  }
  return shiny;
}
std::pair<std::vector<TString>,std::vector<TString> > sig_plotter::li_vars(int which)const{
  //training variables
  std::vector<TString> train_std={"mBB","dRBB","dEtaBB","pTV","pTB1","pTB2","dPhiVBB","dEtaVBB","mLL"};
//   std::vector<TString> train_li={"j0_j1", "j0_l0", "j0_l1", "j1_l0", "j1_l1", "l0_l1", "angle_bb_z", "angle_bbz_bbll"};
  std::vector<TString> train_li={"j0_j1", "j0_l0", "j0_l1", "j1_l0", "j1_l1", "l0_l1", "angle_bb_z", "angle_bbz_bbll"};
  std::vector<TString> train_li_met_full={"inv_bb", "inv_l0", "angle_inv_bb", "angle_inv_l0","inv_lb0","angle_inv_lb0"};
  std::vector<TString> train_li_met_ll={"inv_ll", "angle_inv_ll"},train_li_met_bb={"inv_bb", "angle_inv_bb"},train_li_met_lb0={"inv_lb0","angle_inv_lb0"};
  std::vector<TString> train_met={"MET"},train_common={"MV2c20B1","MV2c20B2"};
  std::vector<TString> vars={"mBB","dRBB"},
    spectators={"dEtaBB","pTV","pTB1","pTB2","dPhiVBB","dEtaVBB","dEtaBB","mLL","MV2c20B1","MV2c20B2",
		"MET","j0_j1", "j0_l0", "j0_l1", "j1_l0", "j1_l1", "l0_l1", "angle_bb_z", "angle_bbz_bbll",
		"inv_bb", "inv_ll", "angle_inv_bb", "angle_inv_ll"};
  if(which>=0&&which<n_li_cases()){
    vars=train_li;spectators=train_std;
    //MET
    if(which==0){
      vars.insert(vars.end(),train_met.begin(),train_met.end());
      spectators.insert(spectators.end(),train_li_met_full.begin(),train_li_met_full.end());
    }
    else{
      spectators.insert(spectators.end(),train_met.begin(),train_met.end());
      //ll inv
      if(which==1){
	vars.insert(vars.end(),train_li_met_ll.begin(),train_li_met_ll.end());
	spectators.insert(spectators.end(),train_li_met_bb.begin(),train_li_met_bb.end());
	spectators.insert(spectators.end(),train_li_met_lb0.begin(),train_li_met_lb0.end());
      }
      //bb inv
      else if(which==2){
	vars.insert(vars.end(),train_li_met_bb.begin(),train_li_met_bb.end());
	spectators.insert(spectators.end(),train_li_met_ll.begin(),train_li_met_ll.end());
	spectators.insert(spectators.end(),train_li_met_lb0.begin(),train_li_met_lb0.end());
      }
      //all inv
      else if(which==3){
	vars.insert(vars.end(),train_li_met_full.begin(),train_li_met_full.end());
      }
      //LI only, no invisibles
      else if(which==4){
	spectators.insert(spectators.end(),train_li_met_full.begin(),train_li_met_full.end());
      }
      //lb0 inv
      else if(which==5){
	vars.insert(vars.end(),train_li_met_lb0.begin(),train_li_met_lb0.end());
	spectators.insert(spectators.end(),train_li_met_ll.begin(),train_li_met_ll.end());
	spectators.insert(spectators.end(),train_li_met_bb.begin(),train_li_met_bb.end());
      }
    }
    vars.insert(vars.end(),train_common.begin(),train_common.end());
  }
  else if(which>=n_std_min()&&which<0){
    vars=train_std;spectators=train_li;
    if(which==-1)vars.insert(vars.end(),train_met.begin(),train_met.end());
    spectators.insert(spectators.end(),train_li_met_full.begin(),train_li_met_full.end());
    vars.insert(vars.end(),train_common.begin(),train_common.end());
  }
  std::pair<std::vector<TString>,std::vector<TString> >derp(sorted_vars(vars),keep_spectators?sorted_vars(spectators):std::vector<TString>());
  return derp;
}

void sig_plotter::auto_limits(vector<TH1D*> loki,double zero)const{
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

TStyle* sig_plotter::set_style() const{
  TStyle *style = AtlasStyle();//new TStyle("Plain","");
  style->SetNumberContours(100);
  style->SetPadRightMargin(0.10);
  style->SetPadLeftMargin(0.16);
  style->SetPadTopMargin(0.075);
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

TLatex sig_plotter::plot_latex() const{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextSize(0.06);
//   l.SetTextFont(132);
  //label text color: black
  l.SetTextColor(kBlack);
  return l;
}

