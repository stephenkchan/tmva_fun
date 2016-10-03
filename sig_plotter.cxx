#include "sig_plotter.h"

sig_plotter::sig_plotter(const string&bdir, int _tr_opt,int _bg,bool ks):bdt_dir(bdir),tr_opt(_tr_opt),bg(_bg),keep_spectators(ks){
  vars=li_vars(tr_opt).first;nv=vars.size();
}

double sig_plotter::sig_in_file(const string&fname)const{
  double sig=0;string word;
  ifstream fin(fname);
  while(fin.good()){
    fin>>word;
    if(word=="Significance"){fin>>sig;break;}
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

TString sig_plotter::li_tag(int which)const{
  if(which==0)return "_li_met";
  else if(which==1) return "_li_ll";
  else if(which==2) return "_li_bb";
  else if(which==3) return "_li_bb_ll";
  else if(which==4) return "_li_only";
  else if(which==-1)return "_std";
  else if(which==-2)return "_std_met";
  else return "_cut";
}
TString sig_plotter::bg_tag(int bg)const{
  if(bg==0)return "_zjets";
  else if(bg==1)return "_ttbar";
  else if(bg==2)return "_zz";
  return "_all";
}
//figure out how to sort variables in the right order for orderly and consistent training
//order should look like:
// LI inner products(ips)/mBB,dRBB, LI angles/other std variables, btag, LI inv ips/MET, LI inv angles

std::pair<TString,int> sig_plotter::var_type(const TString& var)const{
  int iWhite =0,   iBlack =1,   iGray=920,  iRed   =632, iGreen =416, iBlue=600, iYellow=400, iMagenta=616, iCyan=432,  iOrange=800, iSpring=820, iTeal=840, iAzure =860, iViolet =880, iPink=900;
  int colors[]={iPink+7,iRed,iOrange+7,iOrange-2,iSpring-5,iGreen+2,iTeal-5,iCyan+2,iAzure+1,iBlue+2,iViolet-2,iMagenta+2,iPink+6,iRed+2,iOrange-3,iYellow-3,iSpring-2,iGreen+1,iTeal+2,iCyan-2,iAzure-4,iBlue-3,iViolet+1,iMagenta-7};
  int one_cycle=12,skip=1,skip2=one_cycle+skip,step=2;
  if(var=="mBB"||var=="dRBB") return std::pair<TString,int>("00_cut",colors[skip]);
  else if((var.Contains("0_")||var.Contains("1_"))&&(var.Contains("_j")||var.Contains("_l"))&&!var.Contains("inv")&&!var.Contains("angle"))return std::pair<TString,int>("01_li_ip",colors[skip2]);

//   else if(var=="dEtaBB"||var=="pTV"||var=="pTB1"||var=="pTB2"||var=="dPhiVBB"||var=="dEtaVBB"||var=="dEtaBB"||var=="mLL")return std::pair<TString,int>("02_std",iOrange-2);
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
  std::vector<TString> train_std={"mBB","dRBB","dEtaBB","pTV","pTB1","pTB2","dPhiVBB","dEtaVBB","dEtaBB","mLL"};
  std::vector<TString> train_li={"j0_j1", "j0_l0", "j0_l1", "j1_l0", "j1_l1", "l0_l1", "angle_bb_z", "angle_bbz_bbll"};
  std::vector<TString> train_li_met_full={"inv_bb", "inv_ll", "angle_inv_bb", "angle_inv_ll"};
  std::vector<TString> train_li_met_ll={"inv_ll", "angle_inv_ll"},train_li_met_bb={"inv_bb", "angle_inv_bb"},train_met={"MET"},train_common={"MV2c20B1","MV2c20B2"};
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
      }
      //bb inv
      else if(which==2){
	vars.insert(vars.end(),train_li_met_bb.begin(),train_li_met_bb.end());
	spectators.insert(spectators.end(),train_li_met_ll.begin(),train_li_met_ll.end());
      }
      //all inv
      else if(which==3){
	vars.insert(vars.end(),train_li_met_full.begin(),train_li_met_full.end());
      }
      //LI only, no invisibles
      else spectators.insert(spectators.end(),train_li_met_full.begin(),train_li_met_full.end());
    }
  }
  else if(which>=n_std_min()){
    vars=train_std;spectators=train_li;
    if(which==-1)vars.insert(vars.end(),train_met.begin(),train_met.end());
    else if(which==-2)spectators.insert(spectators.end(),train_met.begin(),train_met.end());
    spectators.insert(spectators.end(),train_li_met_full.begin(),train_li_met_full.end());
  }
  vars.insert(vars.end(),train_common.begin(),train_common.end());
  std::pair<std::vector<TString>,std::vector<TString> >derp(sorted_vars(vars),keep_spectators?sorted_vars(spectators):std::vector<TString>());
  return derp;
}

