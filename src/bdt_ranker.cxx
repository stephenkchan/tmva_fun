#include "bdt_ranker.h"

bdt_ranker::bdt_ranker(const string&in,const string&out,const vector<TString>&bs,const vector<TString>&vs,int nt,int nj,int eff_cut,int ratio,bool cts,int lptv,bool btv,bool wmll,bool met):
  bdt_base(in,out,"rank",nt,nj,eff_cut,ratio,cts,lptv,btv,wmll,met),ranked(bs),unranked(vs){}

bdt_ranker::bdt_ranker(const string&in,const string&out,const vector<TString>&bs,const vector<TString>&vs,const bdt_base&base):
  bdt_base(base),ranked(bs),unranked(vs){}

vector<double> bdt_ranker::ranking_front_matter(const string&tag,bool overwrite,int trans_DF){
  string pdir=check_mkdir(tmva_out_dir+"bdt_spikes");
  set_tag(tag+tdf_str(trans_DF));
  cout<<"Begin ranking: "<<log()<<" "<<tmva_out_dir<<" "<<these_cuts()<<endl;
  vector<TString> done(ranked),to_do(unranked);  bdt_validate starter(done,overwrite,(bdt_base)*this);
  vector<double>daisy_chain(1,starter.val_sig(-1,trans_DF).second);
  ostringstream samples,events; pair<double,double>master=starter.evts_S_B(-1);samples<<setw(20)<<"ZHll125"<<setw(20)<<" TOTAL_BG";events<<setw(20)<<master.first<<setw(20)<<master.second;
  for(auto i:{2,3,4,5,6}){samples<<setw(20)<<sam_label(i);events<<setw(20)<<starter.evts_S_B(i).second;}
  ofstream fout(log().c_str());
  fout<<"START BDT ranking of type: "<<identifier<<setprecision(5)<<endl;
  fout<<">-----------------------------------EVENT COUNTS-----------------------------------------------<"<<endl;
  fout<<samples.str()<<endl<<events.str()<<endl;
  fout<<">----------------------------------------------------------------------------------------------<"<<endl;
  fout<<"RANKED VARIABLES: ";
  for(auto a:ranked)fout<<a<<" ";
  fout<<endl<<"UNRANKED VARIABLES: ";
  for(auto a:unranked)fout<<a<<" ";
  fout<<endl<<"Significance starts at "<<daisy_chain.back()<<endl;
  fout<<"------------------------------------------------------------------------------------------------"<<endl;
  fout.close();
  return daisy_chain;
}
void bdt_ranker::do_ranking(const string&tag,bool overwrite,int trans_DF){
  vector<double>daisy_chain=ranking_front_matter(tag,overwrite,trans_DF);
  vector<TString> done(ranked),to_do(unranked);
  ofstream fout(log().c_str(),std::ofstream::app);
  while(!to_do.empty()){
    rank_loop_one_var(to_do,done,daisy_chain,fout,tag,overwrite,trans_DF);
  }
  ranking_back_matter(done,daisy_chain,fout,overwrite,trans_DF);
}

void bdt_ranker::finish_ranking(const string&tag,bool overwrite,int trans_DF){
  set_tag(tag+tdf_str(trans_DF));  string lname=log();
  if(!file_exists(lname)||overwrite){do_ranking(tag,true,trans_DF);return;}
  ifstream hey(lname.c_str());vector<TString>to_do(unranked),done(ranked);string line,prevline;vector<double>daisy_chain;
  map<TString,double>buffer;bool buffer_was_last=false;
  while(hey.good()){
    prevline=line;
    getline(hey,line);
    if(line.find("UNRANKED VARIABLES: ")==0){
      buffer_was_last=false;
      vector<TString>daredevil(unranked);vector<string>nlist=line_to_words(line.substr(line.find(":")+1));
      for(auto word:nlist){
	int killme=-1;
	for(int i=0;i<(int)daredevil.size();i++){
	  if(word==string(daredevil[i].Data())){killme=i;break;}
	}
	if(killme>=0)daredevil.erase(daredevil.begin()+killme,daredevil.begin()+killme+1);
      }
      if(!daredevil.empty()){
	do_ranking(tag,true,trans_DF);
	return;
      }
    }
    else if(line.find("Significance starts at ")==0){
      buffer_was_last=false;
      daisy_chain.push_back(atof(line.substr(line.find("at")+3).c_str()));
    }
    else if(line.find("BUFFER: ")==0){
      buffer_was_last=true;
      vector<string>stuff=line_to_words(line.substr(line.find(":")+1),vector<char>(1,' '));
      for(int buff=0;buff<(int)stuff.size();buff++)buffer[stuff[buff].substr(0,stuff[buff].find(","))]=atof(stuff[buff].substr(stuff[buff].find(",")+1).c_str());
    }
    else if(line.find("Maximum is (")==0){
      buffer_was_last=false;
      TString word(line.substr(line.find("(")+1,line.find(",")-line.find("(")-1)), number(line.substr(line.find(",")+1,line.find(")")-line.find(",")-1));int killme=-1;
      for(int i=0;i<(int)to_do.size();i++){
	if(word==to_do[i]){killme=i;break;}
      }
      if(killme>=0){
	to_do.erase(to_do.begin()+killme,to_do.begin()+killme+1);
	done.push_back(word);daisy_chain.push_back(number.Atof());
	cout<<"MAX found "<<word<<" "<<number<<endl;
      }
    }
    else if(line.find("DAISY_CHAIN:")==0){
      buffer_was_last=false;
      cout<<"Ranking complete; no finishing necessary!  Will check kfold files..."<<endl;
      produce_kfold_files(done,overwrite,trans_DF);
      return;
    }
    else  buffer_was_last=false;
  }
  hey.close();
  cout<<"Last line was "<<line<<" prevline was "<<prevline<<endl;
  bool print_header=!((line.find("Significances for ")!=0)||(prevline.find("Significances for ")!=0&&line=="")||buffer_was_last);//sometimes printing will freeze in the middle of a ranking loop because of spike validation printing; skip printing of the header it it's already there---histogram processing will have a hard time of it otherwise
  ofstream fout(lname.c_str(),std::ofstream::app);
  cout<<"FINISHING ranking in "<<lname<<endl;
  DIR*lana=opendir(tmva_out_dir.c_str());struct dirent*archer;
  if(lana!=NULL){
    TString key=Form("TMVA%s%s-",sam_tag(-1).Data(),ftag(done).Data());
    while((archer=readdir(lana))!=NULL){
      if(TString(archer->d_name).Contains(key))remove((tmva_out_dir+archer->d_name).c_str());
    }
  }
  while(!to_do.empty()){
    rank_loop_one_var(to_do,done,daisy_chain,fout,tag,overwrite,trans_DF,print_header,(buffer_was_last?buffer:map<TString,double>()));
    print_header=true;buffer_was_last=false;
  }
  ranking_back_matter(done,daisy_chain,fout,overwrite,trans_DF);
}

void bdt_ranker::ranking_back_matter(const vector<TString>&done,const vector<double>&daisy_chain,ofstream&fout,bool overwrite,int trans_DF,bool print)const{
  int offset=20;
  if(print){
    //okay, so you don't always start ranking with just one variable
    //when looping through daisy chain (which is the whole shtick),
    //you need to start at an index = ranked.size()-1 (ranked is just the starting set)
    int startdex=ranked.size()-1;
    fout<<"Ranking Complete.  The training path looks like:"<<endl;
    fout<<"Began with: ";
    for(auto a:ranked)fout<<a<<" ";
    fout<<", S/sqrt(S+B)="<<daisy_chain[startdex]<<endl;
    fout<<setw(offset+2)<<"VARIABLE |"<<setw(offset)<<"S/sqrt(S+B)"<<endl<<"----------------------------"<<endl;
    for(int i=startdex+1;i<(int)daisy_chain.size();i++)fout<<setw(offset)<<done[i]<<" |"<<setw(offset-3)<<daisy_chain[i]<<endl;
    fout<<"DAISY_CHAIN: ";
    for(auto a:done)fout<<a<<" ";
    fout<<endl;
  }
  produce_kfold_files(done,overwrite,trans_DF);
}

void bdt_ranker::produce_kfold_files(const vector<TString>&done,bool overwrite,int trans_DF)const{
  TCanvas *can=new TCanvas(TString(tmva_out_dir+identifier),TString(tmva_out_dir+identifier),500,300);
  bdt_trainer evendb(done,overwrite,bdt_base(*this),false,2,true);
  bdt_trainer odddb(done,overwrite,bdt_base(*this),false,1,true);
  bdt_trainer even(done,overwrite,bdt_base(*this),false,2,false);
  bdt_validate(even).print_validation_bdts(can,tmva_out_dir,-1,0,identifier+"_k-fold-even");
  bdt_trainer odd(done,overwrite,bdt_base(*this),false,1,false);
  bdt_validate(odd).print_validation_bdts(can,tmva_out_dir,-1,0,identifier+"_k-fold-odd");
  delete can;
}

void bdt_ranker::rank_loop_one_var(vector<TString>&to_do,vector<TString>&done,vector<double>&daisy_chain,ofstream&fout,const string&tag,bool overwrite,int trans_DF,bool print_header,const map<TString,double>&buffer)const{
  string pdir=check_mkdir(tmva_out_dir+"bdt_spikes");
  TCanvas *can=new TCanvas(TString(pdir+tag),TString(pdir+tag),500,300);  int offset=20;
  map<double,TString>candidates;
  if(print_header){
    fout<<"Significances for ";
    for(auto a:done)fout<<a<<", ";
    fout<<"+ :"<<endl;
  }
  fout<<"BUFFER: ";
  for(const auto&it:buffer){
    fout<<it.first<<","<<it.second<<" ";
    candidates[it.second]=it.first;
  }
  for(auto var:to_do){
    if(buffer.find(var)!=buffer.end())continue;
    vector<TString>frodo(done);frodo.push_back(var);double sig=bdt_validate(frodo,overwrite,(bdt_base)*this).val_sig(-1,trans_DF).second;
    fout<<var<<","<<sig<<" ";
    candidates[sig]=var;
  }
  fout<<endl;
  ostringstream labels,values;double last_jump=0.,last_val=0.;int index=0;
  for(auto&c:candidates){
    values<<setw(offset)<<c.first;labels<<setw(offset)<<c.second;
    double new_jump=0.;bool flip=false;
    if(index==0){
      new_jump=c.first-daisy_chain.back();
    }
    else{
      new_jump=c.first-last_val;flip=(signbit(new_jump)!=signbit(last_jump));
    }
    if((abs(new_jump)>0.05 && !flip) || (to_do.size()==1)){//print the weird ones and the last one
      vector<TString>bilbo(done);bilbo.push_back(c.second);ostringstream tagger;tagger<<identifier<<"-"<<bilbo.size()<<"+"<<c.second<<(to_do.size()==1?"_FIN_":"");
      bdt_validate(bilbo,overwrite,(bdt_base)*this).print_validation_bdts(can,pdir,-1,trans_DF,tagger.str());
    }
    last_jump=new_jump;last_val=c.first;
    index++;
  }
  fout<<labels.str()<<endl<<values.str()<<endl;
  TString pushit=candidates.rbegin()->second;double val=candidates.rbegin()->first;index=-1;
  for(int i=0;i<(int)to_do.size();i++){
    if(to_do[i]==pushit)index=i;
  }
  done.push_back(pushit);daisy_chain.push_back(val);to_do.erase(to_do.begin()+index,to_do.begin()+index+1);
  fout<<endl<<"Maximum is ("<<done.back()<<","<<val<<")"<<endl;
  cout<<endl<<"Maximum is ("<<done.back()<<","<<val<<")"<<endl;
  fout<<"------------------------------------------------------------------------------------------------"<<endl;
  delete can;
}

vector<TString> bdt_ranker::print_training_plot(TCanvas *c,const TString&tag,const TString&pdir,int tDF)const{
  TStyle *style=set_style();style->cd();
  ifstream hey(log(string(tag)+tdf_str(tDF)).c_str());vector<TString> all_vars,labels;vector<string>unranked;string line;
  if(debug)cout<<"print_training_plot("<<c<<","<<tag<<","<<pdir<<")---log: "<<log(string(tag)+tdf_str(tDF))<<endl;
  vector<pair<vector<double>,vector<double> > >xy_coords;vector<int>woops;
  while(hey.good()){
    getline(hey,line);
    if(line.find("DAISY_CHAIN:")==0){
      for(auto w:line_to_words(line.substr(line.find(" ")+1)))all_vars.push_back(TString(w));
    }
    else if(line.find("RANKED VARIABLES:")==0)unranked=line_to_words(line.substr(line.find(":")+1));
    else if(line.find("Significance starts at")==0){
      xy_coords.push_back(pair<vector<double>,vector<double> >(vector<double>(1,0.),vector<double>(1,atof(line.substr(line.find("at")+4).c_str()))));
    }
    else if(line.find("Significances for")==0){
      string lv,ln;getline(hey,lv);getline(hey,ln);
      if(debug)cout<<"processing lines "<<lv<<" "<<ln<<endl;
      vector<string> words=line_to_words(lv);vector<float> numbers=line_to_floats(ln);
      map<double,string>ordered;
      for(int i=0;i<(int)words.size();i++)ordered[numbers[i]]=words[i];
      vector<double>xs,ys;
      xs.push_back(xy_coords.back().first.back());ys.push_back(xy_coords.back().second.back());
      int index=xs.back();woops.push_back(index);
      for(auto&it:ordered){
	index++;
	if(debug)cout<<it.second<<" "<<it.first<<endl;
	xs.push_back(index);ys.push_back(it.first);
	labels.push_back(TString(it.second));
      }
      xy_coords.push_back(pair<vector<double>,vector<double> >(xs,ys));
    }
  }
  TString name="training-"+tag+tdf_tstr(tDF)+".pdf";
  TMultiGraph*zihuatanejo=new TMultiGraph(name,name);
  for(int i=0;i<(int)xy_coords.size();i++){
    if(debug)cout<<"graph with "<<xy_coords[i].first.size()<<" into "<<name<<endl;
    TGraph*dasboot=new TGraph(xy_coords[i].first.size(),&xy_coords[i].first[0],&xy_coords[i].second[0]);
    dasboot->SetLineColor(color(i));
    dasboot->SetMarkerColor(color(i));
    dasboot->SetMarkerStyle(kFullDotMedium);
    zihuatanejo->Add(dasboot);
  }
  c->cd();
  zihuatanejo->SetMinimum(0.85*xy_coords.front().second.front());
  zihuatanejo->SetMaximum(1.15*xy_coords.back().second.back());
  zihuatanejo->Draw("ALP");
  TAxis *xax=zihuatanejo->GetXaxis();
  TLatex l=plot_latex();l.SetTextSize(0.015);
  xax->SetLabelSize(0.025);
  zihuatanejo->GetYaxis()->SetTitle("S/#sqrt{S+B}, VALIDATION");
  zihuatanejo->GetYaxis()->SetLabelSize(0.03);
  zihuatanejo->GetYaxis()->SetTitleSize(0.04);
  zihuatanejo->GetYaxis()->SetTitleOffset(0.85);
  for(int i=0;i<(int)labels.size();i++) xax->SetBinLabel(xax->FindBin(i*1.+1.),labels[i]);
  double min=xax->GetXmin(),max=xax->GetXmax(),lo=style->GetPadLeftMargin(),hi=1.-style->GetPadRightMargin();
  int nunranked=unranked.size();
  for(int il=0;il<(int)all_vars.size()-nunranked;il++){
    TString label,label1;label+=(nunranked+il+1);label+=" vars";label1+="(+";label1+=all_vars[nunranked+il];label1+=")";
    double coord=lo+(woops[il]-min)/(max-min)*(hi-lo);
    TLine*ligne=new TLine();ligne->SetLineWidth(1);ligne->SetLineStyle(7);ligne->DrawLineNDC(coord,style->GetPadBottomMargin(),coord,1.-style->GetPadTopMargin());
    l.SetTextColor(color(il+1));
    l.DrawLatex(coord+0.005,0.75,label);
    l.DrawLatex(coord+0.005,0.72,label1);
  }
  l.SetTextColor(kBlack);
  l.SetTextSize(0.03);ostringstream pedado;pedado<<"Rank: ";
  for(auto var:all_vars)pedado<<var.Data()<<", ";
  pedado<<"sig="<<setprecision(3)<<xy_coords.back().second.back();
  l.DrawLatex(0.08,0.95,pedado.str().c_str());
  c->Print(check_dir(pdir)+name);
  return all_vars;
}

vector<float> bdt_ranker::line_to_floats(const string&line)const{
  vector<float> floats;
  for(auto yup:line_to_words(line))floats.push_back(atof(yup.c_str()));
  return floats;
}

int bdt_ranker::find_next_whitespace(const string&line,int pos,const vector<char>whites)const{
  vector<int>nexts;
  for(auto c:whites)nexts.push_back(line.find(c,pos)==string::npos?-99:line.find(c,pos));
  return *min_element(nexts.begin(),nexts.end());
}

bool bdt_ranker::word_is_white(const string&word,const vector<char>whites)const{
  if(word=="")return false;
  for(auto c:whites){
    if(word==c)return true;
  }
  return false;
}

vector<string> bdt_ranker::line_to_words(const string&line,const vector<char>whites)const{
  int pos=0;vector<string>words;
  while(pos>=0&&pos<(int)line.length()){
    int space=find_next_whitespace(line,pos,whites);
    if(space<0){
      string new_word=line.substr(pos);
      if(!word_is_white(new_word,whites))words.push_back(new_word);
      break;
    }
    string new_word=line.substr(pos,space-pos);
    if(!word_is_white(new_word,whites)&&new_word.length()>0)words.push_back(new_word);
    pos=space+1;
  }
  return words;
}
