#ifndef BDT_RANKER_H
#define BDT_RANKER_H
#include "bdt_validate.h"

class bdt_ranker : public bdt_base{
 public:
  bdt_ranker(const string&in,const string&out,const vector<TString>&bs,const vector<TString>&vs,int nt=2,int nj=2,int eff_cut=70,int ratio=20,bool cts=true,bool lptv=true,bool btv=true,bool wmll=true,bool met=false);
  bdt_ranker(const string&in,const string&out,const vector<TString>&bs,const vector<TString>&vs,const bdt_base&base);
  void do_ranking(const string&tag,bool overwrite=false,int trans_DF=0);
  void finish_ranking(const string&tag,bool overwrite,int trans_DF=0);
  vector<TString> print_training_plot(TCanvas *c,const TString&tag,const TString&pdir,int tDF)const;
  vector<string> line_to_words(const string&line,const vector<char>whites=vector<char>(1,' '))const;
  vector<string> line_to_words(const string&line,char white)const{return line_to_words(line,vector<char>(1,white));}
  vector<float> line_to_floats(const string&line)const;

 protected:
  void ranking_back_matter(const vector<TString>&done,const vector<double>&daisy_chain,ofstream&fout,bool overwrite,int trans_DF)const;
  void produce_kfold_files(const vector<TString>&done,bool overwrite,int trans_DF)const;
  void rank_loop_one_var(vector<TString>&to_do,vector<TString>&done,vector<double>&daisy_chain,ofstream&fout,const string&tag,bool overwrite,int trans_DF,bool print_header=true)const;
  string log(const string&tag)const{return tmva_out_dir+"ranking_log_"+tag+".txt";}
  string log()const{return tmva_out_dir+"ranking_log_"+identifier+".txt";}
  int find_next_whitespace(const string&line,int pos,const vector<char>whites)const;
  bool word_is_white(const string&word,const vector<char>whites)const;
  
  vector<TString> ranked,unranked;
};

#endif
