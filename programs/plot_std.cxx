#include "std_plotter.h"
#include "TMVA/correlations.h"

int main(int argc,char**argv){
  string in="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20160309_srli/",out=in+"output/",pdir="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/plots/tmva_fun/20160309_srli/";
  int tF=0;bool scattershot=false;
  vector<int> cases={tF};//0,1,2};
  if(argc==2){
    string tag=argv[1];
    if(tag[tag.length()-1]!='/')tag+="/";
    in="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/"+tag;out=in+"output/";pdir="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/plots/tmva_fun/"+tag;
  }
  else{
    if(argc>1)in   =argv[1];
    if(argc>2)out  =argv[2];
    if(argc>3)pdir =argv[3];
    if(argc>4)tF   =atoi(argv[4]);
    vector<int> cases={tF};
  }
  std_plotter navis(in,out);vector<TString>tmva_files;
  for(auto tf:cases)tmva_files=navis.all_plots(pdir,tf);
  if(!scattershot)return 0;
  for(auto f:tmva_files){
    if(TFile::Open(f)==NULL)continue;
    bool li=f.Contains("j0_j1");
    string file=string(f);int start=file.find("TMVA",file.find_last_of("/"))+4,endtag=(file.find("-MV2")),end=(li?file.find("-j0_j1"):file.find("-mBB"));
    string ftag=file.substr(start,endtag-start);vector<string>vars=navis.line_to_words(file.substr(end,file.find(".root",end)-end),'-');vector<TString>vs;
    for(auto v:vars)vs.push_back(TString(v));
    bdt_trainer acqua(in,out,"plotting",vs,false,2,(ftag.find("3jet")==string::npos?2:3));
    ofstream fout(pdir+"tmva_corr_std_files.txt");//save files to run correlations over if this is something folks want...
    for(auto i:{-1,1,2}){
      cout<<"Correlations for "<<acqua.filename(vs,i)<<endl;
      if(TFile::Open(acqua.filename(vs,i))==NULL)acqua.train(i);
      fout<<acqua.filename(vs,i)<<endl;
      /*
      correlations(acqua.filename(vs,i));
      for(auto letter:vector<string>{"S","B"}){
	system(("epstopdf plots/CorrelationMatrix"+letter+".eps").c_str());
	system(("mv plots/CorrelationMatrix"+letter+".pdf "+pdir+"CorrelationMatrix"+letter+ftag+acqua.sam_tag(i)+(tF>0?(tF>1?"_tD":"_tF"):"")+".pdf").Data());
	}*/
    }
  }
  return 0;
}
