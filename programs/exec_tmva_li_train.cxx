#include <vector>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TColor.h"
#include "correlations.C"

#include "TSystemDirectory.h"
#include <utility>
#include "sig_plotter.h"

//Dear Lord, I'm sorry that this is horrendous stylistically...
void tmva_li_train();
bool keepspect(){return false;}

 #if not defined(__CINT__) || defined(__MAKECINT__)
 // needs to be included when makecint runs (ACLIC)
 #include "TMVA/Factory.h"
 #include "TMVA/MethodBDT.h"
 #include "TMVA/Tools.h"
 #endif


int main(int argc,char**argv){
  TString tmva_dat_dir="",tag="output";int tr_opt=2, bg=-1, through_var_i=-1;
  if(argc>1)tmva_dat_dir=TString(argv[1]);
  if(argc>2)tr_opt=atoi(argv[2]);
  if(argc>3)bg=atoi(argv[3]);
  if(argc>4)tag=TString(argv[4]);
  if(argc>5)through_var_i=atoi(argv[5]);
  
   sig_plotter parser(keepspect());
   /*
   TString tmva_dir(TString(gRootDir) + "/tmva");
   if(gSystem->Getenv("TMVASYS"))       tmva_dir = TString(gSystem->Getenv("TMVASYS"));
   gROOT->SetMacroPath(tmva_dir + "/test/:" + gROOT->GetMacroPath() );
   */
   //some variables and configurable parameters to be used later

   
   //determine variables to use, including how many of them
    TString ntag="";
    std::pair<std::vector<TString>,std::vector<TString> > vars=parser.li_vars(tr_opt);
    if(through_var_i>0&&through_var_i<vars.first.size()){
      vars.second.insert(vars.second.end(),vars.first.begin()+through_var_i,vars.first.end());
      vars.first=std::vector<TString>(vars.first.begin(),vars.first.begin()+through_var_i);
      ntag+="_";ntag+=through_var_i;
    }
    //input files directory
   if(tmva_dat_dir=="") tmva_dat_dir="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/default/";
   else if(tmva_dat_dir=="RETURN_N_VARS"){
     std::cout<<vars.first.size()<<std::endl;
     return 0;
   }
   else if(tmva_dat_dir=="RETURN_VARS"){
     for(auto ok:vars.first)std::cout<<ok<<" ";
     std::cout<<std::endl;
     return 0;
   }
   else if(string(tmva_dat_dir.Data()).substr(string(tmva_dat_dir.Data()).length()-1)!="/")tmva_dat_dir+="/";

   TString ftag=parser.li_tag(tr_opt)+parser.bg_tag(bg)+ntag;
   //decide whether to ues standard variables for training or LI variables
   //if LI variable training, use standard MET variable or LI MET candidates
   TString outfileName=tmva_dat_dir+"/"+tag+"/"+"TMVA"+ftag+".root";

   // Apply additional cuts on the signal and background samples (can be different)
   //write now the event selection is hi pTV, 2 jet, no MET
//    TCut mycuts = ((tr_opt<0||tr_opt>=parser.n_li_cases())?"pTV / 1000. > 120 && nJ == 2 && mLL / 1000. > 83 && mLL / 1000. < 99.":"pTV / 1000. > 120 && nJ == 2 && mLL / 1000. > 60 && mLL / 1000. < 110.");
   //2btag cut
   double fcbe30=0.9540,fcbe50=0.7535,fcbe60=0.4496,fcbe70=-0.0436,fcbe77=-0.4434,fcbe80=-0.5911,fcbe85=-0.7887,fcbe90=-0.9185;//short for the FixedCutBEff_x weight cuts
   double bcut=fcbe70;
   ostringstream cuts;cuts<<"pTV / 1000. < 120 && nJ == 2 && mLL / 1000. > 60 && mLL / 1000. < 110. && MV2c20B1>"<<bcut<<" && MV2c20B2>"<<bcut<<" && EventNumber%3!=0";//the last to set aside a third of events for validation
   TCut mycuts = cuts.str().c_str();
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

    // This loads the library
    TMVA::Tools::Instance();

    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
 
    // --- Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;

    // --- Here the preparation phase begins

    std::vector<TFile*> openfiles;
    std::srand(std::time(0));
    TH1F* histSig = new TH1F("Significances", "Significances", 30, .4, 1.3);
    TH1F* histPoisSig = new TH1F("New Sig", "New Sig", 30, 0, 1.5);
    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    // Create the factory object. Later you can choose the methods
    // whose performance you'd like to investigate. The factory is 
    // the only TMVA object you have to interact with
    //
    // The first argument is the base of the name of all the
    // weightfiles in the directory weight/
    //
    // The second argument is the output file for the training results
    // All TMVA output can be suppressed by removing the "!" (not) in
    // front of the "Silent" argument in the option string
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification"+ftag, outputFile,
						"!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );

    for(auto var:  vars.first)factory->AddVariable(  var,'F');
    for(auto spec:vars.second)factory->AddSpectator(spec,'F');

    // global event weights per tree (see below for setting event-wise weights)
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;

    //    TSystemDirectory dir("taggedReadProper/fetch/data-MVATree/", "taggedReadProper/fetch/data-MVATree/");
    //    TSystemDirectory dir("TaggedReadAll/fetch/data-MVATree/", "TaggedReadAll/fetch/data-MVATree/");
    TSystemDirectory dir(tmva_dat_dir,tmva_dat_dir);
    //TSystemDirectory dir("NoPileUp/fetch/data-MVATree/", "NoPileUp/fetch/data-MVATree/");
    TList* files = dir.GetListOfFiles();
    TFile* inputFile;
    for (TObject* obj : *files) {
      TString fileName = obj->GetName();
      bool signal=fileName.Contains("ZHll125.root");
      if (!signal&&!parser.is_bkg(fileName,bg)) continue;
      cout<<"Looking at file "<<fileName<<endl;
      fileName.Prepend(dir.GetName());
      inputFile = TFile::Open(fileName);
      openfiles.push_back(inputFile);
      TTree* tree = (TTree*) inputFile->Get("Nominal");
      if (tree==NULL)continue;
      if (tree->GetEntries() == 0) continue;
      if (signal) factory->AddSignalTree(tree, signalWeight);
      else factory->AddBackgroundTree(tree, backgroundWeight);
    }
//    factory->SetBackgroundWeightExpression( "EventWeight" );
   factory->SetWeightExpression("EventWeight");

   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut,
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   TString splitOpt = "SplitMode=random:!V:NormMode=None:SplitSeed=";
   splitOpt += 0;//CHANGE BACK LATER! std::rand();
   factory->PrepareTrainingAndTestTree( mycuts, splitOpt);

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   factory->BookMethod( TMVA::Types::kBDT, "BDT",
			"!H:!V:NTrees=200:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.15:SeparationType=GiniIndex:nCuts=100:PruneMethod=NoPruning:MinNodeSize=1.5%" );

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // ---- STILL EXPERIMENTAL and only implemented for BDT's ! 
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
   double sig = dynamic_cast<TMVA::MethodBDT* >(factory->GetMethod("BDT"))->GetSignificance();
   std::cout << "Significance " << sig << std::endl;
   for (TFile* f : openfiles) {
     f->Close();
   }
   std::cout << "Deleting factory..." << std::endl;
   delete factory;
   return 0;
}
