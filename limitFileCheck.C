#include <fstream>
#include <iostream>
#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TROOT.h>
#include <TKey.h>
#include <TClass.h>
#include <algorithm>
#include <functional>
#include <array>

using namespace std;

// Tool for comparing two VH limit input files. Run with:
// root -b limitFileCheck.C+

// To run on CxAODReader output, set m_isRun2 to true.
//In general need to specify reference and test files, sample to run systematic on, names for lables, etc. below.

// general configuration
bool m_isRun2 = true;  //the naming conventions changed slightly this flag should take care to catch the differences
bool m_allRegions = false; // loop on all regions found in input file
bool m_nominalOnly = true; //produces only plots for nominal if m_drawCompare is set to true
bool m_printToLogFile = false; // printing systematic tables to log files (WARNING: uses same directory as input files)
bool m_drawSingle = false; // ploting and saving all systematics for each input file
bool m_drawGroup = false; // ploting and saving all systematics in groups for each input file
bool m_drawCompare = true; // ploting and saving all systematics comparing both input files
string date = "20160705/0lep";
string m_plotDir = "/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/plots/closure/"+date+"/dtag_ttag/"; // directory for saving plots - will be automatically created
string m_graphicsExt = ".png"; // plot file format
double m_tsize = 0.05; // text size for plots. 0.05 for single, 0.07 for compare
double m_limitXup = 250; // if (m_limitXup > 0) set x-axis range to (0, m_limitXup)

// for printout
int m_width = 10; // column width
string m_testName = "test"; // table headers 
string m_refName = "ref";
double m_scale = 100; // scaling for systematic shifts (100 to have all in %)
bool m_quiet = false; // less printing

// file reading
//int m_rebin = 25; // rebin for all histograms
int m_rebin = 2; // if Run2

// thresholds for flag printing
double m_maxYieldDev = 0.05; // maximum yield deviation (0.1 = 10%)
double m_maxSysNormDev = 0.2; // maximum deviation in normalization shift from systematics
double m_treshSmall = 0.005 * m_scale; // maximum normalization shift for being 'small' (0.1 * m_scale = 10%)
double m_treshZero = 1e-9 * m_scale; // maximum normalization shift for being 'zero' (0.1 * m_scale = 10%)
double m_treshShapeSmall = 0.001 * m_scale; // maximum shape shift for being 'small' (0.1 * m_scale = 10%)
double m_treshNormAsym = 0.05; // minimum difference in normalization shift for being 'asymmetric' (0.1 = 10%)
double m_treshShapeAsym = 0.001; // maximum chi2 in shape comparison for being 'shape asymmetric'

// for reading norm files from WSMaker
bool m_normFile = false;

void setStyleAtlas() {

  Int_t icol = 0;
  gStyle->SetFrameBorderMode(icol);
  gStyle->SetFrameFillColor(icol);
  gStyle->SetCanvasBorderMode(icol);
  gStyle->SetCanvasColor(icol);
  gStyle->SetPadBorderMode(icol);
  gStyle->SetPadColor(icol);
  gStyle->SetStatColor(icol);
  //gStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  gStyle->SetPaperSize(20, 26);

  // set margin sizes
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font = 42; // Helvetica
  Double_t tsize = m_tsize; // 0.05 in ATLAS style
  gStyle->SetTextFont(font);

  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font, "x");
  gStyle->SetTitleFont(font, "x");
  gStyle->SetLabelFont(font, "y");
  gStyle->SetTitleFont(font, "y");
  gStyle->SetLabelFont(font, "z");
  gStyle->SetTitleFont(font, "z");

  gStyle->SetLabelSize(tsize, "x");
  gStyle->SetTitleSize(tsize, "x");
  gStyle->SetLabelSize(tsize, "y");
  gStyle->SetTitleSize(tsize, "y");
  gStyle->SetLabelSize(tsize, "z");
  gStyle->SetTitleSize(tsize, "z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth((Width_t) 2.);
  gStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

  // get rid of X error bars 
  //gStyle->SetErrorX(0.001);
  // get rid of error bar caps
  //gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1111); 
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}

class SysParam {
public:
  string name;
  double shiftSysUp;
  double shiftSysDo;
  double intSysCrossUp;
  double intSysCrossDo;
  double shapeChi2Ndf;
  TH1* h_nominal;
  TH1* h_sysUp;
  TH1* h_sysDo;
  TH1* h_sysUpNorm;
  TH1* h_sysDoNorm;

  SysParam() {
  }

  bool flagBothSmall() {
    return (fabs(shiftSysUp) < m_treshSmall && fabs(shiftSysDo) < m_treshSmall);
  }

  bool flagAnyZero() {
    return (fabs(shiftSysUp) < m_treshZero || fabs(shiftSysDo) < m_treshZero);
  }

  bool flagSameSign() {
    return (shiftSysUp * shiftSysDo > 0);
  }

  bool flagShapeSmall() {
    return (intSysCrossUp < m_treshShapeSmall && intSysCrossDo < m_treshShapeSmall);
  }

  bool flagNormAsym() {
    return (fabs(shiftSysUp / shiftSysDo + 1) > m_treshNormAsym);
  }

  bool flagShapeAsym() {
    return (shapeChi2Ndf > m_treshShapeAsym);
  }
};

class Yield {
public:
  string sample;
  double integral;
  double err;

  Yield() {
  }
};

bool createOutputDir(){
  Long_t flags = 0;
  gSystem->GetPathInfo(m_plotDir.c_str(), (Long_t*)0, (Long_t*)0, &flags, (Long_t*)0);
  if (flags & 2) {
    // directory already exists
    cout << "Output directory for plots "<< m_plotDir << " already exists." << endl;
    return true;
  } 
  //create directory
  if (0 == gSystem->mkdir(m_plotDir.c_str(),true)){
    cout << "Created output directory for plots: "<<m_plotDir << endl;
    return true;
  }
  else{
    cout << "Failed to create output directory "<<m_plotDir<<" for plots." << endl;
    return false;
  }
}

void setEmptyBinsToOne(TH1* ratio){
  unsigned int nBins = ratio->GetNbinsX();
  for (unsigned int i = 0; i < nBins+2; i++) {
    if(!ratio->GetBinContent(i)){
      ratio->SetBinContent(i,1);
    }
  }
}

TH1* getHisto(TDirectory* source, string histoName, string subFolder = "", string systematic = "", string direction = "") {

 if(direction == ""){
    if(m_isRun2) direction = "1up";
    else direction = "Up";
  }

  // vpt integration

  TString histNameTString = histoName.c_str();
  if (histNameTString.Contains("allPTV") || histNameTString.Contains("highPTV")) {

    vector<string> regions;
    regions.clear();
    if (histNameTString.Contains("allPTV")) {
      regions.push_back("vpt0_90");
      regions.push_back("vpt90_120");
    }
    regions.push_back("vpt120_160");
    regions.push_back("vpt160_200");
    regions.push_back("vpt200");

    TH1* histSum = 0;

    for (unsigned int j = 0; j < regions.size(); j++) {

      TString histNameTemp = histNameTString;
      histNameTemp.ReplaceAll("allPTV", regions[j].c_str());
      histNameTemp.ReplaceAll("highPTV", regions[j].c_str());
      if (histSum == 0) {
        histSum = getHisto(source, histNameTemp.Data(), subFolder, systematic, direction);
      } else {
        histSum -> Add(getHisto(source, histNameTemp.Data(), subFolder, systematic, direction));
      }
    }
    if (histSum)
      histSum -> SetName(histNameTString);
    return histSum; 
  }

  // actual get histogram

  source -> cd();
  if (subFolder != "")
    source -> cd(subFolder.c_str());
  string name = histoName;
  TString sysName = systematic;
  if (m_isRun2) sysName.ReplaceAll("1up", direction);
  else          sysName.ReplaceAll("Up", direction);
  
  if (systematic != "" && systematic != "Nominal") {
    if(subFolder==""){
      if(!source->cd(sysName)){
	cout<<"I'm called"<<endl;
	return 0;//systematics are automatically stored in folder with their own name; if not found return null
      }
    }
    if (!m_normFile)
      name += "_" + sysName;
    gDirectory -> cd("Systematics");
  }

  
  if(sysName.Contains("JER"))cout << "getting " << name << " for systematic " << sysName << endl;

  TObject* histObj = gDirectory -> Get(name.c_str());
  if (histObj == 0) {
    if (!m_quiet)
      cout << "WARNING: histogram '" << name.c_str() << "' not found!" << endl;
    return 0;
  }

  TH1* hist = (TH1*) histObj -> Clone();
  delete histObj;

  hist -> SetDirectory(0);
  hist -> SetName(name.c_str());
  hist -> SetLineWidth(2);
  hist -> SetMarkerSize(0);
  hist -> SetFillStyle(0);
  hist -> GetXaxis() -> SetTitleSize(m_tsize);
  hist -> GetYaxis() -> SetTitleSize(m_tsize);
  hist -> GetXaxis() -> SetTitleOffset(1.0);
  hist -> GetYaxis() -> SetTitleOffset(0.8);
  hist -> GetXaxis() -> SetLabelSize(m_tsize);
  hist -> GetYaxis() -> SetLabelSize(m_tsize);


  // hack for MVA to get same binning in mBB as cut-based
  //  int nbins = hist -> GetNbinsX();
  //  hist -> SetBins(nbins, 0, 500);
  //  for (int i = 0; i <= nbins + 1; i++) {
  //    double cont  = ((TH1*) histObj) -> GetBinContent(i);
  //    double err  = ((TH1*) histObj) -> GetBinError(i);
  //    hist -> SetBinContent(i, cont);
  //    hist -> SetBinError(i, err);
  //  }
  //  hist -> Rebin(nbins / 20);
  
  if(name.find("NSigJets") == string::npos && name.find("NFwdJets") == string::npos){
    hist -> Rebin(m_rebin);
  }
  
  // set upper x axis limit
  if (m_limitXup > 0) {
    double binWidth = hist -> GetBinWidth(1);
    double xMin = hist -> GetBinLowEdge(1);
    hist -> GetXaxis() -> SetRangeUser(xMin, m_limitXup - 0.5 * binWidth);
  }

  return hist;
}

bool pushUnique(vector<string> &vec, string &str) {
  bool found = false;
  for (unsigned int i = 0; i < vec.size(); i++)
    found |= (vec[i] == str);
  if (!found) vec.push_back(str);
  return found;
}

vector<string> getListOfSystematics(TDirectory* source, string subFolder = "") {
  
  source -> cd();
  if (subFolder != "")
    source -> cd(subFolder.c_str());
  source -> cd("Systematics") ;
  vector<string> systematics;
  systematics.push_back("Nominal");
  TKey* key;
  TIter nextkey(gDirectory -> GetListOfKeys());
  while ((key = (TKey*) nextkey())) {
    TClass* cl = gROOT -> GetClass(key -> GetClassName());
    if (!cl) continue;
    //  if (cl -> InheritsFrom(TDirectory::Class())) {
      TString sysName = key -> GetName();
      // TODO this works for one sided systemtics if only up is present,
      // what about only down?
      /* if(sysName.BeginsWith("Sys") && !sysName.Contains("JER") && !(sysName.EndsWith("Do") || sysName.EndsWith("1down"))) systematics.push_back(sysName.Data());// JER is buggy right now// this was for before putting all systematics in one directory
	 else */ 
      if (!m_nominalOnly && sysName.Contains("_Sys")
              && !(sysName.EndsWith("Do") || sysName.EndsWith("1down"))) {
        //cout << "Found systematic: " << sysName.Data() << endl;
  	 sysName.Remove(0,sysName.Index("_Sys")+1);
	systematics.push_back(sysName.Data());
      }
      //    }
  }
  std::sort( systematics.begin(), systematics.end() );
  systematics.erase( std::unique( systematics.begin(), systematics.end() ), systematics.end() );
  cout << "Found " << systematics.size() << " systematics (including Nominal)." << endl;
  return systematics;
}

vector<string> getListOfRegions(TDirectory* source, string subFolder = "", string sample = "data") {
  source -> cd();
  if (subFolder != "")
    source -> cd(subFolder.c_str());
  vector<string> histoList;
  TKey* key;
  TIter nextkey(gDirectory -> GetListOfKeys());
  while ((key = (TKey*) nextkey())) {
    TClass* cl = gROOT -> GetClass(key -> GetClassName());
    if (!cl) continue;
    if (cl -> InheritsFrom(TH1::Class())){
      TString histName = key -> GetName();
      if(!histName.BeginsWith((sample + "_").c_str())) continue;
      histName.Remove(0, sample.size()+1);
      cout << "Found region: "<< histName.Data() << endl;
      histoList.push_back(histName.Data());
    }
  }
  cout << "Found " << histoList.size() << " regions." << endl;
  return histoList;
}

vector<string> getListOfSamples(TDirectory* source, string subFolder = "") {
  source -> cd();
  if (subFolder != "")
    source -> cd(subFolder.c_str());
  vector<string> samples;
  TKey* key;
  TIter nextkey(gDirectory -> GetListOfKeys());
  while ((key = (TKey*) nextkey())) {
    TClass* cl = gROOT -> GetClass(key -> GetClassName());
    if (!cl) continue;
    if (cl -> InheritsFrom(TH1::Class())){
      TString histName = key -> GetName();
      int length = histName.First("_");
      if (length < 1) continue;
      histName.Remove(length, histName.Length());
      if (histName == "EventLoop") continue;
      string histNameStr = histName.Data();
      bool found = pushUnique(samples, histNameStr);
      if (!found) {
        cout << "Found sample: "<< histName.Data() << endl;
      }
    }
  }
  cout << "Found " << samples.size() << " samples." << endl;
  return samples;
}

void scaleAndShiftHist(TH1* hist, double delta, double minHist) {
  double maxNom = hist -> GetMaximum();
  hist -> Scale(delta / maxNom);
  for (int iBin = 0; iBin <= hist -> GetNbinsX() + 1; iBin++) {
    double cont = hist -> GetBinContent(iBin);
    cont += minHist;
    hist -> SetBinContent(iBin, cont);
  }
}

double integralAndSysToShift(TH1* h_sys, TH1* h_nominal, bool normalize = false, double scale = 100) {

  double intSys = h_sys -> Integral();
  int zeroCross = 1;
  if (normalize) {
    double intNom = h_nominal -> Integral();
    h_sys -> Scale(intNom / intSys);
  }
  h_sys -> Add(h_nominal, -1);
  h_sys -> Divide(h_nominal);
  h_sys -> Scale(scale);

  if (normalize) {
    intSys = 0;
    zeroCross = 0;
    double prevCont = 0;
    for (int iBin = 1; iBin < h_sys -> GetNbinsX(); iBin++) {
      double cont = h_sys -> GetBinContent(iBin);
      intSys += TMath::Abs(cont);
      if (prevCont * cont < 0)
        zeroCross++;
      prevCont = cont;
    }
    zeroCross *= h_sys -> GetNbinsX();
    if (zeroCross == 0) {
      intSys = 0;
      zeroCross = 1;
    }
  }
  return intSys / zeroCross;
}

double computeChi2Ndf(TH1* hist1, TH1* hist2) {

  double chi2 = 0;
  int ndf = 0;

  for (int i = 0; i <= hist1 -> GetNbinsX() + 1; ++i) {

    double nData = hist1->GetBinContent(i);
    double nResult = hist2->GetBinContent(i);
    double dataErr = hist1->GetBinError(i);
    double predErr = hist2->GetBinError(i);
    if (nResult > 0 && nData > 0) {
      chi2 += pow((nData - nResult) / sqrt(dataErr * dataErr + predErr * predErr), 2);
      ndf++;
    }
  }
  return chi2 / ndf;
}

double getShapeChi2Ndf(TH1* h_sysUpNorm, TH1* h_sysDoNorm) {
  h_sysDoNorm -> Scale(-1);
  double Chi2Ndf = computeChi2Ndf(h_sysDoNorm, h_sysUpNorm);
  h_sysDoNorm -> Scale(-1);
  if (Chi2Ndf == Chi2Ndf)
    return Chi2Ndf;
  return 0;
}

string formatString(double number) {
  return TString::Format("%.4f", number).Data();
}

void printFlag(string flag, bool print) {
  if (print)
    cout << setw(m_width) << flag;
  else
    cout << setw(m_width) << "-";
}

TFile* m_file = 0;
string m_prevFileName = "";

TFile* getFile(string fileName) {
  if (m_prevFileName != fileName) {
    if (m_file) {
      m_file -> Close();
      delete m_file;
    }
    m_file = 0;
    m_prevFileName = fileName;
  }
  if (!m_file)
    m_file = TFile::Open(fileName.c_str(), "READ");
  return m_file;
}

bool sampleAvailable(string fileName, string sample, string region, string subFolder) {
  TFile* file = getFile(fileName);
  string histoName = sample + "_" + region;
  TH1* h_nominal = getHisto(file, histoName, subFolder);
  if (!h_nominal)
    return false;
  delete h_nominal;
  return true;
}

bool sampleZero(string fileName, string sample, string region, string subFolder) {
  TFile* file = getFile(fileName);
  string histoName = sample + "_" + region;
  TH1* h_nominal = getHisto(file, histoName, subFolder);
  if (!h_nominal)
    return true;
  int nBins = h_nominal -> GetNbinsX();
  double err;
  double integral = h_nominal -> IntegralAndError(0, nBins + 1, err);
  bool zero = (integral == 0);
  delete h_nominal;
  return zero;
}

bool sampleLowStats(string fileName, string sample, string region, string subFolder) {
  TFile* file = getFile(fileName);
  string histoName = sample + "_" + region;
  TH1* h_nominal = getHisto(file, histoName, subFolder);
  if (!h_nominal)
    return true;
  int nBins = h_nominal -> GetNbinsX();
  double err;
  double integral = h_nominal -> IntegralAndError(0, nBins + 1, err);
  // if less than Nmin raw events then lowStats
  double Nmin = 100;
  bool lowStats = (fabs(err / integral) > 1 / sqrt(Nmin)) || integral == 0;
  delete h_nominal;
  return lowStats;
}

vector<Yield> readYields(string fileName, vector<string> samples, string region, string subFolder) {

  vector<Yield> yields;
  TFile* file = getFile(fileName);

  for (unsigned int i = 0; i < samples.size(); i++) {
    string histoName = samples[i] + "_" + region;
    TH1* h_nominal = getHisto(file, histoName, subFolder);
    if (h_nominal == 0)
      continue;
    int nBins = h_nominal -> GetNbinsX();
    double err;
    double integral = h_nominal -> IntegralAndError(0, nBins + 1, err);
    delete h_nominal;

    Yield yield;
    yield.sample = samples[i];
    yield.integral = integral;
    yield.err = err;
    yields.push_back(yield);
  }
  return yields;
}

vector<SysParam*> readSysFromLimitFile(string fileName, string sample, string region, string subFolder) {

  string histoName = sample + "_" + region;
  if (m_normFile)
    histoName = sample;

  TFile* file = getFile(fileName);
  ofstream logFile;
  if (m_printToLogFile) {
    string logFileName = fileName + "." + histoName + ".txt";
    cout << "Printing to log file: " << logFileName << endl;
    logFile.open(logFileName.c_str());
  } else {
    if (!m_quiet)
      cout << "Printing to log file is OFF!" << endl;
  }
  
  vector<string> systematics = getListOfSystematics(file, subFolder);
  // vector<string> systematics = {"JET"};
  vector<SysParam*> sysParams;
  sysParams.clear();

  TH1* h_nominal = getHisto(file, histoName, subFolder);

  if (h_nominal == 0) {
    return sysParams;
  }

  double intNom = h_nominal -> Integral();

  //  h_nominal -> SetLineStyle(2);
  //h_nominal -> SetLineWidth(4);
  h_nominal -> SetLineColor(1);


  if (m_drawSingle || m_drawGroup || m_drawCompare) {
    cout << "Drawing is ON!" << endl;
  } else {
    if (!m_quiet)
      cout << "Drawing is OFF!" << endl;
  }

  for (unsigned int i = 0; i < systematics.size(); i++) {

    string sysUpName = "Up";
    string sysDownName = "Do";
    if(m_isRun2){
      sysUpName = "1up";
      sysDownName = "1down";
    }

    TH1* h_sysUp = getHisto(file, histoName, subFolder, systematics[i], sysUpName);
    TH1* h_sysDo = getHisto(file, histoName, subFolder, systematics[i], sysDownName);
    TH1* h_sysUpNorm = getHisto(file, histoName, subFolder, systematics[i], sysUpName);
    TH1* h_sysDoNorm = getHisto(file, histoName, subFolder, systematics[i], sysDownName);

    // for one sided systematics
    if (h_sysUp != 0 && h_sysDo == 0) {
      h_sysDo = (TH1*)h_sysUp->Clone();
      h_sysDoNorm = (TH1*)h_sysUpNorm->Clone();
      string nom=h_sysDo->GetName(),nnom=h_sysDoNorm->GetName();
      cout<<"So I'm called "<<h_sysDo->GetName()<<" "<<h_sysDoNorm->GetName()<<endl;
    }

    if (h_sysUp == 0 || h_sysDo == 0)
      continue;

    double intSysUp = integralAndSysToShift(h_sysUp, h_nominal, false, m_scale);
    double intSysDo = integralAndSysToShift(h_sysDo, h_nominal, false, m_scale);
    double intSysCrossUp = integralAndSysToShift(h_sysUpNorm, h_nominal, true, m_scale);
    double intSysCrossDo = integralAndSysToShift(h_sysDoNorm, h_nominal, true, m_scale);

    double shiftSysUp = m_scale * (intSysUp / intNom - 1);
    double shiftSysDo = m_scale * (intSysDo / intNom - 1);
    if (intNom == 0) {
      shiftSysUp = 0;
      shiftSysDo = 0;
    }
    double shapeChi2Ndf = getShapeChi2Ndf(h_sysUpNorm, h_sysDoNorm);

    SysParam* sysParam = new SysParam();
    sysParam -> name = systematics[i];
    sysParam -> shiftSysUp = shiftSysUp;
    sysParam -> shiftSysDo = shiftSysDo;
    sysParam -> intSysCrossUp = intSysCrossUp;
    sysParam -> intSysCrossDo = intSysCrossDo;
    sysParam -> shapeChi2Ndf = shapeChi2Ndf;

    sysParams.push_back(sysParam);

    if (m_printToLogFile) {
      logFile << setw(2 * m_width) << systematics[i];
      logFile << setw(2 * m_width) << shiftSysUp;
      logFile << setw(2 * m_width) << shiftSysDo;
      logFile << setw(2 * m_width) << intSysCrossUp;
      logFile << setw(2 * m_width) << intSysCrossDo;
      logFile << setw(2 * m_width) << shapeChi2Ndf;
      logFile << endl;
    }

    if (m_drawSingle || m_drawGroup || m_drawCompare) {
      h_sysUp -> SetLineColor(kRed);
      h_sysDo -> SetLineColor(kAzure - 2);
      h_sysUpNorm -> SetLineColor(kRed);
      h_sysDoNorm -> SetLineColor(kAzure - 2);
      sysParam -> h_nominal = h_nominal;
      sysParam -> h_sysUp = h_sysUp;
      sysParam -> h_sysDo = h_sysDo;
      sysParam -> h_sysUpNorm = h_sysUpNorm;
      sysParam -> h_sysDoNorm = h_sysDoNorm;
    } else {
      delete h_sysUp;
      delete h_sysDo;
      delete h_sysUpNorm;
      delete h_sysDoNorm;
    }
  }

  if (!(m_drawSingle || m_drawGroup || m_drawCompare)) {
    delete h_nominal;
  }
  if (m_printToLogFile) {
    logFile.close();
  }
  return sysParams;
}

vector<SysParam> readSysFromLogFile(string fileName, string sample, string region) {

  string histoName = sample + "_" + region;
  string logFileName = fileName + "." + histoName + ".txt";

  ifstream logFile;
  logFile.open(logFileName.c_str());

  vector<SysParam> sysParams;
  sysParams.clear();

  while (!logFile.eof()) {
    SysParam sysParam;
    logFile >> sysParam.name
            >> sysParam.shiftSysUp
            >> sysParam.shiftSysDo
            >> sysParam.intSysCrossUp
            >> sysParam.intSysCrossDo
            >> sysParam.shapeChi2Ndf;
    if (sysParam.name != "")
      sysParams.push_back(sysParam);
  }
  logFile.close();

  return sysParams;
}

void printSysParam(vector<SysParam*> sysParams) {

  cout << setw(4 * m_width) << "Systematic";
  cout << setw(m_width) << "NormUp[%]";
  cout << setw(m_width) << "NormDo[%]";
  cout << setw(m_width) << "ShapUp[%]";
  cout << setw(m_width) << "ShapDo[%]";
  cout << setw(m_width) << "ShapeChi2";
  cout << endl;

  for (unsigned int i = 0; i < sysParams.size(); i++) {

    if (!sysParams[i])
      continue;

    SysParam sysParam = *sysParams[i];

    cout << setw(4 * m_width) << sysParam.name;
    cout << setw(m_width) << formatString(sysParam.shiftSysUp);
    cout << setw(m_width) << formatString(sysParam.shiftSysDo);
    cout << setw(m_width) << formatString(sysParam.intSysCrossUp);
    cout << setw(m_width) << formatString(sysParam.intSysCrossDo);
    cout << setw(m_width) << formatString(sysParam.shapeChi2Ndf);

    printFlag("BothSmall", sysParam.flagBothSmall());
    printFlag("AnyZero", sysParam.flagAnyZero());
    printFlag("SameSign", sysParam.flagSameSign());
    printFlag("ShapSmall", sysParam.flagShapeSmall());
    printFlag("NormAsym", sysParam.flagNormAsym());
    printFlag("ShapeAsym", sysParam.flagShapeAsym());
    cout << endl;
  }
}

int findSys(string name, vector<SysParam*> sysParams) {

  for (unsigned int i = 0; i < sysParams.size(); i++) {
    if (sysParams[i] -> name == name)
      return i;
  }
  return -1;
}

void checkForMissingSys(vector<SysParam*> sysParamsTest, vector<SysParam*> sysParamsRef) {

  for (unsigned int i = 0; i < sysParamsRef.size(); i++) {
    string name = sysParamsRef[i] -> name;
    int index = findSys(name, sysParamsTest);
    if (index < 0) {
      cout << "  " << name << endl;
    }
  }
}

bool getRatioFlag(double ratio, double maxDev) {

  return ratio < (1 - maxDev) || ratio > 1 / (1 - maxDev);
}

void printCompare(double test, double ref, bool print) {
  if (print || true) {
    cout << setw(m_width) << formatString(test);
    cout << setw(m_width) << formatString(ref);
  } else {
    cout << setw(2 * m_width) << " ";
  }
}

void compareYields(vector<Yield> yieldsTest, vector<Yield> yieldsRef) {

  cout << setw(m_width) << " ";
  cout << setw((int) (1.5 * m_width)) << m_testName;
  cout << setw(2 * m_width) << m_refName;
  cout << endl;

  cout << setw(m_width) << "sample";
  cout << setw(m_width) << "yield";
  cout << setw(m_width) << "unc.";
  cout << setw(m_width) << "yield";
  cout << setw(m_width) << "unc.";
  cout << setw(m_width) << "ratio";
  cout << setw(m_width) << "unc.";
  cout << endl;

  for (unsigned int i = 0; i < yieldsTest.size(); i++) {

    cout << setw(m_width) << yieldsTest[i].sample;
    cout << setw(m_width) << TString::Format("%.2f", yieldsTest[i].integral).Data();
    cout << setw(m_width) << TString::Format("%.2f", yieldsTest[i].err).Data();

    bool found = false;
    for (unsigned int j = 0; j < yieldsRef.size(); j++) {
      if (yieldsTest[i].sample == yieldsRef[j].sample) {
        double ratio = yieldsTest[i].integral / yieldsRef[j].integral;
        double relErrTest = yieldsTest[i].err / yieldsTest[i].integral;
        double relErrRef = yieldsRef[j].err / yieldsRef[j].integral;
        double ratioErr = ratio * sqrt(pow(relErrTest, 2) + pow(relErrRef, 2));
        cout << setw(m_width) << TString::Format("%.2f", yieldsRef[j].integral).Data();
        cout << setw(m_width) << TString::Format("%.2f", yieldsRef[j].err).Data();
        cout << setw(m_width) << TString::Format("%.3f", ratio).Data();
        cout << setw(m_width) << TString::Format("%.3f", ratioErr).Data();
        bool flagLargeDev = getRatioFlag(ratio, m_maxYieldDev);
        printFlag("LargeDev", flagLargeDev);
        found = true;
      }
    }
    if (!found) {
      cout << setw(m_width) << "-";
      cout << setw(m_width) << "-";
      cout << setw(m_width) << "-";
      cout << setw(m_width) << "-";
      cout << setw(m_width) << "-";
    }
    cout << endl;
  }
  for (unsigned int j = 0; j < yieldsRef.size(); j++) {
    bool found = false;
    for (unsigned int i = 0; i < yieldsTest.size(); i++) {
      if (yieldsTest[i].sample == yieldsRef[j].sample) {
        found = true;
      }
    }
    if (!found) {
      cout << setw(m_width) << yieldsRef[j].sample;
      cout << setw(m_width) << "-";
      cout << setw(m_width) << "-";
      cout << setw(m_width) << TString::Format("%.2f", yieldsRef[j].integral).Data();
      cout << setw(m_width) << TString::Format("%.2f", yieldsRef[j].err).Data();
      cout << setw(m_width) << "-";
      cout << setw(m_width) << "-";
      cout << setw(m_width) << "-";
      cout << endl;
    }
  }
}

void compareSys(vector<SysParam*> sysParamsTest, vector<SysParam*> sysParamsRef) {

  cout << "Systematics with difference in normalization shift > " << m_maxSysNormDev * 100 << "%" << endl;

  cout << setw(2 * m_width) << " ";
  cout << setw((int) (1.5 * m_width)) << "NormUp[%]";
  cout << setw(2 * m_width) << "NormDo[%]";
  //  cout << setw(2 * m_width) << "ShapUp[%]";
  //  cout << setw(2 * m_width) << "ShapDo[%]";
  //  cout << setw(2 * m_width) << "ShapeChi2";
  cout << endl;

  cout << setw(2 * m_width) << "Systematic";
  cout << setw(m_width) << m_testName << setw(m_width) << m_refName;
  cout << setw(m_width) << m_testName << setw(m_width) << m_refName;
  //  cout << setw(m_width) << m_testName << setw(m_width) << m_refName;
  //  cout << setw(m_width) << m_testName << setw(m_width) << m_refName;
  //  cout << setw(m_width) << m_testName << setw(m_width) << m_refName;
  cout << endl;

  for (unsigned int i = 0; i < sysParamsTest.size(); i++) {
    string name = sysParamsTest[i] -> name;
    int index = findSys(name, sysParamsRef);
    if (index >= 0) {
      SysParam test = *sysParamsTest[i];
      SysParam ref = *sysParamsRef[index];


      bool flagShiftSysUp = getRatioFlag(test.shiftSysUp / ref.shiftSysUp, m_maxSysNormDev);
      bool flagShiftSysDo = getRatioFlag(test.shiftSysDo / ref.shiftSysDo, m_maxSysNormDev);
      //      bool flagIntSysCrossUp = getRatioFlag(test.intSysCrossUp / ref.intSysCrossUp, maxDev);
      //      bool flagIntSysCrossDo = getRatioFlag(test.intSysCrossDo / ref.intSysCrossDo, maxDev);
      //      bool flagShapeChi2Ndf = getRatioFlag(test.shapeChi2Ndf / ref.shapeChi2Ndf, maxDev);

      bool anyFlag = flagShiftSysUp || flagShiftSysDo; // || flagIntSysCrossUp || flagIntSysCrossDo || flagShapeChi2Ndf;

      if (anyFlag) {
        cout << setw(2 * m_width) << name;
        printCompare(test.shiftSysUp, ref.shiftSysUp, flagShiftSysUp);
        printCompare(test.shiftSysDo, ref.shiftSysDo, flagShiftSysDo);
        //        printCompare(test.intSysCrossUp, ref.intSysCrossUp, flagIntSysCrossUp);
        //        printCompare(test.intSysCrossDo, ref.intSysCrossDo, flagIntSysCrossDo);
        //        printCompare(test.shapeChi2Ndf, ref.shapeChi2Ndf, flagShapeChi2Ndf);
        cout << endl;
      }
    }
  }
}

void drawHistograms(vector<TH1*> histos, bool drawErrorBars, double limitYrange = -1) {

  double max = -999.;
  double min = 999.;

  //if ratio plot
  bool isRatio = false;
  if(histos.size() == 1) isRatio = true;

  for (unsigned int i = 0; i < histos.size(); i++) {
    if (i == 0)
      histos[i] -> Draw("hist");
    else
      histos[i] -> Draw("hist same");
    if (drawErrorBars) {
      histos[i] -> Draw("same");
    }
    double max0 = histos[i] -> GetBinContent(histos[i] -> GetMaximumBin());
    double min0 = histos[i] -> GetBinContent(histos[i] -> GetMinimumBin());
    if (max < max0)
      max = max0;
    if (min > min0)
      min = min0;
  }
  double factor = 0.2;
  double delta = max - min;
  if (max != 0)
    max += factor * delta;
  if (min != 0)
    min -= factor * delta;
  if (limitYrange > 0 && max > limitYrange)
    max = limitYrange;
  if (limitYrange > 0 && -min > limitYrange)
    min = -limitYrange;
 
  //to account for ratio plots with straight line at one...
  if(!delta){
    max += factor * 1.;
    min -= factor * 1.;
  }

  histos[0] -> GetYaxis() -> SetRangeUser(min, max);

  //add line to ratio plot
  if(isRatio) {
    uint nBins = histos[0] ->GetNbinsX();
    TH1* line = (TH1*)histos[0]->Clone("line");
    for (unsigned int i = 0; i < nBins+2; i++){
      line->SetBinContent(i,1.);
      line->SetBinError(i,0.);
    }
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->Draw("same HIST");
    gPad->RedrawAxis();
  }

}

void plotSys(vector<SysParam*> sysParams, string fileTag) {

  m_tsize = 0.05;
  setStyleAtlas();
  gStyle->SetPadRightMargin(0.16);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(0);

  TCanvas * canvas = new TCanvas("canvas", "canvas", 600, 400);

  for (unsigned int i = 0; i < sysParams.size(); i++) {

    SysParam param = *sysParams[i];

    param.h_sysUpNorm -> SetLineStyle(2);
    param.h_sysDoNorm -> SetLineStyle(2);
    param.h_sysUp -> Draw("hist");
    param.h_sysDo -> Draw("hist same");
    param.h_sysUpNorm -> Draw("hist same");
    param.h_sysDoNorm -> Draw("hist same");


    double maxFactor = 1.1;
    double maxNom = param.h_nominal -> GetMaximum();
    double maxUp = param.h_sysUp -> GetMaximum();
    double maxDo = param.h_sysDo -> GetMaximum();
    double minUp = param.h_sysUp -> GetMinimum();
    double minDo = param.h_sysDo -> GetMinimum();
    double max = maxUp;
    double min = minUp;
    if (max < maxDo)
      max = maxDo;
    if (min > minDo)
      min = minDo;
    double delta = max - min;
    double maxHist = min + maxFactor * delta;
    double minHist = max - maxFactor * delta;
    param.h_sysUp -> GetYaxis() -> SetRangeUser(minHist, maxHist);
    param.h_sysUp -> GetYaxis() -> SetTitle("Shift [%]");
    TH1* h_nomScaled = (TH1*) param.h_nominal -> Clone();
    scaleAndShiftHist(h_nomScaled, delta, minHist);
    h_nomScaled -> Draw("hist same");

    TLatex* latex = new TLatex();
    latex -> DrawTextNDC(0.3, 0.85, param.name.c_str());
    latex -> DrawTextNDC(0.3, 0.77, param.h_nominal -> GetName());

    double maxX = param.h_nominal -> GetXaxis() -> GetXmax();

    TGaxis* axis = new TGaxis(maxX, minHist, maxX, maxHist, 0, (2 * maxFactor - 1) * maxNom, 510, "+L");
    axis -> SetLabelSize(0.05);
    axis -> SetLabelFont(42);
    axis -> SetTitleSize(0.05);
    axis -> SetTitleFont(42);
    axis -> SetTitle("Events");
    axis -> SetTitleOffset(1.5);
    axis -> Draw();

    TLegend* leg = new TLegend(0.7, 0.75, 0.8, 0.9);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(m_tsize);
    leg -> AddEntry(param.h_sysUp, "Up", "l");
    leg -> AddEntry(param.h_sysDo, "Down", "l");
    leg -> Draw();

    string plotFileName = m_plotDir + fileTag + "_" + param.h_nominal -> GetName() + "_" + param.name + m_graphicsExt;
    canvas -> Print(plotFileName.c_str());
  }
  delete canvas;
}

bool testGroup(string group, TString sysName) {
  bool isGroup = false;
  if (group == "SysLepton") {
    isGroup |= sysName.Contains("SysElec");
    isGroup |= sysName.Contains("SysMuon");
  } else if (group == "SysMbb") {
    isGroup |= sysName.Contains("Mbb");
  } else if (group == "SysDPhi") {
    isGroup |= sysName.Contains("DPhi");
  } else if (group == "SysPt") {
    isGroup |= sysName.Contains("Pt") && !sysName.Contains("SysJet");
  } else if (group == "SysModel") {
    isGroup |= sysName.Contains("Mbb");
    isGroup |= sysName.Contains("DPhi");
    isGroup |= sysName.Contains("Pt") && !sysName.Contains("SysJet");
  } else {
    isGroup |= sysName.Contains(group.c_str());
  }
  return isGroup;
}

void plotGroupSys(vector<SysParam*> sysParamsAll, string fileTag) {

  m_tsize = 0.05;
  setStyleAtlas();
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.4);
  
  vector<string> groups;
  groups.push_back("SysBTagB");
  groups.push_back("SysBTagC");
  groups.push_back("SysBTagL");
  groups.push_back("SysJet");
//  groups.push_back("SysLepton"); // -> other
//  groups.push_back("SysMbb"); // -> model
//  groups.push_back("SysPt"); // -> model
//  groups.push_back("SysDPhi"); // -> model
  groups.push_back("SysModel");
  groups.push_back("SysOther");

  for (unsigned int iG = 0; iG < groups.size(); iG++) {

    vector<SysParam* > sysParams;

    for (unsigned int i = 0; i < sysParamsAll.size(); i++) {
      TString name = sysParamsAll[i] -> name;
      bool isGroup = testGroup(groups[iG], name);
      if (groups[iG] == "SysOther") {
        bool isKnown = false;
        for (unsigned int jG = 0; jG < groups.size(); jG++) {
          if (iG != jG) {
            isKnown |= testGroup(groups[jG], name);
          }
        }
        isGroup = !isKnown;
      }
      if (isGroup)
        sysParams.push_back(sysParamsAll[i]);
    }

    if (sysParams.size() == 0)
      continue;

    const unsigned int N = sysParams.size();
    Int_t MyPalette[N];
    Double_t r[] = {0.0, 0.0, 1.0, 1.0, 0.9, 0.0};
    Double_t g[] = {1.0, 0.0, 0.0, 0.0, 0.9, 1.0};
    Double_t b[] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    Double_t stop[] = {0., .20, .40, .60, .80, 1.0};
    Int_t FI = TColor::CreateGradientColorTable(6, stop, r, g, b, N);
    for (unsigned int i = 0; i < N; i++) MyPalette[i] = FI + i;
    gStyle->SetPalette(N, MyPalette);

    vector<string> dirs;
    dirs.push_back("Up");
    dirs.push_back("Do");

    for (unsigned int iD = 0; iD < dirs.size(); iD++) {

      TCanvas * canvas = new TCanvas("canvas", "canvas", 1200, 600);
      TH1* firstHisto = 0;
      double max = 0;
      double min = 0;
      TLegend* leg = new TLegend(0.65, 0.01, 0.99, 0.99);
      leg -> SetFillStyle(0);
      leg -> SetBorderSize(0);
      leg -> SetTextSize(m_tsize);
      TLatex* latex = new TLatex();

      string region = "";

      for (unsigned int i = 0; i < N; i++) {

        SysParam param = *sysParams[i];
        region = param.h_nominal -> GetName();
        TH1* histo = param.h_sysUp;
        if (dirs[iD] == "Do")
          histo = param.h_sysDo;
        if (i == 0) {
          histo -> Draw("hist");
          firstHisto = histo;
        } else {
          histo -> Draw("hist same");
        }
        //      param.h_sysDo -> Draw("hist same");

        histo -> SetLineColor(MyPalette[i]);
        //      param.h_sysDo -> SetLineColor(MyPalette[i]);

        //      param.h_sysDo -> SetLineStyle(2);

        double maxUp = histo -> GetMaximum();
        //      double maxDo = param.h_sysDo -> GetMaximum();
        double minUp = histo -> GetMinimum();
        //      double minDo = param.h_sysDo -> GetMinimum();
        if (max < maxUp) max = maxUp;
        if (min > minUp) min = minUp;
        //      if (max < maxDo) max = maxDo;
        //      if (min > minDo) min = minDo;

        latex -> DrawTextNDC(0.2, 0.85, param.h_nominal -> GetName());

        leg -> AddEntry(histo, (param.name + dirs[iD]).c_str(), "l");
      }

      float limitYrange = 20;
      if (limitYrange > 0 && max > limitYrange)
        max = limitYrange;
      if (limitYrange > 0 && -min > limitYrange)
        min = -limitYrange;
      double maxFactor = 1.1;
      double delta = max - min;
      double maxHist = min + maxFactor * delta;
      double minHist = max - maxFactor * delta;
      firstHisto -> GetYaxis() -> SetRangeUser(minHist, maxHist);
      firstHisto -> GetYaxis() -> SetTitle("Shift [%]");

      leg -> Draw();

      string plotFileName = m_plotDir + fileTag + "_" + region
              + "_compare_" + groups[iG] + dirs[iD] + m_graphicsExt;
      canvas -> Print(plotFileName.c_str());

      delete canvas;
      delete latex;
    }
  }
}

void plotCompareSys(vector<SysParam*> sysParamsTest, vector<SysParam*> sysParamsRef) {

  m_tsize = 0.07;
  setStyleAtlas();
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.10);

  TCanvas * canvas;

  for (unsigned int i = 0; i < sysParamsTest.size(); i++) {
    string name = sysParamsTest[i] -> name;
    int index = findSys(name, sysParamsRef);
    if (index >= 0) {
      cout<<"Running "<<name<<endl;
      SysParam test = *sysParamsTest[i];
      SysParam ref = *sysParamsRef[index];

      bool nom_only = 0;
      if (ref.name == "Nominal" && test.name == "Nominal") nom_only = 1;
      if (!nom_only) {
        canvas = new TCanvas("canvas", "canvas", 600, 800);
        canvas -> Divide(1, 3);
      }
      else {
        canvas = new TCanvas("canvas", "canvas", 600, 533);
        canvas -> Divide(1, 2);
      }

      test.h_nominal -> SetLineStyle(1);
      test.h_sysUp -> SetLineStyle(1);
      test.h_sysDo -> SetLineStyle(1);
      test.h_sysUpNorm -> SetLineStyle(1);
      test.h_sysDoNorm -> SetLineStyle(1);

      ref.h_nominal -> SetLineStyle(2);
      ref.h_sysUp -> SetLineStyle(2);
      ref.h_sysDo -> SetLineStyle(2);
      ref.h_sysUpNorm -> SetLineStyle(2);
      ref.h_sysDoNorm -> SetLineStyle(2);

      TLatex* latex = new TLatex();
      vector<TH1*> histos;

      canvas -> cd(1);
      histos.clear();
      histos.push_back(ref.h_nominal);
      histos.push_back(test.h_nominal);
      drawHistograms(histos, true);
      ref.h_nominal -> GetYaxis() -> SetTitle("Events");
      latex -> DrawTextNDC(0.3, 0.88, ref.h_nominal -> GetName());

      TLegend* leg = new TLegend(0.75, 0.8, 0.9, 0.95);
      leg -> SetFillStyle(0);
      leg -> SetBorderSize(0);
      leg -> SetTextSize(m_tsize);
      leg -> AddEntry(test.h_nominal, m_testName.c_str(), "l");
      leg -> AddEntry(ref.h_nominal, m_refName.c_str(), "l");
      leg -> Draw();

      canvas -> cd(2);
      histos.clear();
      if (!nom_only) {
        histos.push_back(ref.h_sysUp);
        histos.push_back(ref.h_sysDo);
        histos.push_back(test.h_sysUp);
        histos.push_back(test.h_sysDo);
        drawHistograms(histos, false, 20);
        ref.h_sysUp -> GetYaxis() -> SetTitle("Shift [%]");
        latex -> DrawTextNDC(0.3, 0.88, (ref.name + " (full)").c_str());

        leg = new TLegend(0.8, 0.8, 0.9, 0.95);
        leg -> SetFillStyle(0);
        leg -> SetBorderSize(0);
        leg -> SetTextSize(m_tsize);
        leg -> AddEntry(test.h_sysUp, "Up", "l");
        leg -> AddEntry(test.h_sysDo, "Down", "l");
        leg -> Draw();
      }
      else {
        TH1* h_ratio = (TH1F*)test.h_nominal->Clone();
        h_ratio->Divide(ref.h_nominal);
	//avoid that ratio goes to zero in case of empty bins
	//(and increases the y-axis range unnecessarily)
	setEmptyBinsToOne(h_ratio); 
        histos.push_back(h_ratio);
        drawHistograms(histos, false, 20);
        h_ratio -> GetYaxis() -> SetTitle(Form("%s / %s",m_testName.c_str(),m_refName.c_str()));
      }

      if (!nom_only) {
        canvas -> cd(3);
        histos.clear();
        histos.push_back(ref.h_sysUpNorm);
        histos.push_back(ref.h_sysDoNorm);
        histos.push_back(test.h_sysUpNorm);
        histos.push_back(test.h_sysDoNorm);
        drawHistograms(histos, false, 20);
        ref.h_sysUpNorm -> GetYaxis() -> SetTitle("Shift [%]");
        latex -> DrawTextNDC(0.3, 0.88, (ref.name + " (shape)").c_str());
        leg -> Draw();
      }

      string plotFileName = m_plotDir + ref.h_nominal -> GetName() + "_" + ref.name + m_graphicsExt;
      canvas -> Print(plotFileName.c_str());

      delete canvas;
    }
  }
}

void limitFileCheck() {

  // define files
  //string dir = "root://eosatlas//eos/atlas/atlascerngroupdisk/phys-higgs/HSG5/Run1Paper/";
  //string fileNameTest = dir + "Tsukuba/LimitHistograms.all.LP2013.8TeV.CUT.TSUKUBA.noSys.v0.1.root";
  //string fileNameRef = dir + "CPPM/LimitHistograms.vvbb.8TeV.CUT.CPPM.v0.0_StatOnlyDirectTag.root";
  //string fileNameTest = dir + "CERN/LimitHistograms.lvbb.LHCP2013.8TeV.CERN.v5.2.mJJ.root";
  //string fileNameRef = dir + "CERN/LimitHistograms.lvbb.LHCP2013.8TeV.CERN.v5.1.mJJ.root";
  //string fileNameTest = dir + "CPPM/LimitHistograms.lvbb.8TeV.MVA_WH125.CPPM.v4.3.root";
  //string fileNameRef = dir + "CPPM/LimitHistograms.lvbb.8TeV.MVA.CPPM.v3.1.root";
  //string fileNameRef = dir + "Iowa/LimitHistograms.lvbb.8TeV.MVA_WlvH125.Iowa.v8.5.4.root";
  //string fileNameRef = dir + "Iowa/LimitHistograms.lvbb.8TeV.MVA_WlvH125.Iowa.v8.3.root";

  //Run2
  //-----
  m_testName = "ttag";
  m_refName = "dtag";
  string dir1 = "/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/cxrdr_hists/"+date+"/"+m_refName+"/";
  string dir2 = "/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/cxrdr_hists/"+date+"/"+m_testName+"/";
  //  string fileNameRef = dir1 + "LimitHistograms.AZh.vvbb.Morinond2016.full.13TeV.CUT.LiverpoolBmham.v04c.root"; 
  string fileNameTest = dir2 + "inputsFile-all.root";//LimitHistograms.vvbb.13TeV.CUT_VHbb.AcSinica.v1.0.EXCL.root";
  string fileNameRef = dir1+"inputsFile-all.root";//LimitHistograms.AZh.vvbb.Morinond2016.full.13TeV.CUT.LiverpoolBmham.v05.root";
  if(m_nominalOnly){
    dir1 = "/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/cxrdr_hists/"+date+"/nom/"+m_refName+"/";
    dir2 = "/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/cxrdr_hists/"+date+"/nom/"+m_testName+"/";
    fileNameTest = dir2 + "inputsFile_nom.root";//LimitHistograms.vvbb.13TeV.CUT_VHbb.AcSinica.v1.0.EXCL.root";
    fileNameRef = dir1+"inputsFile_nom.root";//LimitHistograms.AZh.vvbb.Morinond2016.full.13TeV.CUT.LiverpoolBmham.v05.root";
  } 
  
  // define subfolders, set e.g. to "OneLepton" for the merged "FinalInputs"
  //string subFolderTest = "OneLepton";
  string subFolderTest = "";
  string subFolderRef = subFolderTest;

  // table headers for printout and legend labels for plotting
  // set to something meaningful, use less than m_width = 10 characters
  // m_testName = "CPPM_m43";
  // m_refName = "Iow_m854";
 


  //m_plotDir = "./plots/" + m_refName + "_" + m_testName + "_";
  if(!createOutputDir()) return;

  // samples for yield comparison (not all have to be in the limit file)
  vector<string> samples ={"Zbb","Zbc","Zbl","Zcc","Zcl","Zl","ZZ","qqZllH125","ttbar"};
  samples ={"Zbb","qqZvvH125"};
  // samples={"qqWlvH125","Wv22bb"};
  // "regions" to look at
  // vector<string> tags={"2tag"},jets={"2jet","3pjet"},ptvs={"0_120ptv","120ptv"},regs={"topemucr","SR"},vars={"mBB","pTB1","pTB2","dEtaBB","dRBB"},regions;//Run 1
  vector<string> tags={"2tag"},jets={"2jet","3jet","4pjet"},ptvs={"0_150ptv","150_200ptv","200ptv"},regs={/*"topemucr"*/"highmBBcr","lowmBBcr","SR"},vars={"mBB","dRBB"},regions;//mVH
  regs={"SR"};
  // if(!m_nominalOnly)vars={"mBB","mVH"};
  for(auto t:tags){
    for(auto j:jets){
      for(auto p:ptvs){
	for(auto r:regs){
	  for(auto v:vars)regions.push_back(t+j+"_"+p+"_"+r+"_"+v);
	}
      }
    }
  }
  string regionTest, regionRef;
//   regions.push_back(regionTest);


   for (unsigned int j = 0; j < regions.size(); j++) {
     regionTest = regions[j];
     regionRef = regions[j];

     for (unsigned int i = 0; i < samples.size(); i++) {
           string sample = samples[i];
	   // if(sample=="ttbar" && regions[j].find("topemu")==string::npos)continue; 

  //   m_plotDir = "./plots/full8v28tt_" + subFolderRef + "_" + subFolderTest + "_";

  cout << "================================================================================" << endl;
  cout << "Reading yields from test file: " << fileNameTest << endl;
  vector<Yield> yieldsTest = readYields(fileNameTest, samples, regionTest, subFolderTest);
  cout << "Reading yields from reference file: " << fileNameRef << endl;
  vector<Yield> yieldsRef = readYields(fileNameRef, samples, regionRef, subFolderRef);
  cout << "================================================================================" << endl;
  cout << "Reading systematics on " << sample << " in " << regionTest << " from test file: " << fileNameTest << endl;
  vector<SysParam*> sysParamsTest = readSysFromLimitFile(fileNameTest, sample, regionTest, subFolderTest);
  cout << "Reading systematics on " << sample << " in " << regionRef << " from reference file: " << fileNameRef << endl;
  vector<SysParam*> sysParamsRef = readSysFromLimitFile(fileNameRef, sample, regionRef, subFolderRef);
  //vector<SysParam> sysParamsRef = readSysFromLogFile(fileNameRef, sample, region);
  if(!m_quiet){
    cout << "================================================================================" << endl;
    cout << "Printing systematics on " << sample << " in " << regionTest << " from test file: " << fileNameTest << endl;
    printSysParam(sysParamsTest);
    cout << "================================================================================" << endl;
    cout << "Printing systematics on " << sample << " in " << regionRef << " from reference file: " << fileNameRef << endl;
    printSysParam(sysParamsRef);
    cout << "================================================================================" << endl;
    cout << "Additional systematics on " << sample << " in " << regionRef << " in reference file: " << fileNameRef << endl;
    checkForMissingSys(sysParamsTest, sysParamsRef);
    cout << "================================================================================" << endl;
    cout << "Additional systematics on " << sample << " in " << regionTest << " in test file: " << fileNameTest << endl;
    checkForMissingSys(sysParamsRef, sysParamsTest);
    cout << "================================================================================" << endl;
    cout << "Comparing systematics on " << sample << " in " << regionTest << endl;
    compareSys(sysParamsTest, sysParamsRef);
    cout << "================================================================================" << endl;
    cout << "Yield comparison in " << regionTest << ": " << endl;
    compareYields(yieldsTest, yieldsRef);
    cout << "================================================================================" << endl;
  }
  if (m_drawSingle) {
    cout << "Drawing systematics from test file..." << endl;
    plotSys(sysParamsTest, m_testName);
    cout << "Drawing systematics from reference file..." << endl;
    plotSys(sysParamsRef, m_refName);
  }
  if (m_drawGroup) {
    cout << "Drawing systematics groups from test file..." << endl;
    plotGroupSys(sysParamsTest, m_testName);
  }
  if (m_drawCompare) {
    cout << "Drawing systematics comparing both files..." << endl;
    plotCompareSys(sysParamsTest, sysParamsRef);
  }
  if (m_drawSingle || m_drawCompare) {
    cout << "================================================================================" << endl;
  }
  
  if (m_drawSingle || m_drawGroup || m_drawCompare) {
    for (unsigned int i = 0; i < sysParamsTest.size(); i++) {
      SysParam* param = sysParamsTest.at(i);
      if (i==0) {
        delete param->h_nominal;
      }
      delete param->h_sysUp;
      delete param->h_sysDo;
      delete param->h_sysUpNorm;
      delete param->h_sysDoNorm;
      delete param;
    }
    for (unsigned int i = 0; i < sysParamsRef.size(); i++) {
      SysParam* param = sysParamsRef.at(i);
      if (i==0) {
        delete param->h_nominal;
      }
      delete param->h_sysUp;
      delete param->h_sysDo;
      delete param->h_sysUpNorm;
      delete param->h_sysDoNorm;
      delete param;
    }
  }

  
     }
  //  }
   }
  gSystem->Exit(0);
}
