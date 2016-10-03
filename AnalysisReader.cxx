#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "EventLoop/OutputStream.h"
#include <CxAODReader/AnalysisReader.h>

#include "xAODEventInfo/EventInfo.h"
#include "xAODBTagging/BTagging.h"

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "CxAODTools/DummyEvtSelection.h"

// this is needed to distribute the algorithm to the workers
ClassImp(AnalysisReader)

AnalysisReader::AnalysisReader () :
  m_eventInfoReader(nullptr),
  m_METReader(nullptr),
  m_MPTReader(nullptr),
  m_electronReader(nullptr),
  m_photonReader(nullptr),
  m_muonReader(nullptr),
  m_tauReader(nullptr),
  m_jetReader(nullptr),
  m_fatJetReader(nullptr),
  m_trackJetReader(nullptr),
  m_truthParticleReader(nullptr),
  m_truthMuonReader(nullptr),
  m_truthElectronReader(nullptr),
  m_truthNeutrinoReader(nullptr),
  m_truthWZJetReader(nullptr),
  m_mmcReader(nullptr),
  m_eventInfo(nullptr),
  m_mmc(nullptr),
  m_electrons(nullptr),
  m_photons(nullptr),
  m_muons(nullptr),
  m_taus(nullptr),
  m_jets(nullptr),
  m_fatJets(nullptr),
  m_trackJets(nullptr),
  m_metCont(nullptr),
  m_met(nullptr),
  m_met_soft(nullptr),
  m_mptCont(nullptr),
  m_mpt(nullptr),
  m_truthParts(nullptr),
  m_truthMuons(nullptr),
  m_truthElectrons(nullptr),
  m_truthNeutrinos(nullptr),
  m_truthWZJets(nullptr),
  m_bTagTool(nullptr),
  m_eventCounter(0),
  m_eventCountPassed(0),
  m_isMC(false),
  m_mcChannel(-999),
  m_weight(1.),
  m_jetCorrectionOld(true),
  m_maxEvents(-1),
  m_isSherpaVJets(false),
  m_isSherpaPt0WJets(false),
  m_isSherpaPt0ZJets(false),
  m_randomRunNumber(-1),
  m_config(nullptr),
  m_debug(false),
  m_validation(false),
  m_passThroughOR(false),
  m_applyVarOR(false),
  m_applyEventPreSelection(false),
  m_applySherpaTruthPtCut(false),
  m_applyPowhegTruthMttCut(false),
  m_usePowhegInclFraction(-1),
  m_recomputePUWeight(false),
  m_checkEventDuplicates(false),
  m_failEventDuplicates(false),
  m_putAllSysInOneDir(false),
  m_recomputeMuTrigSF(false),
  m_readOROutputLabels(false),
  m_useOverlapRegister(false),
  m_sumOfWeightsFile(""),
  m_histSvc(nullptr),
  m_histNameSvc(nullptr),
  m_overlapRemoval(nullptr),
  m_eventSelection(nullptr),
  m_eventPostSelection(nullptr),
  m_xSectionProvider(nullptr),
  m_sumOfWeightsProvider(nullptr),
  m_triggerTool(nullptr),
  m_puReweightingTool(nullptr),
  m_corrsAndSysts(nullptr),
  m_overlapRegAcc(nullptr),
  m_myNNLO(nullptr),
  m_doQCD(false)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}

EL::StatusCode AnalysisReader::setupJob (EL::Job &job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  if (m_debug) Info("setupJob()", "Setting up job.");

  job.useXAOD();

  // check if ConfigStore is set
  if (!m_config) {
    Error("setupJob()", "ConfigStore not setup! Remember to set it via setConfig()!");
    return EL::StatusCode::FAILURE;
  }

  bool writeMVATree = false,sysonly=false;

  EL::OutputStream out("MVATree");
  m_config->getif<bool>("writeMVATree", writeMVATree);
  m_config->getif<bool>("systematicsOnly", m_sysonly);
  Info("setupJob()", "sysonly is %i",m_sysonly);
  writeMVATree&=!m_sysonly;

  if (writeMVATree) job.outputAdd(out);

  bool writeOSTree = false;
  EL::OutputStream outOS("OSTree");
  m_config->getif<bool>("writeOSTree", writeOSTree);
  if (writeOSTree) job.outputAdd(outOS);
 
  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init("MyxAODAnalysis").ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
} // setupJob

EL::StatusCode AnalysisReader::histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  Info("histInitialize()", "Initializing histograms.");

  // histogram manager
  m_histSvc     = new HistSvc();
  m_histNameSvc = new HistNameSvc();
  m_histSvc->SetNameSvc(m_histNameSvc);
  m_histSvc->SetWeightSysts(&m_weightSysts);

  bool fillHists = true;
  m_config->getif<bool>("writeHistograms", fillHists);
  m_histSvc->SetFillHists(fillHists);

  return EL::StatusCode::SUCCESS;
} // histInitialize

EL::StatusCode AnalysisReader::initializeNNLORW (int DSID)
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader::fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  // std::cout << "fileExecute Input file is " << wk()->inputFile()->GetName() << std::endl;


  return EL::StatusCode::SUCCESS;
} // fileExecute

EL::StatusCode AnalysisReader::changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  // It seems PROOF-lite needs the file initialization here

  TFile  *inputfile = wk()->inputFile();
  TString filename(inputfile->GetName());

  Info("changeInput()", "Processing file '%s'", filename.Data());

  m_isSherpaPt0WJets = ((filename.Contains("Sherpa_CT10_W") || filename.Contains("Sh_CT10_W")) && (filename.Contains("Pt0") && !filename.Contains("Pt70")));
  m_isSherpaPt0ZJets = ((filename.Contains("Sherpa_CT10_Z") || filename.Contains("Sh_CT10_Z")) && (filename.Contains("Pt0") && !filename.Contains("Pt70")));

  // general Sherpa flag
  m_isSherpaVJets = ((filename.Contains("Sherpa") || filename.Contains("Sh_CT10")) && filename.Contains("Pt"));

  return EL::StatusCode::SUCCESS;
} // changeInput

EL::StatusCode AnalysisReader::initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  m_config->getif<bool>("systematicsOnly", m_sysonly);

  EL_CHECK("initialize()",initializeEvent());
  EL_CHECK("initialize()",initializeReader());
  EL_CHECK("initialize()",initializeVariations());
  EL_CHECK("initialize()",initializeSelection());
  EL_CHECK("initialize()",initializeTools());
  EL_CHECK("initialize()",initializeIsMC());
  EL_CHECK("initialize()",initializeSumOfWeights());
  EL_CHECK("initialize()",initializeCorrsAndSysts());

  m_config->getif<bool>("validation", m_validation);
  m_config->getif<int>("maxEvents", m_maxEvents);
  m_config->getif<bool>("jetCorrectionOld", m_jetCorrectionOld);

  if (m_validation) {
    EL_CHECK("initialize()", initializeValidationSelection());
  }

  EL_CHECK("changeInput()", initializeChannel());
  if(m_mcChannel==410000){ // Initalize the NNLO reweighting just for nominal ttbar
    EL_CHECK("initialize()",initializeNNLORW(m_mcChannel));
  }

  return EL::StatusCode::SUCCESS;
} // initialize

EL::StatusCode AnalysisReader::initializeEvent ()
{
  Info("initializeEvent()", "Initialize event.");

  m_event = wk()->xaodEvent();

  Info("initializeEvent()", "Number of events on first file = %lli", m_event->getEntries()); // print long long int

  m_eventCounter = 0;

  m_config->getif<bool>("debug", m_debug);

  // luminosity (for rescaling of MC to desired lumi, in fb-1)
  // default is 0.001 fb-1, which means no rescaling of the MC inputs (already scaled to 1 pb-1)
  m_luminosity = 0.001;
  m_config->getif<float>("luminosity", m_luminosity);
  Info("initializeEvent()", "Luminosity for normalisation of MC = %f fb-1", m_luminosity);

  m_config->getif<bool>("passThroughOR", m_passThroughOR);
  //m_config->getif<bool>("applyVarOR", m_applyVarOR); //New Variable Muon-trackJet OverlapRemoval

  if (m_passThroughOR) Warning("AnalysisReader :: AnalysisReader ()", "Set OR to pass-through mode!");

  m_config->getif<bool>("applyEventPreSelection", m_applyEventPreSelection);
  m_config->getif<bool>("applySherpaTruthPtCut", m_applySherpaTruthPtCut);
  m_config->getif<bool>("applyPowhegTruthMttCut", m_applyPowhegTruthMttCut);
  m_config->getif<float>("usePowhegInclFraction", m_usePowhegInclFraction);
  m_config->getif<bool>("recomputePUWeight", m_recomputePUWeight);
  m_config->getif<bool>("checkEventDuplicates", m_checkEventDuplicates);
  m_config->getif<bool>("failEventDuplicates", m_failEventDuplicates);
  m_config->getif<bool>("putAllSysInOneDir", m_putAllSysInOneDir);
  m_config->getif<bool>("compute_muTrigSF", m_recomputeMuTrigSF);

  std::string whichData = "combined";
  m_config->getif<std::string>("whichData", whichData);
  if(whichData!="combined" && !m_recomputePUWeight){
    Error("PUReweightingTool::initialize()","Running on %s data only -> need to recompute PU weight. Set recomputePUWeight to true in config! Exiting.",whichData.c_str());
    return EL::StatusCode::FAILURE;
  }

  if(m_recomputeMuTrigSF && !m_recomputePUWeight){
    Error("PUReweightingTool::initialize()","Trying to recompute muon trigger SF -> need to recompute PU weight. Set recomputePUWeight to true in config! Exiting.");
    return EL::StatusCode::FAILURE;
  } 
  
  return EL::StatusCode::SUCCESS;
} // initializeEvent

EL::StatusCode AnalysisReader::initializeReader ()
{
  Info("initializeReader()", "Initialize object reader.");

  // note: the names of the readers determine which collections are
  //       read, see also framework-read.cfg

  std::string dummy="eventInfo";m_config->getif<std::string>("eventInfoContainer",dummy);    
  m_eventInfoReader     = registerReader<xAOD::EventInfo>(dummy);
  dummy="MET";                  m_config->getif<std::string>("METContainer",dummy);          
  m_METReader           = registerReader<xAOD::MissingETContainer>("MET");
  dummy="MPT";                  m_config->getif<std::string>("MPTContainer",dummy);          
  m_MPTReader           = registerReader<xAOD::MissingETContainer>("MPT");
  dummy="electron";             m_config->getif<std::string>("electronContainer",dummy);     
  m_electronReader      = registerReader<xAOD::ElectronContainer>(dummy);
  dummy="photon";               m_config->getif<std::string>("photonContainer",dummy);       
  m_photonReader        = registerReader<xAOD::PhotonContainer>(dummy);
  dummy="muon";                 m_config->getif<std::string>("muonContainer",dummy);         
  m_muonReader          = registerReader<xAOD::MuonContainer>(dummy);
  dummy="tau";                  m_config->getif<std::string>("tauContainer",dummy);          
  m_tauReader           = registerReader<xAOD::TauJetContainer>(dummy);
  dummy="jet";                  m_config->getif<std::string>("jetContainer",dummy);          
  m_jetReader           = registerReader<xAOD::JetContainer>("jet");
  dummy="fatJet";               m_config->getif<std::string>("fatJetContainer",dummy);       
  m_fatJetReader        = registerReader<xAOD::JetContainer>("fatJet");
  dummy="trackJet";             m_config->getif<std::string>("trackJetContainer",dummy);     
  m_trackJetReader      = registerReader<xAOD::JetContainer>("trackJet");
  dummy="truthParticle";         m_config->getif<std::string>("truthParticleContainer",dummy);
  m_truthParticleReader = registerReader<xAOD::TruthParticleContainer>("truthParticle");
  m_truthMuonReader     = registerReader<xAOD::TruthParticleContainer>("truthMuon");
  m_truthElectronReader = registerReader<xAOD::TruthParticleContainer>("truthElectron");
  m_truthNeutrinoReader = registerReader<xAOD::TruthParticleContainer>("truthNeutrino");
  m_truthWZJetReader    = registerReader<xAOD::JetContainer>("truthWZJet");
  m_mmcReader           = registerReader<xAOD::EventInfo>("ditau");
  return EL::StatusCode::SUCCESS;
} // initializeReader

EL::StatusCode AnalysisReader::initializeSelection ()
{
  Info("initializeSelection()", "Initialize selection.");

  // set event selection: use a dummy here
  m_eventSelection = new DummyEvtSelection();
  m_eventPostSelection = new DummyEvtSelection();
  m_fillFunction   = std::bind(&AnalysisReader::fill_example, this);

  return EL::StatusCode::SUCCESS;
} // initializeSelection

EL::StatusCode AnalysisReader::initializeCorrsAndSysts ()
{
  if(!m_isMC) return EL::StatusCode::SUCCESS;
  Info("initializeCorrsAndSysts()", "Initialize CorrsAndSysts.");
  m_corrsAndSysts = new CorrsAndSysts("8TeV_ZeroLepton"); // dummy : re-initialized in AnalysisReader_VHbb                                                                                                           
  return EL::StatusCode::SUCCESS;
} // initializeCorrsAndSysts                                                                                                                                                                                         

EL::StatusCode AnalysisReader::initializeTools ()
{
  Info("initializeTools()", "Initialize tools.");

  // overlap removal
  // ---------------
  m_overlapRegAcc = nullptr;
  m_metMJCalc=0;
  m_config->getif<int>("metMJCalc",m_metMJCalc);
  m_jetPtCut=25.0e3;
  m_config->getif<float>("jetPtCut",m_jetPtCut);
  m_readOROutputLabels=false;
  m_config->getif<bool>("readOROutputLabels",m_readOROutputLabels);
  m_useOverlapRegister = false;
  m_config->getif<bool>("useOverlapRegister",m_useOverlapRegister);
  if ( m_useOverlapRegister ) {
    // instantiate overlapRegisterAccessor in read mode
    m_overlapRegAcc = new OverlapRegisterAccessor(OverlapRegisterAccessor::READ);
  }
  else {
    // soon obsolete!!
    m_overlapRemoval = new OverlapRemoval(*m_config);
    EL_CHECK("initializeTools()", m_overlapRemoval->initialize());
  }

  // b-tagging
  // ---------
  bool use2DbTagCut = false;
  m_config->getif<bool>("use2DbTagCut", use2DbTagCut);
  
  std::vector<std::string> bTagToolConfigs;
  if (use2DbTagCut) m_config->getif<std::vector<std::string> >("bTagToolConfigs2D", bTagToolConfigs);
  else m_config->getif<std::vector<std::string> >("bTagToolConfigs", bTagToolConfigs);
  std::vector<std::string> bTagVariationConfigs;
  m_config->getif<std::vector<std::string> >("weightVariations", bTagVariationConfigs);
  std::string reductionScheme = "Medium";
  bool doWeightVar            = false;

  for (auto var : bTagVariationConfigs) {
    if (var.find("BTAG") != std::string::npos) {
      doWeightVar = true;

      if (var.find("LOOSE") != std::string::npos) reductionScheme = "Loose";

      if (var.find("MEDIUM") != std::string::npos) reductionScheme = "Medium";

      if (var.find("TIGHT") != std::string::npos) reductionScheme = "Tight";
    }
  }

  doWeightVar&=!m_sysonly;//sysonly means we're only doing systematics--so don't do btagging systematics, which are done on the nominal
  if (bTagToolConfigs.size() >= 4) {
    BTaggingTool::Config_t args {
      { "TaggerName", bTagToolConfigs[0] },
        { "OperatingPoint", bTagToolConfigs[1] },
          { "JetAuthor", bTagToolConfigs[2] },
            { "Scheme", bTagToolConfigs[3] },
              { "rScheme", reductionScheme }
    };

    if (!m_bTagTool) m_bTagTool = new BTaggingTool();
    EL_CHECK("initializeTools()", m_bTagTool->initialize(args, use2DbTagCut));
    m_bTagTool->setWeightVar(doWeightVar);
  } else {
    Warning("initializeTools()",
            "Could not initialize BTaggingTool due to invalid bTagToolConfigs in config!");
  }


  // pileup reweighting tool
  // ---------------------------
  m_puReweightingTool = new PUReweightingTool(*m_config);
  if(m_recomputePUWeight) EL_CHECK("EventInfoHandler::initialize()", m_puReweightingTool->initialize());

  return EL::StatusCode::SUCCESS;
} // initializeTools

EL::StatusCode AnalysisReader::initializeValidationSelection ()
{
  Info("initializeValidationSelection()", "Initialize validation selection");

  m_validationFillFunction = std::bind(&AnalysisReader::fill_validation, this);
  return EL::StatusCode::SUCCESS;
} // initializeValidationSelection

EL::StatusCode AnalysisReader::initializeIsMC ()
{
  Info("initializeIsMC()", "Initialize isMC.");

  // get nominal event info
  // -------------------------------------------------------------
  const xAOD::EventInfo *eventInfo = m_eventInfoReader->getObjects("Nominal");

  if (!eventInfo) return EL::StatusCode::FAILURE;

  // get MC flag - different info on data/MC files
  // -----------------------------------------------
  m_isMC = Props::isMC.get(eventInfo);
  Info("initializeIsMC()", "isMC = %i", m_isMC);

  // COM energy
  std::string comEnergy = m_config->get<std::string>("COMEnergy");

  std::string xSectionFile = gSystem->Getenv("ROOTCOREBIN");
  xSectionFile      += "/data/FrameworkSub/XSections_";
  xSectionFile      += comEnergy;
  xSectionFile      += ".txt";
  // if xSectionFile is given in config file, replace all we just did
  m_config->getif<string>("xSectionFile", xSectionFile);
  m_xSectionProvider = new XSectionProvider(xSectionFile);

  if (!m_xSectionProvider) {
    Error("initializeXSections()", "XSection provider not initialized!");
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
} // initializeIsMC

EL::StatusCode AnalysisReader::initializeChannel (){
  //set MC channel number
  //----------------------
  const xAOD::EventInfo *eventInfo = m_eventInfoReader->getObjects("Nominal");

  if (!eventInfo) return EL::StatusCode::FAILURE;
  
  if (!m_xSectionProvider) {
    Error("initializeChannel()", "XSection provider not initialized!");
    return EL::StatusCode::FAILURE;
  }

  if(m_isMC) m_mcChannel = (int)eventInfo->mcChannelNumber(); 
  if(!m_isMC) m_histNameSvc->set_sample("data");
  else m_histNameSvc->set_sample(m_xSectionProvider->getSampleName(m_mcChannel));

  Info("initializeChannel()", Form("Initialize channel %d",m_mcChannel));

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader::initializeSumOfWeights ()
{
  Info("initializeSumOfWeights()", "Initialize sum of weights.");

  if (!m_isMC) {
    return EL::StatusCode::SUCCESS;
  }

  std::string sumOfWeightsFile;
  bool generateYieldFile = false;
  m_config->getif<bool>("generateYieldFile", generateYieldFile);
  bool mc15b = false;
  m_config->getif<bool>("mc15b", mc15b);
  if (generateYieldFile && (m_sumOfWeightsFile != "")) {
    sumOfWeightsFile = m_sumOfWeightsFile;
  } else {
    std::string comEnergy    = m_config->get<std::string>("COMEnergy");
    std::string analysisType = m_config->get<std::string>("analysisType");
    sumOfWeightsFile  = gSystem->Getenv("ROOTCOREBIN");
    sumOfWeightsFile += "/data/FrameworkSub/yields." + analysisType + ".";
    sumOfWeightsFile += comEnergy;
    if (mc15b) sumOfWeightsFile += "_mc15b";
    sumOfWeightsFile  += ".txt";
  }

  m_sumOfWeightsProvider = new sumOfWeightsProvider(sumOfWeightsFile);

  if (!m_sumOfWeightsProvider) {
    Error("initializeSumOfWeights()", "SumOfWeights not initialized!");
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
} // initializeSumOfWeights

EL::StatusCode AnalysisReader::initializeVariations ()
{
  Info("initializeVariations()", "Initialize systematic variations.");

  m_variations.clear();

  // retrieve variations from config if set
  std::vector<std::string> variations; //exp syst stored in CxAODs
  std::vector<std::string> csVariations; //Corrs & Syst modelling syst (weights!)
  std::vector<std::string> weightVariations; //b-tagging syst (weights!)
  //weight variations are dealt with elsewhere!

  bool nominalOnly            = false;
  bool autoDiscoverVariations = true;
  m_config->getif<std::vector<std::string> >("variations", variations);
  m_config->getif<std::vector<std::string> >("csVariations", csVariations);
  m_config->getif<std::vector<std::string> >("weightVariations", weightVariations);
  m_config->getif<bool>("nominalOnly", nominalOnly);
  m_config->getif<bool>("autoDiscoverVariations", autoDiscoverVariations);
  m_config->getif<bool>("doQCD", m_doQCD);
  if(!nominalOnly && m_doQCD) m_config->getif<std::vector<std::string> >("FakeFactorSyst", variations);
  // read available variations in input file
  TTree *collectionTree = dynamic_cast<TTree*>(wk()->inputFile()->Get("CollectionTree"));

  for (ObjectReaderBase *reader : m_objectReader) {
    reader->discoverVariations(collectionTree);
  }


  //need to specify whether to run nominalOnly or also syst variations
  //abort otherwise
  //--------------------
  if ((nominalOnly && (autoDiscoverVariations || variations.size() || csVariations.size() || weightVariations.size() ))
      || (!nominalOnly && !(autoDiscoverVariations || variations.size() || csVariations.size() || weightVariations.size())) ){
    Error("initializeVariations()", "Need to specify in config whether to run nominal only or systematic variations, too. Exiting.");
    return EL::StatusCode::FAILURE;
  }

  // nominal only
  // ---------------
  if (nominalOnly) {
    m_variations.push_back("Nominal");
    Info("initializeVariations()", "Running Nominal only!");
  }
  // if do syst variations
  // ----------------------
  else if (autoDiscoverVariations || variations.size()) {
    // create list - without double counting syst variations
    std::vector<std::string> availableVariations;
    if(!m_sysonly)availableVariations.push_back("Nominal");//for specialized running 
    if(!nominalOnly && m_doQCD){
      for (const std::string &varName : variations) {
     	availableVariations.push_back(varName);
      }
    }
    
    for (ObjectReaderBase *reader : m_objectReader) {
      if (m_debug) Info("initializeVariations()", "Container %s with %i variations:", reader->getContainerName().c_str(), (int)reader->getVariations().size());

      for (const std::string &varName : reader->getVariations()) {
        if (m_debug) Info("initializeVariations()", varName.c_str());

        // ensure to use each systematic once
        // would like to use set, but need to run Nominal first for OR (bug in EDM?)
        bool found = false;

        for (const std::string &varName1 : availableVariations) {
          if (varName1 == varName) {
            found = true;
            break;
          }
        }

        if (!found) availableVariations.push_back(varName);
        else if (m_debug) Info("initializeVariations()", "Already in list!");
      }
    }

    // print variations found in input file to screen
    // -----------------------------------------------
    Info("initializeVariations()", "Variations found in input file:");

    for (const std::string &varName1 : availableVariations) {
      Info("initializeVariations()", varName1.c_str());
    }

    // autoDiscoverVariations
    // ----------------------
    if (autoDiscoverVariations) {
      m_variations = availableVariations;
      Info("initializeVariations()", "Running all variations found in input file.");
    }

    // get variations from config file
    // --------------------------------
    else {
      // make sure nominal is run in any case
      if(!m_sysonly)variations.push_back("Nominal");

      Info("initializeVariations()", "Running only variations specified in config file (+ Nominal):");

      for (std::string varName : availableVariations) {
        // check if specified in config
        for (std::string systName : variations) {
          if (varName.find(systName) != std::string::npos) {
            m_variations.push_back(varName.c_str());
            break;
          }
        }
      }

      if (m_debug) {
        // check if all variations specified in config are available in input
        // and if not, print which
        for (std::string systName : variations) {
          for (std::string varName : m_variations) {
            if (varName.find(systName) == std::string::npos) {
              Warning("initializeVariations()", "Systematic variations %s not found in input file. Skipping it.", systName.c_str());
            }
          }
        }
      }

      // print variations to screen
      // ----------------------------
      for (std::string varName : m_variations) {
        std::cout << varName << std::endl;
      }
    }
  }
  //just doing weight variations
  else{
    m_variations.push_back("Nominal");
    Info("initializeVariations()", "Running only weight variations specified in config file (+ Nominal).");
  }
 
  return EL::StatusCode::SUCCESS;
} // initializeVariations

EL::StatusCode AnalysisReader::finalizeTools () {
  Info("finalizeTools()", "Finalizing tools.");

  // overlap removal
  // ---------------
  // EL_CHECK("AnalysisReader::finalizeTools()", m_overlapRemoval->release());
  if (m_overlapRemoval) {
    delete m_overlapRemoval;
    m_overlapRemoval = nullptr;
  }

  // b-tagging
  // ---------
  // EL_CHECK("AnalysisReader::finalizeTools()", m_bTagTool->release());
  if (m_bTagTool) {
    delete m_bTagTool;
    m_bTagTool = nullptr;
  }

  // pileup reweighting
  // --------------------
  if (m_puReweightingTool) {
    delete m_puReweightingTool;
    m_puReweightingTool = nullptr;
  }


  // trigger tool
  // -------------
  if (m_triggerTool) {
    delete m_triggerTool;
    m_triggerTool = nullptr;
  }

  // overlapRegisterAccessor
  if (m_overlapRegAcc) {
    delete m_overlapRegAcc;
    m_overlapRegAcc = nullptr;
  }

  return EL::StatusCode::SUCCESS;
} // finalizeTools

EL::StatusCode AnalysisReader::setObjectsForOR (const xAOD::ElectronContainer*,
                                                const xAOD::PhotonContainer*,
                                                const xAOD::MuonContainer*,
                                                const xAOD::TauJetContainer*,
                                                const xAOD::JetContainer*,
                                                const xAOD::JetContainer*) {
  Error("AnalysisReader::setObjectsForOR",
        "No generic implementation available, please see https://its.cern.ch/jira/browse/CXAOD-176");
  return EL::StatusCode::FAILURE;
} // setObjectsForOR

EL::StatusCode AnalysisReader::copyConeTruthLabels (const xAOD::JetContainer *jets) {
  for (const xAOD::Jet *jet : *jets) {
    if (!Props::HadronConeExclTruthLabelID.exists(jet)) {
      Error("AnalysisReader::copyConeTruthLabels", "HadronConeExclTruthLabelID does not exist in input! Exiting.");
      return EL::StatusCode::FAILURE;
    }
    int flav = Props::HadronConeExclTruthLabelID.get(jet);
    Props::TruthLabelID.set(jet, flav);
  }
  return EL::StatusCode::SUCCESS;
} // copyConeTruthLabels

PROPERTY(PropsEmptyContFix, float, pt)
PROPERTY(PropsEmptyContFix, float, px)
EL::StatusCode AnalysisReader::clearEmptyContainersWithNonZeroSize () {
  // empty containers for replacement
  static const xAOD::ElectronContainer empty_electrons;
  static const xAOD::PhotonContainer   empty_photons;
  static const xAOD::MuonContainer     empty_muons;
  static const xAOD::TauJetContainer   empty_taus;
  static const xAOD::JetContainer empty_jets;
  static const xAOD::TruthParticleContainer empty_truthParts;

  // replace invalid particle containers with empty ones
  if (m_electrons) if (m_electrons->size()) if (!PropsEmptyContFix::pt.exists(m_electrons->at(0))) m_electrons = &empty_electrons;

  if (m_photons) if (m_photons->size()) if (!PropsEmptyContFix::pt.exists(m_photons->at(0))) m_photons = &empty_photons;

  if (m_muons) if (m_muons->size()) if (!PropsEmptyContFix::pt.exists(m_muons->at(0))) m_muons = &empty_muons;

  if (m_taus) if (m_taus->size()) if (!PropsEmptyContFix::pt.exists(m_taus->at(0))) m_taus = &empty_taus;

  if (m_jets) if (m_jets->size()) if (!PropsEmptyContFix::pt.exists(m_jets->at(0))) m_jets = &empty_jets;

  if (m_fatJets) if (m_fatJets->size()) if (!PropsEmptyContFix::pt.exists(m_fatJets->at(0))) m_fatJets = &empty_jets;

  if (m_trackJets) if (m_trackJets->size()) if (!PropsEmptyContFix::pt.exists(m_trackJets->at(0))) m_trackJets = &empty_jets;

  if (m_truthParts) if (m_truthParts->size()) if (!PropsEmptyContFix::px.exists(m_truthParts->at(0))) m_truthParts = &empty_truthParts;

  return EL::StatusCode::SUCCESS;
} // clearEmptyContainersWithNonZeroSize

EL::StatusCode AnalysisReader::execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  if (((m_eventCounter % 10000) == 0) || m_debug) {
    Info("execute()", "Event number = %li", m_eventCounter);
  }
  m_eventCounter++;

  // b-tagging requires special index for MC-to-MC efficiencies
  // if (m_isMC) m_bTagTool->setMCIndex(m_bTagTool->indexMCEfficiencyFromFileName(wk()->inputFile()->GetName()));
  // if (m_isMC) m_bTagTool->setMCIndex(m_bTagTool->indexMCEfficiencyFromChannel(m_mcChannel, *m_xSectionProvider), *m_xSectionProvider);
  if (m_isMC) m_bTagTool->setMCIndex(m_mcChannel, *m_xSectionProvider);

  // store random number from pu tool -> working with nominal is enough
  if(m_recomputePUWeight){
    m_randomRunNumber = m_puReweightingTool->GetRandomRunNumber( m_eventInfoReader->getObjects("Nominal") );
  }

  // loop on systematic variations - or at least run "Nominal"!
  // ----------------------------------------------------------
  for (std::string varName : m_variations) {
    bool isNominal = (varName == "Nominal");
     
    // HistNameSvc settings
    // ---------------------
    m_histNameSvc -> reset();
    m_histNameSvc -> set_variation(varName);
    m_histNameSvc -> set_doOneSysDir(m_putAllSysInOneDir); //put all syst in one directory
    m_weightSysts.clear();  // prepare list of weight variations
    
    m_currentVar = varName;

    if (m_debug) Info("execute ()", "Running variation %s", varName.c_str());

    // retrieve containers
    // --------------------
    if (m_debug) Info("execute ()", "Retrieving input containers...");

    if (m_eventInfoReader) m_eventInfo = m_eventInfoReader->getObjects(varName);

    if (m_electronReader) m_electrons = m_electronReader->getObjects(varName);

    if (m_photonReader) m_photons = m_photonReader->getObjects(varName);

    if (m_muonReader) m_muons = m_muonReader->getObjects(varName);

    if (m_tauReader) m_taus = m_tauReader->getObjects(varName);

    if (m_jetReader) m_jets = m_jetReader->getObjects(varName);

    if (m_fatJetReader) m_fatJets = m_fatJetReader->getObjects(varName);

    if (m_trackJetReader) m_trackJets = m_trackJetReader->getObjects(varName);

    if (m_truthParticleReader) m_truthParts = m_truthParticleReader->getObjects(varName);
    if (m_truthMuonReader) m_truthMuons = m_truthMuonReader->getObjects(varName);
    if (m_truthElectronReader) m_truthElectrons = m_truthElectronReader->getObjects(varName);
    if (m_truthNeutrinoReader) m_truthNeutrinos = m_truthNeutrinoReader->getObjects(varName);
    if (m_truthWZJetReader) m_truthWZJets = m_truthWZJetReader->getObjects(varName);
    if (m_mmcReader) m_mmc = m_mmcReader->getObjects(varName);

    m_met = nullptr;
    m_met_soft = nullptr;

    if (m_METReader) {
      m_metCont = m_METReader->getObjects(varName);

      // if(m_doQCD&&m_currentVar=="MJ_El_METstr") m_metCont = m_METMJTightReader->getObjects(varName);
      // else m_metCont = m_METReader->getObjects(varName);

      if (m_metCont->size() == 1) {
        m_met = m_metCont->at(0);
      } else if (m_metCont->size() == 2) {
        m_met = m_metCont->at(0);
        m_met_soft = m_metCont->at(1);
      } else {
        Error("execute ()", "MET container '%s' has size != 1 or 2",
              m_METReader->getContainerName().c_str());
        return EL::StatusCode::FAILURE;
      }
    }
    m_mpt = nullptr;

    if (m_MPTReader) {
      m_mptCont = m_MPTReader->getObjects(varName);

      if (m_mptCont->size() > 0) {
        m_mpt = m_mptCont->at(0);
      } else {
        Error("execute ()", "MET container '%s' has size != 1",
              m_MPTReader->getContainerName().c_str());
        return EL::StatusCode::FAILURE;
      }
    }

    m_isMC = Props::isMC.get(m_eventInfo);
    if(m_isMC) m_mcChannel = (int)m_eventInfo->mcChannelNumber();
    if(!m_isMC) m_histNameSvc->set_sample("data");
    else m_histNameSvc->set_sample(m_xSectionProvider->getSampleName(m_mcChannel));
      
    // get event weight - fill m_weight, decorate eventInfo with individual weights
    // also get corrected mu in data
    // ---------------------------------
    if (m_debug) Info("execute ()", "Calculating event weight...");
    EL_CHECK("AnalysisReader::execute()", setEventWeight());

    // Event duplicates
    if (m_checkEventDuplicates) {
      EL_CHECK("execute()", checkEventDuplicates());
    }

    // Sherpa2.2 Njet weight (from Corrs&Sys)
    // ----------------------------------
    CAS::EventType type = CAS::NONAME;
    if ( ((m_mcChannel >= 363102) && (m_mcChannel <= 363122 ))            //Zvv
      || ((m_mcChannel >= 363361) && (m_mcChannel <= 363435 )) )          //Zll
       type = CAS::Z;
    else if ( ((m_mcChannel >= 363331) && (m_mcChannel <= 363354))        //Wtauv
           || ((m_mcChannel >= 363436) && (m_mcChannel <= 363483)) )      //Wlv
       type = CAS::W;
    float sherpa22weight = m_corrsAndSysts->Get_BkgNJetCorrection(type,(Props::NTruthWZJets20.exists(m_eventInfo) ? Props::NTruthWZJets20.get(m_eventInfo) : -1));
    m_weight *= sherpa22weight; 
    // -------------------------------------

    // SherpaTruthPt
    // --------------
    if ((varName == "Nominal") && m_isMC && m_isSherpaVJets && m_applySherpaTruthPtCut) {
      bool  passSherpaTruthPt;
      float sherpaTruthPt;
      EL_CHECK("AnalysisReader::execute()", checkSherpaTruthPt(passSherpaTruthPt, sherpaTruthPt));

      if (!passSherpaTruthPt) return EL::StatusCode::SUCCESS;

      m_histSvc->BookFillHist("VptTruth", 100, 0, 1000, sherpaTruthPt / 1000., m_weight);
    }

    // PowhegTruthMtt
    // --------------
    if (m_isMC && m_applyPowhegTruthMttCut) {
      bool  passPowhegTruthMtt;
      float powhegTruthMtt;
      EL_CHECK("AnalysisReader::execute()", checkPowhegTruthMtt(passPowhegTruthMtt, powhegTruthMtt));

      if (!passPowhegTruthMtt) return EL::StatusCode::SUCCESS;

      m_histSvc->BookFillHist("MttTruth", 100, 0, 5000, powhegTruthMtt / 1000., m_weight);
    }

    if (m_debug) Info("execute ()", "Event = %llu ; Run = %u", m_eventInfo->eventNumber(), m_eventInfo->runNumber());

    // fixes - TODO: check if still needed
    // -------
    // protection against empty containers with size > 0
    // TODO how can this happen?
    if (m_debug) Info("execute ()", "Clearing empty containers...");
    clearEmptyContainersWithNonZeroSize();

    // copy truth label flags (for BTagEfficiencyMapReader)
    if (m_isMC){
      if (m_debug) Info("execute ()", "Copying truth labels...");
      if (m_jets) EL_CHECK("AnalysisReader::execute()", copyConeTruthLabels(m_jets));
      if (m_trackJets) EL_CHECK("AnalysisReader::execute()", copyConeTruthLabels(m_trackJets));
    }

    // fill validation histograms
    // ----------------------------
    if (m_validation) {
      EL_CHECK("AnalysisReader::execute()", m_validationFillFunction());
    }

    // Overlap removal
    // -------------------
    // set the passPreSel flags (needed as input flag for the OR)
    if (m_debug) Info("execute ()", "Preparing objects for overlap removal...");
    EL_CHECK("AnalysisReader::execute()", setObjectsForOR(m_electrons, m_photons, m_muons, m_taus, m_jets, m_fatJets));

    if (m_passThroughOR) {
      // let all particles pass -> only valid if CxAOD was produced w/o kinematic variations
      OverlapPassThrough(m_electrons, m_photons, m_muons, m_taus, m_jets, m_fatJets);
    } else {
      // proper overlap removal
      if (m_useOverlapRegister) {
        // load overlap register container from TEvent
        EL_CHECK("AnalysisReader::execute()",m_overlapRegAcc->loadRegister(m_event));
        // decorate objects with OR decisions from register
	if(m_doQCD)varName = "Nominal";     
	OverlapRegisterAccessor::Containers containers;
        containers.variation = varName;
        containers.jets      = m_jets;
        containers.fatjets   = m_fatJets;
        containers.muons     = m_muons;
        containers.electrons = m_electrons;
        containers.taus      = m_taus;
        containers.photons   = m_photons;
        EL_CHECK("AnalysisReader::execute()",m_overlapRegAcc->decorateObjects( containers ));
      }else if(m_readOROutputLabels){ m_overlapRemoval->getOROutputLabels(m_electrons, m_photons, m_muons, m_taus, m_jets, m_fatJets); }
      else {
        // soon obsolete!!
        EL_CHECK("AnalysisReader::execute()", m_overlapRemoval->removeOverlap(m_electrons, m_photons, m_muons, m_taus, m_jets, m_fatJets));
      }
    
    }

    EL_CHECK("AnalysisReader::execute()", executePreEvtSel());

    // event selection
    // ------------------
    if (m_debug) Info("execute ()", "Checking the event selection...");
    bool passSel = true;

    if (m_eventSelection) {
      // all pointers are zero initialised in default constructor of SelectionContainers
      // so leaving out unused containers is okay (the corresponding pointer will not be dangling)
      SelectionContainers containers;
      containers.evtinfo   = m_eventInfo;
      containers.met       = m_met;
      containers.met_soft  = m_met_soft;
      containers.electrons = m_electrons;
      containers.photons   = m_photons;
      containers.muons     = m_muons;
      containers.taus      = m_taus;
      containers.jets      = m_jets;
      containers.fatjets   = m_fatJets;
      containers.trackjets = m_trackJets;
      containers.truthParticles = m_truthParts;
      containers.truthMuons     = m_truthMuons;
      containers.truthElectrons = m_truthElectrons;
      containers.truthNeutrinos = m_truthNeutrinos;
      containers.truthAntiKt4TruthWZJets = m_truthWZJets;

      if (m_applyEventPreSelection) {
	m_eventSelection->setJetPtCut(m_jetPtCut);
        passSel &= m_eventSelection->passPreSelection(containers, !isNominal);
        // Modify the met for the MJ
        m_eventSelection->setMetMJCalc(m_metMJCalc);
      }// end preselection
      passSel &= m_eventSelection->passSelection(containers, !isNominal);
    }

    if (!passSel) continue;

    if (isNominal) m_eventCountPassed++;

    // fill histograms
    // ----------------
    if (m_debug) Info("execute ()", "Calling the fill function...");
    EL_CHECK("AnalysisReader::execute()", m_fillFunction());
  }


  return EL::StatusCode::SUCCESS;
} // execute

EL::StatusCode AnalysisReader::setEventWeight () {
  // reset event weight
  m_weight = 1.0;

  // generator weight
  if (m_isMC) m_weight *= Props::MCEventWeight.get(m_eventInfo);

  // lumi weight
  EL_CHECK("AnalysisReader::setEventWeight()", applyLumiWeight());

  // PU weight (and get (corrected) mu in data...)
  EL_CHECK("AnalysisReader::setEventWeight()", applyPUWeight());
  return EL::StatusCode::SUCCESS;
} // setEventWeight

EL::StatusCode AnalysisReader::applyLumiWeight () {
  double weight = 1.;

  if (!m_isMC) return EL::StatusCode::SUCCESS;

  // return if weight not requested
  bool applyWeight = true;
  m_config->getif<bool>("applyLumiWeight", applyWeight);

  if (!applyWeight) {
    Props::LumiWeight.set(m_eventInfo, weight);
    return EL::StatusCode::SUCCESS;
  }

  // Query - reading the cross section every event.
  // This is because on the first event fileExecute and changeInput are called before initialize
  // so there is no event available to read dataset
  //m_mcChannel = (int)m_eventInfo->mcChannelNumber();
  // get sum of weights from text file
  if (!m_sumOfWeightsProvider) {
    Error("applyLumiWeight()", "SumOfWeights provider not initialized!");
    return EL::StatusCode::FAILURE;
  }
  double sumOfWeights = m_sumOfWeightsProvider->getsumOfWeights(m_mcChannel);

  //If we have a max number of events to run on we should scale the sample weight by that over the total available entries
  if(m_maxEvents > 0) {
    double availableEntries = m_sumOfWeightsProvider->getNEntriesSelectedOut(m_mcChannel);
    sumOfWeights *= m_maxEvents / availableEntries;
  }

  if (!m_xSectionProvider) {
    Error("applyLumiWeight()", "XSection provider not initialized!");
    return EL::StatusCode::FAILURE;
  }

  float sigmaEff = m_xSectionProvider->getXSection(m_mcChannel);

  if (m_debug) Info("applyLumiWeight()", "Cross section times eff. for dataset id %i = %f", m_mcChannel, sigmaEff);

  // we are normalising to MC lumi, sumOfWeights calculated per m_mcChannel
  weight = (sumOfWeights) ? sigmaEff / sumOfWeights : 1.0;

  // scale to desired luminosity
  weight *= (m_luminosity * 1e3); // the MC is scaled to 1pb-1 but m_luminosity is in fb-1

  // decorate eventInfo with the lumi weight
  Props::LumiWeight.set(m_eventInfo, weight);

  // multiply with global event weight
  m_weight *= weight;

  return EL::StatusCode::SUCCESS;
} // applyLumiWeight

EL::StatusCode AnalysisReader::applyPUWeight () {
  float weight = 1.;

  // take PU weight from CxAOD (if exists)
  if (!m_recomputePUWeight) {
    if (Props::Pileupweight.exists(m_eventInfo)) {
      weight = Props::Pileupweight.get(m_eventInfo);
      Props::PileupweightRecalc.set(m_eventInfo, weight);
      Props::averageInteractionsPerCrossingRecalc.set(m_eventInfo, Props::averageInteractionsPerCrossing.get(m_eventInfo));
    }
    else {
      Error("AnalysisReader::applyPUWeight ()", "Pileupweight not found in input!");
      return EL::StatusCode::FAILURE;
    }
  }
  else { // get PU weight from tool, decorate event info
    EL_CHECK("EventInfoHandler::executeEvent()", m_puReweightingTool->decorateWeight(m_eventInfo));
    weight = Props::PileupweightRecalc.get(m_eventInfo);
  }

  // multiply with global event weight
  bool applyWeight = true;
  m_config->getif<bool>("applyPUWeight", applyWeight);

  if (applyWeight) m_weight *= weight;

  return EL::StatusCode::SUCCESS;
} // applyPUWeight

EL::StatusCode AnalysisReader::checkEventDuplicates() {
  // test nominal only
  if (!m_histNameSvc -> get_isNominal()) {
    return EL::StatusCode::SUCCESS;
  }

  // count the event
  long int run = 0;
  if (m_isMC) run = m_eventInfo->mcChannelNumber();
  else        run = m_eventInfo->runNumber();
  long int evt = m_eventInfo->eventNumber();
  int count = 1;
  if (m_eventCountDuplicates.count(run)) {
    if (m_eventCountDuplicates[run].count(evt)) {
      count += m_eventCountDuplicates[run][evt];
    }
  }
  m_eventCountDuplicates[run][evt] = count;

  // check it
  if (count > 1) {
    Warning("checkEventDuplicates()", "Have %i counts for run %li event %li.", count, run, evt);
    if (m_failEventDuplicates) {
      Error("checkEventDuplicates()", "This is unacceptable! Exiting.");
      return EL::StatusCode::FAILURE;
    }
  }
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader::checkSherpaTruthPt (bool &pass, float &sherpaTruthPt) {
  // determine if it is a Sherpa Pt0 sample, and remove events overlapping with slice (flag set for 13TeV only)

  pass          = true;
  sherpaTruthPt = -999;

  if (!m_isSherpaVJets) {
    return EL::StatusCode::SUCCESS;
  }

  const xAOD::TruthParticleContainer *truthParts = m_truthParticleReader->getObjects("Nominal");

  if (!truthParts) {
    Error("checkSherpaTruthPt", "Did not find truth particles!");
    return EL::StatusCode::FAILURE;
  }

  const xAOD::TruthParticle *lep0 = nullptr;
  const xAOD::TruthParticle *lep1 = nullptr;

  for (const xAOD::TruthParticle *part : *truthParts) {
    // assume leptons are the only status 3 particles in SherpaV+jets
    if (part->status() == 3) {
      if (!lep0) lep0 = part;
      else lep1 = part;
    }
  }

  if (!(lep0 && lep1)) {
    Warning("checkSherpaTruthPt", "Did not find truth leptons!");
    return EL::StatusCode::SUCCESS;
  }

  TLorentzVector V = lep0->p4() + lep1->p4();
  sherpaTruthPt = V.Pt();

  if (m_isSherpaPt0WJets && (sherpaTruthPt > 40000.)) pass = false;

  if (m_isSherpaPt0ZJets && (sherpaTruthPt > 70000.)) pass = false;

  return EL::StatusCode::SUCCESS;
} // checkSherpaTruthPt

EL::StatusCode AnalysisReader::checkPowhegTruthMtt (bool &pass, float &powhegTruthMtt) {
  // determine if it is a Powheg m(ttbar) inclusive sample and
  // remove events overlapping with slices

  pass           = true;
  powhegTruthMtt = -1;

  if (!m_eventInfo) {
    Error("checkPowhegTruthMtt", "Did not find event info!");
    return EL::StatusCode::FAILURE;
  }

  int  channel                = m_eventInfo->mcChannelNumber();
  bool isPowhegTTbarInclusive = (channel == 410000);
  bool isPowhegTTbarMttSlice  = (channel >= 301528 && channel <= 301532);

  if (!(isPowhegTTbarInclusive || isPowhegTTbarMttSlice)) {
    return EL::StatusCode::SUCCESS;
  }

  if (!m_truthParts) {
    Error("checkPowhegTruthMtt", "Did not find truth particles!");
    return EL::StatusCode::FAILURE;
  }

  TLorentzVector ttbarVec(0., 0., 0., 0.);
  int nTopsFound = 0;

  for (auto *part : *m_truthParts) {
    // select status == 3 (before ISR/FSR) top or anti top
    if ((part->status() == 3) && ((part->pdgId() == 6) || (part->pdgId() == -6))) {
      ttbarVec += part->p4();
      nTopsFound++;
    }
  }

  if (nTopsFound != 2) {
    Error("checkPowhegTruthMtt", "Found wrong number of tops: %i", nTopsFound);
    return EL::StatusCode::FAILURE;
  }

  powhegTruthMtt = ttbarVec.M();

  if (isPowhegTTbarInclusive && (powhegTruthMtt > 1.1e6)) pass = false;

  if ((!pass || isPowhegTTbarMttSlice) &&
      (m_usePowhegInclFraction > 0) &&
      (m_usePowhegInclFraction < 1)) {
    float weight = m_usePowhegInclFraction;

    if (isPowhegTTbarMttSlice) weight = 1 - m_usePowhegInclFraction;
    m_weight *= weight;
    pass      = true;
  }

  return EL::StatusCode::SUCCESS;
} // checkPowhegTruthMtt

EL::StatusCode AnalysisReader::fill_example ()
{
  // get selection results
  ResultDummy selectionResult        = ((DummyEvtSelection*)m_eventSelection)->result();
  std::vector<const xAOD::Jet*> jets = selectionResult.jets;

  // there are also a global pointers to the current collections, e.g.
  // const xAOD::JetContainer*       m_jets;
  // but this one is not suspect to re-selection in the selection class

  // book and fill histograms
  m_histSvc->BookFillHist("jetN", 11, -0.5, 10.5, jets.size(), m_weight);

  for (const xAOD::Jet *jet : jets) {
    m_histSvc->BookFillHist("jetPt", 100, 0, 100, jet->pt() / 1e3, m_weight);
  }
  return EL::StatusCode::SUCCESS;
} // fill_example

EL::StatusCode AnalysisReader::fill_validation ()
{
  m_histSvc->BookFillHist("ELECTRON_nelectrons", 11, 0, 10, m_electrons->size(), m_weight);

  for (const xAOD::Electron *electron : *m_electrons) {
    m_histSvc->BookFillHist("ELECTRON_electronE", 100, 0, 200, electron->e() / 1e3, m_weight);
    m_histSvc->BookFillHist("ELECTRON_electronPt", 100, 0, 200, electron->pt() / 1e3, m_weight);
    m_histSvc->BookFillHist("ELECTRON_electronEta", 100, -5, 5, electron->eta(), m_weight);
    m_histSvc->BookFillHist("ELECTRON_electronPhi", 100, 0, 3.15, electron->phi(), m_weight);
  }

  m_histSvc->BookFillHist("MUON_nmuons", 11, 0, 10, m_muons->size(), m_weight);

  for (const xAOD::Muon *muon : *m_muons) {
    m_histSvc->BookFillHist("MUON_muonE", 100, 0, 200, muon->e() / 1e3, m_weight);
    m_histSvc->BookFillHist("MUON_muonPt", 100, 0, 200, muon->pt() / 1e3, m_weight);
    m_histSvc->BookFillHist("MUON_muonEta", 100, -5, 5, muon->eta(), m_weight);
    m_histSvc->BookFillHist("MUON_muonPhi", 100, 0, 3.15, muon->phi(), m_weight);
  }

  m_histSvc->BookFillHist("TAU_ntaus", 11, 0, 10, m_taus->size(), m_weight);

  for (const xAOD::TauJet *tau : *m_taus) {
    m_histSvc->BookFillHist("TAU_tauE", 100, 0, 200, tau->e() / 1e3, m_weight);
    m_histSvc->BookFillHist("TAU_tauPt", 100, 0, 200, tau->pt() / 1e3, m_weight);
    m_histSvc->BookFillHist("TAU_tauEta", 100, -5, 5, tau->eta(), m_weight);
    m_histSvc->BookFillHist("TAU_tauPhi", 100, 0, 3.15, tau->phi(), m_weight);
  }

  m_histSvc->BookFillHist("JET_njets", 11, 0, 10, m_jets->size(), m_weight);

  for (const xAOD::Jet *jet : *m_jets) {
    m_histSvc->BookFillHist("JET_jetE", 100, 0, 200, jet->e() / 1e3, m_weight);
    m_histSvc->BookFillHist("JET_jetPt", 100, 0, 200, jet->pt() / 1e3, m_weight);
    m_histSvc->BookFillHist("JET_jetEta", 100, -5, 5, jet->eta(), m_weight);
    m_histSvc->BookFillHist("JET_jetPhi", 100, 0, 3.15, jet->phi(), m_weight);
  }

  m_histSvc->BookFillHist("MET_met", 100, 0, 200, m_met->met() / 1e3, m_weight);
  m_histSvc->BookFillHist("MET_metPx", 100, 0, 200, m_met->mpx() / 1e3, m_weight);
  m_histSvc->BookFillHist("MET_metPy", 100, 0, 200, m_met->mpy() / 1e3, m_weight);
  if(m_met_soft!=nullptr){
    m_histSvc->BookFillHist("METSoft_met", 100, 0, 200, m_met_soft->met() / 1e3, m_weight);
    m_histSvc->BookFillHist("METSoft_metPx", 100, 0, 200, m_met_soft->mpx() / 1e3, m_weight);
    m_histSvc->BookFillHist("METSoft_metPy", 100, 0, 200, m_met_soft->mpy() / 1e3, m_weight);
  }
  return EL::StatusCode::SUCCESS;
} // fill_validation

EL::StatusCode AnalysisReader::postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  // clear event in object readers
  for (ObjectReaderBase *reader : m_objectReader) {
    reader->clearEvent();
  }

  return EL::StatusCode::SUCCESS;
} // postExecute

EL::StatusCode AnalysisReader::finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  Info("finalize()", "Finalizing job.");

  EL_CHECK("finalize", finalizeTools());
  
  if (m_eventSelection) {
    CutFlowCounter counter = m_eventSelection->getCutFlowCounter();
    TH1D* cutflowHisto_preselection = counter.getCutFlow("CutFlow/");
    wk()->addOutput(cutflowHisto_preselection);
  }

  Info("finalize()", "Processed events         = %li", m_eventCounter);
  Info("finalize()", "Passed nominal selection = %li", m_eventCountPassed);

  return EL::StatusCode::SUCCESS;
} // finalize

EL::StatusCode AnalysisReader::histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.

  Info("histFinalize()", "Finalizing histograms.");

  // a Tree with the luminosity
  m_histSvc->BookFillTree("MetaData", "luminosity", &m_luminosity, "m_luminosity");

  m_histSvc->Write(wk());

  return EL::StatusCode::SUCCESS;
} // histFinalize

EL::StatusCode AnalysisReader::OverlapPassThrough (const xAOD::ElectronContainer *electrons,
                                                   const xAOD::PhotonContainer   *photons,
                                                   const xAOD::MuonContainer     *muons,
                                                   const xAOD::TauJetContainer   *taus,
                                                   const xAOD::JetContainer      *jets,
                                                   const xAOD::JetContainer      *fatJets)  {
  if (m_eventCounter < 10) {
    Warning("OverlapPassThrough()", "Passing all objects to the OR. This method should be overridden.");
  }

  if (electrons) {
    for (const xAOD::Electron *elec : *electrons) {
      Props::passOR.set(elec, true);
    }
  }

  if (photons) {
    for (const xAOD::Photon *photon : *photons) {
      Props::passOR.set(photon, true);
    }
  }

  if (muons) {
    for (const xAOD::Muon *muon : *muons) {
      Props::passOR.set(muon, true);
    }
  }

  if (taus) {
    for (const xAOD::TauJet *tau : *taus) {
      Props::passOR.set(tau, true);
    }
  }

  if (jets) {
    for (const xAOD::Jet *jet : *jets) {
      Props::passOR.set(jet, true);
    }
  }

  if (fatJets) {
    for (const xAOD::Jet *jet : *fatJets) {
      Props::passOR.set(jet, true);
    }
  }

  return EL::StatusCode::SUCCESS;
} // OverlapPassThrough

float AnalysisReader::computeBTagSFWeight (const std::vector<const xAOD::Jet*> &signalJets, const std::string &authorName)
{
  if (m_debug) Info("computeBTagSFWeight()", "Called.");

  if (!m_isMC) {
    return 1;
  }
 
  bool truth_tag=false;
  m_config->getif<bool>("doTruthTagging", truth_tag);
  float weight = 1.0;
  m_bTagTool->setJetAuthor(authorName);
  bool doWeightVar{m_bTagTool->doWeightVar()},
    isNominal{(m_currentVar == "Nominal")};

    if (doWeightVar && !isNominal) m_bTagTool->setWeightVar(false);
    auto btagEffSFs = (truth_tag?m_bTagTool->computeEventWeight_truthTag(signalJets,m_config):m_bTagTool->computeEventWeight(signalJets));
    m_bTagTool->setWeightVar(doWeightVar);

    weight = btagEffSFs["Nominal"];

    //only write syst histos if nominalOnly is set to false
    bool nominalOnly = false;
    m_config->getif<bool>("nominalOnly", nominalOnly);
    if (!nominalOnly && isNominal && doWeightVar) {
      for (auto effSF : btagEffSFs) {
        std::string varName = effSF.first;
        if (varName == "Nominal") continue;
        if (authorName != "") {
          auto n = varName.rfind("__1up");
          if (n == std::string::npos) n = varName.rfind("__1down");
          if (n == std::string::npos) varName += "_" + authorName;
          else varName.insert(n, "_" + authorName);
        }
        m_weightSysts.push_back({ varName, effSF.second / weight });

        // Info("computeBTagSFWeight", "Relative weight for %s = %f", effSF.first.c_str(), effSF.second/weight);
      }
    }

    // Info("computeBTagSFWeight()", "Event weight = %f", weight);
    return weight;
} // computeBTagSFWeight


void AnalysisReader::getNeutrinoPz (double MET, double MET_phi, const TLorentzVector &lepVec, double &nu_pz1, double &nu_pz2) {
  // This code gives results identical (except using a slightly different mW) to:
  // https://svnweb.cern.ch/trac/atlasphys-exa/browser/Physics/Exotic/Analysis/DibosonResonance/
  // Data2015/VV_JJ/Code/trunk/CxAODReader_DB/Root/AnalysisReader_DB.cxx?rev=239154#L812
  double mW = 80385.; // PDG 2014

  double el  = lepVec.E();
  double ptl = lepVec.Pt();
  double pzl = lepVec.Pz();

  TLorentzVector metVec;

  metVec.SetPtEtaPhiM(MET, 0, MET_phi, 0);

  double mu    = 0.5 * mW * mW + ptl * MET * cos(lepVec.DeltaPhi(metVec));
  double delta = (mu * mu * pzl * pzl / (pow(ptl, 4))) - (el * el * MET * MET - mu * mu) / (ptl * ptl);

  if (delta < 0) {
    delta = 0;
  }

  nu_pz1 = (mu * pzl / (ptl * ptl)) + sqrt(delta);
  nu_pz2 = (mu * pzl / (ptl * ptl)) - sqrt(delta);
} // getNeutrinoPz

double AnalysisReader::getNeutrinoPz (double MET, double MET_phi, const TLorentzVector &lepVec, bool min) {
  double nu_pz1;
  double nu_pz2;

  getNeutrinoPz(MET, MET_phi, lepVec, nu_pz1, nu_pz2);

  if (fabs(nu_pz1) < fabs(nu_pz2)) {
    if (min) return nu_pz1;

    return nu_pz2;
  }

  if (min) return nu_pz2;

  return nu_pz1;
} // getNeutrinoPz

TLorentzVector AnalysisReader::getNeutrinoTLV (const TLorentzVector &metVec,const TLorentzVector &lepVec, bool min) {
  double pXNu = metVec.Px();
  double pYNu = metVec.Py();
  double pZNu = getNeutrinoPz(metVec.Pt(), metVec.Phi(), lepVec, min);
  double eNu = sqrt(pXNu*pXNu + pYNu*pYNu + pZNu*pZNu);
  TLorentzVector nuVec(pXNu, pYNu, pZNu, eNu);
  return nuVec;
} // getNeutrinoTLV

int AnalysisReader::getGAHeavyFlavHadronLabels_leadPt(const xAOD::Jet * jet) const {

  std::vector<std::pair<int, float> > hadrons;
  getGAHeavyFlavHadronLabels_PtSort( jet , hadrons );
  
  if ( hadrons.size() > 0 ) return hadrons[0].first;

  return 0;
  
}

void AnalysisReader::getGAHeavyFlavHadronLabels_PtSort(const xAOD::Jet * jet, int & jetflav1, int & jetflav2) const {

  std::vector<std::pair<int, float> > hadrons;
  getGAHeavyFlavHadronLabels_PtSort( jet , hadrons );
  
  if ( hadrons.size() > 1 ) {
    jetflav1 = hadrons[0].first;
    jetflav2 = hadrons[1].first;
  }
  else if ( hadrons.size() > 0 ) {
    jetflav1 = hadrons[0].first;
    jetflav2 = 0; 
  }
  else {
    jetflav1 = 0;
    jetflav2 = 0;
  }
  
}

void AnalysisReader::getGAHeavyFlavHadronLabels_PtSort(const xAOD::Jet * jet, std::vector<std::pair<int, float> > & hadrons) const {
  
  // get hadrons
  const std::string labelB   = "GhostBHadronsFinal";
  const std::string labelC   = "GhostCHadronsFinal";
  const std::string labelTau = "GhostTausFinal"; 
  std::vector<const xAOD::IParticle *> ghostB;
  std::vector<const xAOD::IParticle *> ghostC;
  std::vector<const xAOD::IParticle *> ghostTau;
  jet->getAssociatedObjects<xAOD::IParticle>( labelB   , ghostB   );
  jet->getAssociatedObjects<xAOD::IParticle>( labelC   , ghostC   );
  jet->getAssociatedObjects<xAOD::IParticle>( labelTau , ghostTau ); 

  // fill vector
  hadrons.clear();
  hadrons.reserve( ghostB.size() + ghostC.size() + ghostTau.size() );
  for (const xAOD::IParticle * had : ghostB  ) hadrons.push_back( std::make_pair(  5 , had->pt() ) );
  for (const xAOD::IParticle * had : ghostC  ) hadrons.push_back( std::make_pair(  4 , had->pt() ) );
  for (const xAOD::IParticle * had : ghostTau) hadrons.push_back( std::make_pair( 15 , had->pt() ) );

  // sort vector
  std::sort(hadrons.begin(), hadrons.end(), sortPtGAHeavyFlavHadrons);
  
}


bool AnalysisReader::sortPtGAHeavyFlavHadrons(std::pair<int, float> had1, std::pair<int, float> had2) {

  // sort by hadron pt
  // if two hadrons have same pt, sort by flavour (B before C before Tau)
  if      ( had1.second > had2.second ) return true;
  else if ( had1.second < had2.second ) return false;
  else if ( abs(had1.first) == 5 ) return true;
  else if ( abs(had2.first) == 5 ) return false;
  else if ( abs(had1.first) == 4 ) return true;
  else if ( abs(had2.first) == 4 ) return false;
  return true;

}


EL::StatusCode AnalysisReader::getMuonInJetCorrTLV(const xAOD::Jet * jet, TLorentzVector & veccor, bool isFatJet, const std::string & muJetCorr) {
  
  // We don't store the b-get energy corrected jets for systematic shifts 
  // it is not just muon-in-jet, but also Regression and versions of PtReco
  // so we have to work them out from the nominal
  
  // find jet from nominal container                                                                                                                          
  const xAOD::Jet * jetnom = 0;  

  // retrieve nominal jets
  const xAOD::JetContainer * jets_nominal = isFatJet ? m_fatJetReader->getObjects("Nominal") : m_jetReader->getObjects("Nominal");
  // check if the current variation is in fact the nominal
  if ( (isFatJet && m_fatJets == jets_nominal) || (!isFatJet && m_jets == jets_nominal) ) {
    jetnom=jet;
  } else {
    // the current jet container is not the nominal
    for (const xAOD::Jet * jettest : *jets_nominal) {
      if ( Props::partIndex.get(jettest) == Props::partIndex.get(jet) ) {
	jetnom = jettest;
	break;
      }
    }
  }
  // if we don't have a nominal jet, we return as failure
  if ( ! jetnom ) return EL::StatusCode::FAILURE;

  //veccor is evaluated below as the difference for the nominal (no systematic) jet between the desired correction and the Nominal (GSC+in-situ) correction
  //the code below should have only jetnom->, not jet->
  if(muJetCorr=="OneMu"){
    //muon-in-jet correction
    xAOD::JetFourMom_t HVecMuTmp = jetnom->jetP4(muJetCorr);
    veccor.SetPx(HVecMuTmp.px());
    veccor.SetPy(HVecMuTmp.py());
    veccor.SetPz(HVecMuTmp.pz());
    veccor.SetE(HVecMuTmp.e());
    if(m_jetCorrectionOld){
      veccor -= jetnom->p4();
    } else{
      //do nothing as the new style is to store directly the difference with respect to GSC+in-situ
    }
  } else {
    //correction is other than OneMu
    if (isFatJet){
      Warning("AnalysisReader::getMuonInJetCorrTLV(...)","You have a fat jet and asking for %s correction, which is not defined for fatjet.", muJetCorr.c_str());
    }
    if(muJetCorr=="Regression"){
      //MVA-based regression
      //applied as a multiplicative factor
      //directly on top of GSC(+in-situ), called nominal, which is retrieved with jetnom->pt() directly.
      veccor = jetnom->p4();
      float initialPt=veccor.Pt();
      if(initialPt<0.1){
	Warning("AnalysisReader::getMuonInJetCorrTLV(...)","for Regression initialPt=%f",initialPt);
      }
      float correctedPt=Props::ptCorr.get(jet);
      float factor=correctedPt/initialPt;
      veccor*=factor;
      veccor -= jetnom->p4();
    } else if(muJetCorr.find("PtReco") != std::string::npos){
      //one of the several PtReco corrections
      //applied as a multiplicative factor
      //on top of the OneMu correction, which is retrieved with jetnom->jetP4("OneMu") 
      //in old style it was the correction itself
      //in new style it is the difference with with respect to GSC+in-situ
      xAOD::JetFourMom_t HVecMuTmp = jetnom->jetP4("OneMu");
      veccor.SetPx(HVecMuTmp.px());
      veccor.SetPy(HVecMuTmp.py());
      veccor.SetPz(HVecMuTmp.pz());
      veccor.SetE(HVecMuTmp.e());
      if(!m_jetCorrectionOld){
        veccor += jetnom->p4();
      }
      float initialPt=veccor.Pt();
      if(initialPt<0.1){
	Warning("AnalysisReader::getMuonInJetCorrTLV(...)","for Regression initialPt=%f",initialPt);
      }
      float correctedPt=jetnom->jetP4(muJetCorr).Pt();
      if(!m_jetCorrectionOld){
	correctedPt+=jetnom->p4().Pt();
      }
      float factor=correctedPt/initialPt;
      veccor*=factor;
      veccor -= jetnom->p4();
    } else {
      Warning("AnalysisReader::getMuonInJetCorrTLV(...)","You are asking for %s correction, which is not defined", muJetCorr.c_str());
    }//end if muJetCorr is Regression or PtReco
  }//end if correction is OneMu or other
  
  //veccor is now evaluated below as the difference for the nominal (no systematic) jet between the desired correction and the Nominal (GSC+in-situ) correction
  //we add it on the top of the Nominal (GSC) of the jet with the desired systematic
  veccor += jet->p4();

  //now ready to return success
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode AnalysisReader::applyRtrkUncertFix(const xAOD::Jet * jet) {

  // This is a temporary fix to fat jet variations JET_Rtrk_baseline, JET_Rtrk_Modelling, JET_Rtrk_Tracking,
  // which were accidently applied twice in FatJetHandler since around January 2016 to May 11 2016.
  // This method will undo the second application of the variations by assuming that the same correction factor 
  // was applied twice. There will be non-closure effects in cases where the jet migrates from one bin in the
  // correction map to another after the first correction is applied (which means that two different correction
  // factors were applied in CxAODMaker). 

  // The correction done here is:
  // pT(variation) = pT(nominal) * sqrt( pT(2*variation) / pT(nominal) )
  // m(variation)  = m(nominal)  * sqrt( m(2*variation)  / m(nominal) )

  // The argument in this method is a jet from one of the problematic variations:
  // JET_Rtrk_Baseline 
  // JET_Rtrk_Modelling
  // JET_Rtrk_Tracking
  // ====> The user needs to make sure only fat jets from these variations are passed to this method !!! <=====
  // (The method will check - and return an error if not.)

  // This method will attach the corrected TLorentzVector ("RtrkUncertFix_TLV") which can be retrieved later with:
  //
  // TLorentzVector jetP4 = jet->auxdata<TLorentzVector>("RtrkUncertFix_TLV");
  // 
  // or:
  //
  // static SG::AuxElement::Accessor<TLorentzVector> acc("RtrkUncertFix_TLV");
  // TLorentzVector jetP4 = acc(*jet);                                                                                                                                                                                                 

  // check the current variation 
  if ( m_currentVar.find("JET_Rtrk") == std::string::npos ) {
    Error("applyRtrkUncertFix()","Ooops, make sure to only pass fat jets from variations: JET_Rtrk_[Baseline/Modelling/Tracking] !!");
    return EL::StatusCode::FAILURE;
  }

  // retrieve nominal jets  
  const xAOD::JetContainer * jets_nominal = m_fatJetReader->getObjects("Nominal");

  // find jet from nominal container
  const xAOD::Jet * jetnom = 0;
  for (const xAOD::Jet * jettest : *jets_nominal) {
    if ( Props::partIndex.get(jettest) == Props::partIndex.get(jet) ) {
      jetnom = jettest;
      break;
    }
  }

  // check if nominal jet is found
  if ( ! jetnom ) {
    Error("applyRtrkUncertFix()","Couldn't find nominal jet - this should never happen - something is wrong!!");
    return EL::StatusCode::FAILURE;
  }

  // pT correction squared
  double pTcorr2 = jet->pt()/jetnom->pt();
  if ( pTcorr2 < 0. ) {
    Error("applyRtrkUncertFix()","pT correction factor squared is negative : %f", pTcorr2);
    return EL::StatusCode::FAILURE;
  }

  // mass correction squared
  double Mcorr2 = jet->m()/jetnom->m();
  if ( Mcorr2 < 0. ) {
    Error("applyRtrkUncertFix()","mass correction factor squared is negative : %f", Mcorr2);
    return EL::StatusCode::FAILURE;
  }

  // store corrected values as decorations
  xAOD::JetFourMom_t jetP4_wrong = jet->jetP4();
  TLorentzVector jetP4_correct;
  jetP4_correct.SetPtEtaPhiM( jetnom->pt() * sqrt( pTcorr2 ) , jetP4_wrong.Eta() , jetP4_wrong.Phi() , jetnom->m() * sqrt( Mcorr2 ) );
  static SG::AuxElement::Decorator<TLorentzVector> dec("RtrkUncertFix_TLV");
  dec(*jet) = jetP4_correct;

  // std::cout << "nominal : " << jetnom->pt() << std::endl;
  // std::cout << "before  : " << jet->pt() << std::endl;
  // std::cout << "after   : " << (jet->auxdata<TLorentzVector>("RtrkUncertFix_TLV")).Pt() << std::endl;
  // std::cout << "correction factors (Pt,M) : " << sqrt(pTcorr2) << " " << sqrt(Mcorr2) << std::endl;

  // all is good
  return EL::StatusCode::SUCCESS;

}
