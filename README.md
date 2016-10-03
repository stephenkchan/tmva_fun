----------------------------------------------------------------------------------
                                   TMVA_FUN
                     BDT Training with CxAODReader Outputs
		              S. Chan 2016-03-15
			        s.chan@cern.ch
----------------------------------------------------------------------------------
The tmva_fun package consists of four main worker classes:
    --bdt_trainer
    --bdt_validate
    --bdt_ranker
    --bdt_tester
As well as a variety of wrapper classes for analysis-specific ranking/training/plot making.

----------------------------------------------------------------------------------
I.  SETUP

    1. Setup the local workspace with ROOT; TMVA must be installed:
       a. % setupATLAS; lsetup root #on lxplus
       b. % source ~stephenchan/atlas_setup.sh; lsetup root #on herophysics

    2. Construct the local file architecture if you haven't copied the code properly
       --The architecture can be read off from the Makefile
       --(If it would be convenient for you, contact the developer to make a setup script)

    3. % make

    4. Executables are now ready for use.

----------------------------------------------------------------------------------
II. TRAINING
    Training is handled by bdt_trainer. As the name suggests, training is done with a BDT.

    BDT Settings: "!H:!V:NTrees=200:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.15:SeparationType=GiniIndex:nCuts=100:PruneMethod=NoPruning:MinNodeSize=1.5%"
    Event Splitting: "SplitMode=Alternate" (events with "EventNumber%3==0" are skipped to be set aside for testing)

    Its constructor looks like:

    bdt_trainer(inputs (data-MVA/*root files from the CxAODReader) directory, 
    		outputs directory (where TMVA's .xml/.C weights files and .root files will go),
		identifier for files (to avoid duplication)
		vector of TString's containing variable names as in the "Nominal" TTree of input files,
		int for btagging fixed point efficiency cuts [default 70],
		charm ratio for MV2c (00, 10, 20 currently supported, with working points only available for 20) [default 20],
		bool for continuous b-tagging inputs [default true])

    Output file names are of the form {output_dir}/TMVA{identifier}-MV2c{charm ratio}-beff{btagging working point}-{- separated variable names in order for BDT training}.root (similar for TMVAClassification*weights.(xml|C)).
    Training is done with bdt_trainer::train(int bg type [def -99 for Z+jets+ttbar],bool for setting aside validation events [def true])

----------------------------------------------------------------------------------
III. VALIDATION
     bdt_validate takes the same arguments as and is derived from bdt_trainer plus an additional argument for overwriting [default false].  
     
     bdt_validate's main purpose is to provide significances to bdt_ranker.  It does this by:
         1. Fetching the training output .root file (if non-existent or overwrite set to true, will run training)
         2. Fetching the training distributions "MVA_BDT_Train_[SB]" from the file
	 3. Calculate transformation F bins based on training distributions
	 4. Fill contents of validation (i.e. TMVA "testing") distributions "MVA_BDT_[SB]" into transformation F histograms
	 5. Calculate optimal bin-wise significance

----------------------------------------------------------------------------------
IV. RANKING
    Ranking is done with bdt_ranker; it makes use of bdt_trainer.  Ranking is done iteratively.  In addition to specifying input and output directories, (ordered) starting variables and (unordered) variables to be ranked must be given.

    bdt_ranker(inputs,
 	       outputs (bdt_trainer),
	       TString vector of ordered variables already ranked (e.g. {"mBB","dRBB"}),
	       TString vector of unordered variables to be ranked) //tag is done

    bdt_ranker::do_ranking(string tag for log file) does the ranking.  Log output names look like: {output_dir}/ranking_log{tag}.txt

----------------------------------------------------------------------------------
IV. TESTING
    bdt_tester takes the same arguments of and is derived from bdt_trainer.  Its functionality is similar to bdt_validate.

    The only difference from bdt_validate is that instead of making use of TMVA's testing samples, the .xml weights file from training is used (if this or the .root file does not exist, training is called).  After transformation F, testing BDT distributions are instead filled by using the Reader on the original CxAODReader data-MVA/*root files on events that pass event selection with EventNumber%3==0.


----------------------------------------------------------------------------------
VI. LI/B-TAGGING SPECIFIC CLASSES
*****under construction*****
     li_ranking and li_plotter do ranking and plotting for the different classes of LI variables defined.

     Main executables:
     % exec/ranking_variables [tmva_in] [tmva_out] [variable set int] [overwrite?] #typical run time of about 15-20 minutes for 2+10 variables on pleiades SLURM
     % exec/val_rank_plots [tmva_in] [tmva_out] [plot directory] 

     % python/admiral.py [tmva_in] [identifier tag for output (tmva_out is tmva_in/{tag})] [do_batch] [overwrite?]
     Curently, admiral.py will work on herophysics and lxplus without modification.
	 
