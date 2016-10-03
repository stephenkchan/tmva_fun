#TMVA_FUN Makefile
IDIR = include
SDIR = src
ODIR = obj
EDIR = exec
PDIR = programs
CFLAGS=-c -Wall  #put these things in $(insert here) and it's like you typed them
INCTMVA  = -I$(ROOTSYS)/tmva/test/
ROOTFLAGS = $(shell root-config --cflags) -fPIC
ROOTLIBS = $(shell root-config --libs)  
CC=g++ $(ROOTFLAGS) -g -I$(IDIR) $(INCTMVA)

LIBS = -lGraf -lHistPainter -lTMVA $(ROOTLIBS) #-lRooFit -lRooFitCore -lMinuit  -lgsl

_CMNOBJ = AtlasStyle.o HistoTransform.o bdt_base.o bdt_trainer.o bdt_testing.o bdt_validate.o bdt_ranker.o li_ranking.o li_plotter.o std_ranking.o btag_ranking.o std_plotter.o btag_plotter.o #sig_plotter.o 
CMNOBJ = $(patsubst %,$(ODIR)/%,$(_CMNOBJ))

_DEPS = bdt_base.h bdt_trainer.h bdt_testing.h bdt_validate.h bdt_ranker.h li_ranking.h li_plotter.h std_ranking.h btag_ranking.h std_plotter.h btag_plotter.h #sig_plotter.h headers that might change; others in principle, but those are at this point deprecated
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

SRC_EXT = C cpp cxx

define comp_objects
$(ODIR)/%.o: $$(SDIR)/%.$1 $$(DEPS)
	echo 'Lucky charms!'
	$$(CC) $$(CFLAGS) $$< $$(LIBS) -o $$@
endef
$(foreach EXT,$(SRC_EXT),$(eval $(call comp_objects,$(EXT))))

all: val_rank_plots rank_variables rank_btag plot_btag finish_btag finish_std finish_li tester #exec_tmva_li_train plots rank_std plot_std 

plot_btag: $(PDIR)/plot_btag.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

plot_std: $(PDIR)/plot_std.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

val_rank_plots: $(PDIR)/val_rank_plots.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

finish_li: $(PDIR)/finish_li.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

finish_btag: $(PDIR)/finish_btag.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

finish_std: $(PDIR)/finish_std.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

rank_btag: $(PDIR)/rank_btag.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

rank_std: $(PDIR)/rank_std.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

rank_variables: $(PDIR)/rank_variables.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

tester: $(PDIR)/tester.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

exec_tmva_li_train: $(PDIR)/exec_tmva_li_train.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

plots: $(PDIR)/plots.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

clean:
	rm -rf $(ODIR)/*.o $(EDIR)/* *.o *~ core*