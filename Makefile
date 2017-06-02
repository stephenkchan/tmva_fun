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

_CMNOBJ = AtlasStyle.o HistoTransform.o bdt_base.o bdt_trainer.o bdt_testing.o bdt_validate.o bdt_ranker.o rfli_rank.o rfli_plot.o
CMNOBJ = $(patsubst %,$(ODIR)/%,$(_CMNOBJ))

_DEPS = bdt_base.h bdt_trainer.h bdt_testing.h bdt_validate.h bdt_ranker.h rfli_rank.h rfli_plot.h #sig_plotter.h headers that might change; others in principle, but those are at this point deprecated
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

SRC_EXT = C cpp cxx

define comp_objects
$(ODIR)/%.o: $$(SDIR)/%.$1 $$(DEPS)
	$$(CC) $$(CFLAGS) $$< $$(LIBS) -o $$@
endef
$(foreach EXT,$(SRC_EXT),$(eval $(call comp_objects,$(EXT))))

all: plot_rfli fin_rfli train_rfli rfli_single

plot_rfli: $(PDIR)/plot_rfli.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

fin_rfli: $(PDIR)/fin_rfli.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

rfli_single: $(PDIR)/rfli_single.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

train_rfli: $(PDIR)/train_rfli.cxx $(OGCM) $(CMNOBJ)
	$(CC) -o $(EDIR)/$@ $^ $(LIBS) 

clean:
	rm -rf $(ODIR)/*.o $(EDIR)/* *.o *~ core*