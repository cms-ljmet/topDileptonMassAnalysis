CXX=g++

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)

ROOTIL      = $(ROOTCFLAGS) $(ROOTLIBS)

LIBS            = libObjects.so

TopAnalysis:
	$(CXX) $(ROOTIL) -L ${TOOLSDIR}/lib  -l Objects -l MassAMWT -l LHAPDF -l JetMETObjects -I ${TOOLSDIR} -I ${LHAPDFINC} -O2 -Wall -o TopAnalysis TopAnalysis.C

TopPartnersv4:
	$(CXX) $(ROOTIL) -L ${TOOLSDIR}/lib  -l Objects -l MassAMWT -l LHAPDF -l JetMETObjects -I ${TOOLSDIR} -I ${LHAPDFINC} -O2 -Wall -o TopPartnersv4 TopPartnersv4.C

Yield:
	$(CXX) $(ROOTIL) -L ${TOOLSDIR}/lib  -l Objects -I ${TOOLSDIR} -o Yield Yield.C

YieldSA:
	$(CXX) $(ROOTIL) -l RooFit -l RooFitCore -O2 -Wall -o YieldSA YieldSA.C

YieldSA_pdf:
	$(CXX) $(ROOTIL) -L ${TOOLSDIR}/lib -l RooFit  -l LHAPDF -I ${LHAPDFINC} -O2 -Wall -o YieldSA_pdf YieldSA_pdf.C

combineBins1D:
	$(CXX) $(ROOTIL) -O2 -Wall -o combineBins1D combineBins1D.C

DrawOneHistoFR:
	$(CXX) $(ROOTIL) -O2 -Wall -o DrawOneHistoFR DrawOneHistoFR.C

DrawOneHistoDDB:
	$(CXX) $(ROOTIL) -O2 -Wall -o DrawOneHistoDDB DrawOneHistoDDB.C

DrawOneHistoWithSyst:
	$(CXX) $(ROOTIL) -O2 -Wall -o DrawOneHistoWithSyst DrawOneHistoWithSyst.C

EventYields:
	$(CXX) $(ROOTIL) -O2 -Wall -o EventYields EventYields.C

Systematics:
	$(CXX) $(ROOTIL) -O2 -Wall -o Systematics Systematics.C

PromptLnTProb:
	$(CXX) $(ROOTIL) -O2 -Wall -o PromptLnTProb PromptLnTProb.C

compareTwo:
	$(CXX) $(ROOTIL) -O2 -Wall -o compareTwo compareTwo.C

GetFakeRate:
	$(CXX) $(ROOTIL) -O2 -Wall -o GetFakeRate GetFakeRate.C

IsoCheck:
	$(CXX) $(ROOTIL) -O2 -Wall -o IsoCheck IsoCheck.C

ErrorPropagation:
	$(CXX) $(ROOTIL) -O2 -Wall -o ErrorPropagation ErrorPropagation.C

clean: 
	rm *.o *.so *.d


all: TopPartnersv4 DataDrivenOS LooseTopPartnersv4

.PHONY: TopAnalysis
.PHONY: TopPartnersv4
.PHONY: Yield
.PHONY: YieldSA
.PHONY: YieldSA_pdf
.PHONY: combineBins1D
.PHONY: DrawOneHistoFR
.PHONY: DrawOneHistoDDB
.PHONY: DrawOneHistoWithSyst
.PHONY: EventYields
.PHONY: Systematics
.PHONY: PromptLnTProb
.PHONY: compareTwo
.PHONY: GetFakeRate
.PHONY: IsoCheck
.PHONY: ErrorPropagation
