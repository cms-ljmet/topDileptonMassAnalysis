
import cPickle
import re
import ROOT as root

inFileName = "MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl"
outFileName = "MuonEfficiencies2011.root"

file1 = open(inFileName,"r")

data = cPickle.load(file1)

def myCompare(x,y):
  xnum = float(re.sub(r"_.*","",x))
  ynum = float(re.sub(r"_.*","",y))
  return cmp(xnum,ynum)

def printAll():
    #for studyName in data:
    for studyName in data:
	print("getting study name: "+studyName)
	#study is a study of a particlar effect on a certain dataset
	study = data[studyName]
	for plotName in study:
	    print("plotting: "+plotName)
	    #plot is an efficiency plot
	    plot = study[plotName]
	#  graph = root.TGraphAsymmErrors()
	#  graph.SetName(studyName+plotName)
	#  i = 0
	    for varName in sorted(plot,cmp=myCompare):
	    #var is a binned variable
		var = plot[varName]
		SF = var['data/mc']
        	print varName, plotName, " ", SF['efficiency_ratio'], SF['err_low'], SF['err_hi']
	#	    val = SF[effName]
	#	    print val

def spitOut(studyName):
    #for studyName in data:
    print("getting study name: "+studyName)
    #study is a study of a particlar effect on a certain dataset
    study = data[studyName]
    for plotName in study:
	print("pltting: "+plotName)
	#plot is an efficiency plot
	plot = study[plotName]
    #  graph = root.TGraphAsymmErrors()
    #  graph.SetName(studyName+plotName)
    #  i = 0
	for varName in sorted(plot,cmp=myCompare):
	#var is a binned variable
	    var = plot[varName]
	    SF = var['data/mc']
            print varName, plotName, " ", SF['efficiency_ratio'], SF['err_low'], SF['err_hi']

#spitOut('Loose')
spitOut('combRelIsoPF04dBeta<02_Loose')
#printAll()

#rootFile.Close()
