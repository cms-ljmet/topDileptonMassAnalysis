#!/bin/tcsh
# Parameters of the script:
# directory
# output file name

# source  drawWithSystPlotsLogDDB.csh pass2 final.root

make DrawOneHistoWithSyst
set startDir=$cwd
cd $argv[1]
set rootFile=$cwd/$argv[2]
echo $cwd

#$startDir/DrawOneHistoWithSyst $rootFile Top All Top amwt	0 "AMWT mass" 10 0 440

#$startDir/DrawOneHistoWithSyst $rootFile TwoL2J1BmetZ All 2L_2J1BmetZ amwt	0 "AMWT mass" 10 0 440

#$startDir/DrawOneHistoWithSyst $rootFile TwoL2Jmet All 2L_2Jmet nBJets  1 "Number of b-tagged jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoL2Jmet All 2L_2Jmet amwt	0 "AMWT mass" 20 0 600
#$startDir/DrawOneHistoWithSyst $rootFile TwoL2J All 2L_2J amwt	0 "AMWT mass" 10 0 600
###$startDir/DrawOneHistoWithSyst $rootFile TwoL2J All 2L_2J amwt	0 "AMWT mass" 20 0 800
#$startDir/DrawOneHistoWithSyst $rootFile TwoL2J1BmetZ All 2L_2J1BmetZ MET     1 "MET [GeV]" 2 0 300
##$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS MuMu 2L_AS Lep1Pt  1 "Leading lepton p_{T} [GeV/c]" 5
#
#cd ..
#exit(1)
#
#DrawOneHisto parameter list:
#1: Input directory in root file
#2: Channel (ElEl, ElMu or MuMu   
#3: Prefix for histograms 
#4: Name of histogram
#5: Log or not log 
#6: x-Axis name
#7: Optional: rebin by this number

#mkdir 8TeV_test_SS/
#mkdir AllPDF


##Two good leptons only -- no charge requirement
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS Lep1Pt  1 "Leading lepton p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS Lep2Pt  1 "Second lepton p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS Lep3Pt  1 "Third lepton p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS Jet1Pt  1 "Leading jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS Jet2Pt  1 "Second jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS LepInvM 1 "Invariant Mass of Lepton Pair [GeV/c^{2}]" 4
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS nJets   1 "Number of Jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS nBJets  1 "Number of b-tagged jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS HT      1 "H_{T} [GeV]" 20
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS MET     1 "MET [GeV]"
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS nLept   1 "Number of Leptons"
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS sumPtL  1 "Sum of lepton p_{T} [GeV/c]" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS sumPtJ  1 "Sum of jet p_{T} [GeV/c]" 2
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS minMlb  1 "Min(m_{lb}})" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS Mlb     1 "m_{lb}}" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS dRlept  1 "#delta R (l1,l2)" 2
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS MinLepCATopdR	1 "Min (#delta R (l,top-jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS MinLepCAWdR	1 "Min (#delta R (l,CAW Jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS MinLepAK5dR	1 "Min (#delta R (l,AK5 jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS MinLepJetdR 	1 "Min (#delta R (l,jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS MinLepBJetdR	1 "Min (#delta R (l,b-jet))" 2
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS CAWJet1Pt  1 "Leading CAW Jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS CAWJet2Pt  1 "Second CAW Jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS nCAWJets   1 "Number of CAW Jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS CAWBtag    1 "CAW Jet CSV discriminant"
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS CATopJet1Pt  1 "Leading CATop Jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS CATopJet2Pt  1 "Second CATop Jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS nCATopJets   1 "Number of CATop Jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS CATopBtag    1 "CATop Jet CSV discriminant"
#
##$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsAS All 2L_AS amwt    1 "AMWT mass" 5
#
#
#mkdir TwoLeptons
#mv  *_2L_AS*.eps *_2L_AS*.png            TwoLeptons
#cd TwoLeptons
## diowroot.pl -o index.html -t 'Two leptons, without charge requirements'
#cd ..
##cp indexFiles/i2L_AS.html TwoLeptons/index.html
##mv h_2L_AS*pdf            AllPDF/

#Two same-sign/opposite-sign leptons, with eventual vetos
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV Lep1Pt  1 "Leading lepton p_{T} [GeV/c]" 20 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV Lep2Pt  1 "Second lepton p_{T} [GeV/c]" 20
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV Lep3Pt  1 "Third lepton p_{T} [GeV/c]" 5 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV Jet1Pt  1 "Leading jet p_{T} [GeV/c]" 10 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV Jet2Pt  1 "Second jet p_{T} [GeV/c]" 10 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV LepInvM 1 "Invariant Mass of Lepton Pair [GeV/c^{2}]" 10 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV nJets   1 "Number of Jets" &
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV nBJets  1 "Number of b-tagged jets" 
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV HT      1 "H_{T} [GeV]" 20 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV MET     1 "MET [GeV]"
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV nLept   1 "Number of Leptons"&
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV sumPtL  1 "Sum of lepton p_{T} [GeV/c]" 2
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV minMlb  1 "Min(m_{lb}})" 2&
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV Mlb     1 "m_{lb}}" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV dRlept  1 "#delta R (l1,l2)" 2&
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV MinLepCATopdR	1 "Min (#delta R (l,top-jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV MinLepCAWdR	1 "Min (#delta R (l,CAW Jet))" 2 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV MinLepAK5dR	1 "Min (#delta R (l,AK5 jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV MinLepJetdR 	1 "Min (#delta R (l,jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV MinLepBJetdR	1 "Min (#delta R (l,b-jet))" 2
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV CAWJet1Pt  1 "Leading CAW Jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV CAWJet2Pt  1 "Second CAW Jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV nCAWJets   1 "Number of CAW Jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV CAWBtag    1 "CAW Jet CSV discriminant"
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV CATopJet1Pt  1 "Leading CATop Jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV CATopJet2Pt  1 "Second CATop Jet p_{T} [GeV/c]" 5
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV nCATopJets   1 "Number of CATop Jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV CATopBtag    1 "CATop Jet CSV discriminant"
#
##$startDir/DrawOneHistoWithSyst $rootFile TwoLeptonsZQV All 2L_ZQV amwt    1 "AMWT mass" 5
#
#mkdir                     TwoLeptonsSigned
#
#mv *_2L_ZQV*.eps *_2L_ZQV*.p*            TwoLeptonsSigned
#cd  TwoLeptonsSigned
# diowroot.pl -o index.html -t 'Two leptons, with charge requirements' &
#cd ..


#foreach ch (2J 2Jmet 2J1B 2J1Bmet 2J2Bmet 2J2B 3J 3Jmet 3J1B 3J1Bmet 3J1BmetHT 3J2B 3J2Bmet 3J2BmetHT)

foreach ch (2J 2Jmet 2J1B 2J1Bmet 2J1BmetZ)

$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} Lep1Pt  1 "Leading lepton p_{T} [GeV/c]" 10 &
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} Lep2Pt  1 "Second lepton p_{T} [GeV/c]" 10
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} Lep3Pt  1 "Third lepton p_{T} [GeV/c]" 10 &
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} Jet1Pt  0 "Leading jet p_{T} [GeV/c]" 2 0 500
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} Jet2Pt  0 "Second jet p_{T} [GeV/c]" 2 0 500 &
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} LepInvM 1 "Invariant Mass of Lepton Pair [GeV/c^{2}]" 5 0 250
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} nJets	1 "Number of Jets" &
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} nBJets  1 "Number of b-tagged jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} HT	1 "H_{T} [GeV]" 20 &
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} MET	1 "MET [GeV]" 2 0 500

$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} nLept	1 "Number of Leptons" &
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} sumPtL  1 "Sum of lepton p_{T} [GeV/c]" 5 &

#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} maxMlb  1 "Max(m_{lb})" 2 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} minMlb  1 "min(m_{lb})" 2 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} Mlb	1 "m_{lb}}" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} dRlept  1 "#delta R (l1,l2)" 2 &

#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} MinLepCATopdR	    1 "Min (#delta R (l,top-jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} MinLepCAWdR 1 "Min (#delta R (l,CAW Jet))" 2 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} MinLepAK5dR 1 "Min (#delta R (l,AK5 jet))" 2
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} MinLepJetdR	    1 "Min (#delta R (l,jet))" 2 &
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} MinLepBJetdR	    1 "Min (#delta R (l,b-jet))" 2
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} CAWJet1Pt  1 "Leading CAW Jet p_{T} [GeV/c]" 10
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} CAWJet2Pt  1 "Second CAW Jet p_{T} [GeV/c]" 10
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} nCAWJets   1 "Number of CAW Jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} CAWBtag    1 "CAW Jet CSV discriminant" 5
#
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} CATopJet1Pt  1 "Leading CATop Jet p_{T} [GeV/c]" 10
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} CATopJet2Pt  1 "Second CATop Jet p_{T} [GeV/c]" 10
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} nCATopJets   1 "Number of CATop Jets"
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} CATopBtag    1 "CATop Jet CSV discriminant" 5

$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} amwt	0 "AMWT mass" 10 100 440

#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} leptonMETSum	    1 "leptonMETSum    " 10
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} leptonJetsSum	    1 "leptonJetsSum   " 10
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} jetsMETSum	    1 "jetsMETSum      " 10
#$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} leptonJetsMETSum    1 "leptonJetsMETSum" 10
$startDir/DrawOneHistoWithSyst $rootFile TwoL${ch} All 2L_${ch} nPV    0 "Nbr PV"


wait
mkdir TwoLeptonsSigned${ch}
mv *_2L_${ch}*.p*             TwoLeptonsSigned${ch}
cd  TwoLeptonsSigned${ch}
 diowroot.pl -o index.html  &
cd ..

end

cd $startDir
