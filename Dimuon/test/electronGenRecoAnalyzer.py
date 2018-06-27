import ROOT as rt
import sys,os
import argparse
import math

from DataFormats.FWLite import Events, Handle, Runs

# Function to check for (anti)muons that we are interested in
def getElectron(ele, debug):
    electron = ele
    if (debug):
        print ele.pdgId()
        print ele.status()
        print "Number of muon mothers: %d" %ele.numberOfMothers()
        #    numloopMuons+=1
        print "Muon mother is: %d \t with status: %d" %(ele.mother().pdgId(), ele.mother().status())
    nu = ele

    motherId = nu.pdgId()
    print "nuId: %d" %nu.pdgId()
    while (abs(motherId) > 6):
        nu = nu.mother()
        motherId = nu.pdgId()
        if (args.debug):
            print "motherId: %d" %motherId
        if (abs(motherId) > 13):
            electron = 0
        if (abs(motherId) < 7):
            electron = ele
    return electron

#Function to print Muon Event information
def printGenInfo(numLoopTot, numElectrons, numAntiElectrons):
    pairFound = 0
    print "The number of events is: %d" %numLoopTot
    if (numLoopTot != numAntiElectrons):
        print "Look HERE!"
        numLoopTot+=1
        print "The number of events is: %d" %numLoopTot
    if (numElectrons != numAntiElectrons):
        print "number of Muons different from antimuons!!!"
    else:
        pairFound = 1
    return pairFound
    print"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print "                                          pt:eta:phi:id:status"
    print "The electron     information in this event is: %f:%f:%f:%d:%d" %(electron.pt(),electron.eta(),electron.phi(),electron.pdgId(),electron.status())
    print "The antielectronon information in this event is: %f:%f:%f:%d:%d" %(antielectron.pt(),antielectron.eta(),antielectron.phi(),antielectron.pdgId(),antielectron.status())
    print"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
            #        numMuons+=1
#            print "The number of muons so far are: %d" %numMuons                    


parser = argparse.ArgumentParser()
#parser.add_argument("infile", help="Make plots for cosmics",type=str)
#parser.add_argument("--local", help="local file",action="store_true")
parser.add_argument("-d", "--debug", help="debugging information",action="store_true")
parser.add_argument("-i", help="Input miniAOD File", type=str)

args = parser.parse_args()


filename = args.i
with open(filename, "r") as ins:
    array = []
    for line in ins:
        array.append(line.rstrip('\n'))
eossrc = "root://cmsxrootd.fnal.gov//"

rt.gROOT.SetBatch(True)
# load FWLite C++ libraries
rt.gSystem.Load("libFWCoreFWLite.so");
rt.gSystem.Load("libDataFormatsFWLite.so");
rt.FWLiteEnabler.enable()

#Create output ROOT file

of = rt.TFile("genOutput.root","recreate")



#Instantiate Histograms

#histTkLayWMeas = rt.TH1D("histTkLayWMeas","TrackerLayersWithMeasurement",25,-0.5,24.5)
#histVdPxHits = rt.TH1D("histVdPxHits", "ValidPixelHits",25,-0.5,24.5)
#histVdMuHits = rt.TH1D("histVdMuHits", "ValidMuonHits", 25,-0.5,54.5);
#histPtRes = rt.TH1D("histPtRes", "Ptresolution", 25, -0.15, 0.80); 
#histDxy = rt.TH1D("histDxy", "Dxy", 25, -0.70, 0.70)
#histMatchMuStation = rt.TH1D("histMatchMuStation", "MatchedMuonStations", 20, -0.5, 19.5)
histPt =rt.TH1D("histPt", "Pt", 200, 0.0, 2000)
histEta = rt.TH1D("histEta", "#Eta", 50, -4.0, 4.0)
histPhi = rt.TH1D("histPhi", "#Phi", 30, -3.2, 3.2)


histEleMinusGenPt = rt.TH1D("histEleMinusGenPt","Pt",100,0.,500.)
histEleMinusGenEta = rt.TH1D("histEleMinusGenEta", "Eta",50,-3.5,3.5)
histEleMinusGenPhi = rt.TH1D("histEleMinusGentPhi", "Phi", 50,-3.15,3.15);
histElePlusGenPt = rt.TH1D("histElePlusGenPt","Pt",100,0.,500.)
histElePlusGenEta = rt.TH1D("histElePlusGenEta", "Eta",50,-3.5,3.5)
histElePlusGenPhi = rt.TH1D("histElePlusGentPhi", "Phi", 50,-3.15,3.15);
histInvariantMass = rt.TH1D("histInvariantMass","Mass Dist",300,0.,3000)

#Counters to track number of muons

totalPairs = 0
totalElectrons = 0
totalAntiElectrons = 0

#Loop over ROOT file list

for x in array:
    print "\n Now Reading file:"
    print x
    print "\n"
    f = rt.TNetXNGFile("%s/%s"%(eossrc,x),"read")
    events  = Events("%s/%s"%(eossrc,x))
    runs = Runs("%s/%s"%(eossrc, x))
                     
               
#Get Events Tree and create handles for GRN and RECO
    t = f.Get("Events")
    genHandle = Handle('std::vector<reco::GenParticle>')
    genLabel  = "prunedGenParticles"
    eleHandle = Handle('std::vector<pat::Electron>')
    eleLabel  = "slimmedElectrons"


    rt.TH1.SetDefaultSumw2()

#Counters to track number of muons in event and store particle object
    numloopElectrons = 0
    numElectrons = 0
    numloopAntiElectrons = 0
    numAntiElectrons = 0
    electron = 0
    antielectron = 0
#tmpMuCount = 0
#tmpAmuCount = 0
    numLoopTot = 0
#Event Loop
    for i,ev in enumerate(events):
        ev.getByLabel(genLabel, genHandle)
        genParticles = genHandle.product()
        ev.getByLabel(eleLabel, eleHandle)
        electrons = eleHandle.product()

#Initialize muon objects
        electron = 0
        antielectron = 0
#        print "Starting loop over particle list!"
#Begin GEN looping
#Loop through Particle list
        for ele in genParticles:
#            print "New particle!"
#            print "particle pdgID = ", mu.pdgId()
            #Check for muon
            if (ele.pdgId() == 11 and ele.status() == 1 and electron == 0 and (ele.mother().status() != 1)):#(abs(mu.mother().pdgId()) < 6) ):
                electron = getElectron(ele,args.debug)
                if electron!=0:
                    numElectrons+=1
                    totalElectrons+=1
                print "The number of electrons so far are: %d" %numElectrons
#Check for antimuon
            if (ele.pdgId() == -11 and ele.status() == 1 and antielectron == 0 and (ele.mother().status() != 1)):#(abs(mu.mother().pdgId()) < 6) ):

                antielectron = getElectron(ele, args.debug)
                if antielectron != 0:
                    numAntiElectrons+=1
                    totalAntiElectrons+=1
                print "The number of antielectrons so far are: %d" %numAntiElectrons
#        print "=================================================================================="



#Fill histograms
        numLoopTot+=1
        pairFound = printGenInfo(numLoopTot, numElectrons, numAntiElectrons)
        totalPairs += pairFound
        histEleMinusGenPt.Fill(electron.pt())
        histEleMinusGenEta.Fill(electron.eta())
        histEleMinusGenPhi.Fill(electron.phi())
        histElePlusGenPt.Fill(antielectron.pt())
        histElePlusGenEta.Fill(antielectron.eta())
        histElePlusGenPhi.Fill(antielectron.phi())
        histInvariantMass.Fill(math.sqrt(2*electron.pt()*antielectron.pt() * (math.cosh(electron.eta() - antielectron.eta()) - math.cos(electron.phi() - antielectron.phi() ) ) ) )

#Begin RECO looping
        for ele in electrons:
#            if mu.globalTrack().isNonnull():
#                if args.debug:
#                    print("Found muon with {} tracker layers with measurement".format(mu.globalTrack().hitPattern().trackerLayersWithMeasurement()))
#                    pass
#                histTkLayWMeas.Fill(mu.globalTrack().hitPattern().trackerLayersWithMeasurement())
#                histVdPxHits.Fill(mu.globalTrack().hitPattern().numberOfValidPixelHits())
#                histVdMuHits.Fill(mu.globalTrack().hitPattern().numberOfValidMuonHits())
#                histPtRes.Fill(mu.tunePMuonBestTrack().ptError()/mu.tunePMuonBestTrack().pt())
#                histDxy.Fill(mu.tunePMuonBestTrack().dxy())
#                histMatchMuStation.Fill(mu.numberOfMatchedStations())
            histPt.Fill(ele.pt())
#            print "The pt and track pt are: %f:%f" %(mu.pt(),mu.p4().pt())
            print "The kinematic properties are: %d:%f:%f:%f" %(ele.gsfTrack().charge(),ele.pt(),ele.superCluster().eta(),ele.superCluster().phi())
            histEta.Fill(ele.superCluster().eta())
            histPhi.Fill(ele.superCluster().phi())

print "The number of muons:antimuons are: %d:%d" %(numElectrons, numAntiElectrons)
print "The number times through the loop muon:antimuon: %d:%d" %(numloopElectrons, numloopAntiElectrons)
print "The total number of muons:antimuons:pairs are: %d:%d:%d" %(totalElectrons, totalAntiElectrons, totalPairs)
of.cd()
histInvariantMass.Write()
histEleMinusGenPt.Write()
histEleMinusGenEta.Write()
histEleMinusGenPhi.Write()
histElePlusGenPt.Write()
histElePlusGenEta.Write()
histElePlusGenPhi.Write()


#Runs loop to get cross section
t = f.Get("Runs")
xsHandle = Handle('GenRunInfoProduct')
xsLabel = "generator"

for i,run in enumerate(runs):
    runs.getByLabel(xsLabel, xsHandle)
    gen = xsHandle.product()

    print gen.internalXSec().value()
    pass

#Write to file and close; end program
of.cd()
#histTkLayWMeas.Write()
#histVdPxHits.Write()
#histVdMuHits.Write()
#histPtRes.Write()
#histDxy.Write()
#histMatchMuStation.Write()
histPt.Write()
histEta.Write()
histPhi.Write()
of.Close()
