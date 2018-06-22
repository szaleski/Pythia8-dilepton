import ROOT as rt
import sys,os
import argparse
import math

from DataFormats.FWLite import Events, Handle, Runs

# Function to check for (anti)muons that we are interested in
def getMuon(mu, debug):
    muon = mu
    if (debug):
        print mu.pdgId()
        print mu.status()
        print "Number of muon mothers: %d" %mu.numberOfMothers()
        #    numloopMuons+=1
        print "Muon mother is: %d \t with status: %d" %(mu.mother().pdgId(), mu.mother().status())
    nu = mu

    motherId = nu.pdgId()
    print "nuId: %d" %nu.pdgId()
    while (abs(motherId) > 6):
        nu = nu.mother()
        motherId = nu.pdgId()
        if (args.debug):
            print "motherId: %d" %motherId
        if (abs(motherId) > 13):
            muon = 0
        if (abs(motherId) < 7):
            muon = mu
    return muon

#Function to print Muon Event information
def printGenInfo(numLoopTot, numMuons, numAntiMuons):
    pairFound = 0
    print "The number of events is: %d" %numLoopTot
    if (numLoopTot != numAntiMuons):
        print "Look HERE!"
        numLoopTot+=1
        print "The number of events is: %d" %numLoopTot
    if (numMuons != numAntiMuons):
        print "number of Muons different from antimuons!!!"
    else:
        pairFound = 1
    return pairFound
    print"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print "                                          pt:eta:phi:id:status"
    print "The muon     information in this event is: %f:%f:%f:%d:%d" %(muon.pt(),muon.eta(),muon.phi(),muon.pdgId(),muon.status())
    print "The antimuon information in this event is: %f:%f:%f:%d:%d" %(antimuon.pt(),antimuon.eta(),antimuon.phi(),antimuon.pdgId(),antimuon.status())
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

histTkLayWMeas = rt.TH1D("histTkLayWMeas","TrackerLayersWithMeasurement",25,-0.5,24.5)
histVdPxHits = rt.TH1D("histVdPxHits", "ValidPixelHits",25,-0.5,24.5)
histVdMuHits = rt.TH1D("histVdMuHits", "ValidMuonHits", 25,-0.5,54.5);
histPtRes = rt.TH1D("histPtRes", "Ptresolution", 25, -0.15, 0.80); 
histDxy = rt.TH1D("histDxy", "Dxy", 25, -0.70, 0.70)
histMatchMuStation = rt.TH1D("histMatchMuStation", "MatchedMuonStations", 20, -0.5, 19.5)
histPt =rt.TH1D("histPt", "Pt", 200, 0.0, 2000)
histEta = rt.TH1D("histEta", "#Eta", 50, -4.0, 4.0)
histPhi = rt.TH1D("histPhi", "#Phi", 30, -3.2, 3.2)


histMuMinusGenPt = rt.TH1D("histMuMinusGenPt","Pt",100,0.,500.)
histMuMinusGenEta = rt.TH1D("histMuMinusGenEta", "Eta",50,-3.5,3.5)
histMuMinusGenPhi = rt.TH1D("histMuMinusGentPhi", "Phi", 50,-3.15,3.15);
histMuPlusGenPt = rt.TH1D("histMuPlusGenPt","Pt",100,0.,500.)
histMuPlusGenEta = rt.TH1D("histMuPlusGenEta", "Eta",50,-3.5,3.5)
histMuPlusGenPhi = rt.TH1D("histMuPlusGentPhi", "Phi", 50,-3.15,3.15);
histInvariantMass = rt.TH1D("histInvariantMass","Mass Dist",300,0.,3000)

#Counters to track number of muons

totalPairs = 0
totalMuons = 0
totalAntiMuons = 0

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
    muHandle = Handle('std::vector<pat::Muon>')
    muLabel  = "slimmedMuons"


    rt.TH1.SetDefaultSumw2()

#Counters to track number of muons in event and store particle object
    numloopMuons = 0
    numMuons = 0
    numloopAntiMuons = 0
    numAntiMuons = 0
    muon = 0
    antimuon = 0
#tmpMuCount = 0
#tmpAmuCount = 0
    numLoopTot = 0
#Event Loop
    for i,ev in enumerate(events):
        ev.getByLabel(genLabel, genHandle)
        genMuons = genHandle.product()
        ev.getByLabel(muLabel, muHandle)
        muons = muHandle.product()

#Initialize muon objects
        muon = 0
        antimuon = 0
#        print "Starting loop over particle list!"
#Begin GEN looping
#Loop through Particle list
        for mu in genMuons:
#            print "New particle!"
#            print "particle pdgID = ", mu.pdgId()
            #Check for muon
            if (mu.pdgId() == 13 and mu.status() == 1 and muon == 0 and (mu.mother().status() != 1)):#(abs(mu.mother().pdgId()) < 6) ):
                muon = getMuon(mu,args.debug)
                if muon!=0:
                    numMuons+=1
                    totalMuons+=1
                print "The number of muons so far are: %d" %numMuons
#Check for antimuon
            if (mu.pdgId() == -13 and mu.status() == 1 and antimuon == 0 and (mu.mother().status() != 1)):#(abs(mu.mother().pdgId()) < 6) ):

                antimuon = getMuon(mu, args.debug)
                if antimuon != 0:
                    numAntiMuons+=1
                    totalAntiMuons+=1
                print "The number of antimuons so far are: %d" %numAntiMuons
#        print "=================================================================================="



#Fill histograms
        numLoopTot+=1
        pairFound = printGenInfo(numLoopTot, numMuons, numAntiMuons)
        totalPairs += pairFound
        histMuMinusGenPt.Fill(muon.pt())
        histMuMinusGenEta.Fill(muon.eta())
        histMuMinusGenPhi.Fill(muon.phi())
        histMuPlusGenPt.Fill(antimuon.pt())
        histMuPlusGenEta.Fill(antimuon.eta())
        histMuPlusGenPhi.Fill(antimuon.phi())
        histInvariantMass.Fill(math.sqrt(2*muon.pt()*antimuon.pt() * (math.cosh(muon.eta() - antimuon.eta()) - math.cos(muon.phi() - antimuon.phi() ) ) ) )

#Begin RECO looping
        for mu in muons:
            if mu.globalTrack().isNonnull():
                if args.debug:
                    print("Found muon with {} tracker layers with measurement".format(mu.globalTrack().hitPattern().trackerLayersWithMeasurement()))
                    pass
                histTkLayWMeas.Fill(mu.globalTrack().hitPattern().trackerLayersWithMeasurement())
                histVdPxHits.Fill(mu.globalTrack().hitPattern().numberOfValidPixelHits())
                histVdMuHits.Fill(mu.globalTrack().hitPattern().numberOfValidMuonHits())
                histPtRes.Fill(mu.tunePMuonBestTrack().ptError()/mu.tunePMuonBestTrack().pt())
                histDxy.Fill(mu.tunePMuonBestTrack().dxy())
                histMatchMuStation.Fill(mu.numberOfMatchedStations())
                histPt.Fill(mu.tunePMuonBestTrack().pt())
                histEta.Fill(mu.tunePMuonBestTrack().eta())
                histPhi.Fill(mu.tunePMuonBestTrack().phi())

print "The number of muons:antimuons are: %d:%d" %(numMuons, numAntiMuons)
print "The number times through the loop muon:antimuon: %d:%d" %(numloopMuons, numloopAntiMuons)
print "The total number of muons:antimuons:pairs are: %d:%d:%d" %(totalMuons, totalAntiMuons, totalPairs)
of.cd()
histInvariantMass.Write()
histMuMinusGenPt.Write()
histMuMinusGenEta.Write()
histMuMinusGenPhi.Write()
histMuPlusGenPt.Write()
histMuPlusGenEta.Write()
histMuPlusGenPhi.Write()


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
histTkLayWMeas.Write()
histVdPxHits.Write()
histVdMuHits.Write()
histPtRes.Write()
histDxy.Write()
histMatchMuStation.Write()
histPt.Write()
histEta.Write()
histPhi.Write()
of.Close()
