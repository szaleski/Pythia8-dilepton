import ROOT as rt
import sys,os
import argparse
import math

from DataFormats.FWLite import Events, Handle, Runs

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Make plots for cosmics",type=str)
parser.add_argument("--local", help="local file",action="store_true")
parser.add_argument("-d", "--debug", help="debugging information",action="store_true")

args = parser.parse_args()

rt.gROOT.SetBatch(True)
# load FWLite C++ libraries
rt.gSystem.Load("libFWCoreFWLite.so");
rt.gSystem.Load("libDataFormatsFWLite.so");
rt.FWLiteEnabler.enable()
if args.local:
    f = rt.TFile("%s"%(args.infile),"read")
    events  = Events("%s"%(args.infile))
    runs = Runs("%s"%(args.infile))
else:
    eossrc = "root://eoscms.cern.ch//eos/cms"
    eossrc = "root://cmseos.fnal.gov//"
    eossrc = "root://cmsxrootd.fnal.gov//"
    f = rt.TNetXNGFile("%s/%s"%(eossrc,args.infile),"read")
    events  = Events("%s/%s"%(eossrc,args.infile))
    runs = Runs("%s/%s"%(eossrc,args.infile))
    pass

t = f.Get("Events")
genHandle = Handle('std::vector<reco::GenParticle>')
genLabel  = "prunedGenParticles"



rt.TH1.SetDefaultSumw2()
of = rt.TFile("genOutput.root","recreate")

histMuMinusGenPt = rt.TH1D("histMuMinusGenPt","Pt",100,0.,500.)
histMuMinusGenEta = rt.TH1D("histMuMinusGenEta", "Eta",50,-3.5,3.5)
histMuMinusGenPhi = rt.TH1D("histMuMinusGentPhi", "Phi", 50,-3.15,3.15);
histMuPlusGenPt = rt.TH1D("histMuPlusGenPt","Pt",100,0.,500.)
histMuPlusGenEta = rt.TH1D("histMuPlusGenEta", "Eta",50,-3.5,3.5)
histMuPlusGenPhi = rt.TH1D("histMuPlusGentPhi", "Phi", 50,-3.15,3.15);
histInvariantMass = rt.TH1D("histInvariantMass","Mass Dist",300,0.,3000)

numloopMuons = 0
numMuons = 0
numloopAntiMuons = 0
numAntiMuons = 0
muon = 0
antimuon = 0
for i,ev in enumerate(events):
    ev.getByLabel(genLabel, genHandle)
    genMuons = genHandle.product()
#    if args.debug:
#        print("Event {0:8d} had {1:3d} muons".format(i,len(genMuons)))
#        pass
#    print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    for mu in genMuons:
#        print mu.pdgId()
        if (mu.pdgId() == 13 and mu.status() == 1 and (mu.mother().status() != 1)):#(abs(mu.mother().pdgId()) < 6) ):
            print mu.pdgId()
            print mu.status()
            print "Number of muon mothers: %d" %mu.numberOfMothers()
            numloopMuons+=1
            print "Muon mother is: %d \t with status: %d" %(mu.mother().pdgId(), mu.mother().status())
            nu = mu

            motherId = nu.pdgId()
            print "nuId: %d" %nu.pdgId()
            while (abs(motherId) > 6):
                nu = nu.mother()
                motherId = nu.pdgId()
                print "motherId: %d" %motherId
                if (abs(motherId) > 13):
                    break
                if (abs(motherId) < 7):
                    muon = mu
                    histMuMinusGenPt.Fill(mu.pt())
                    histMuMinusGenEta.Fill(mu.eta())
                    histMuMinusGenPhi.Fill(mu.phi())
                    numMuons+=1
                    
#            print mu.pt()
#            print mu.eta()
#            print mu.phi()
#            print "The number of muons so far are: %d" %numMuons
        if (mu.pdgId() == -13 and mu.status() == 1 and (mu.mother().status() != 1)):#(abs(mu.mother().pdgId()) < 6) ):
            print mu.pdgId()
            print mu.status()
            numloopAntiMuons+=1
            print "Antimuon mother is: %d \t with status: %d" %(mu.mother().pdgId(), mu.mother().status())
            nu = mu

            motherId = nu.pdgId()
            print "nuId: %d" %nu.pdgId()
            while (abs(motherId) > 6):
                nu = nu.mother()
                motherId = nu.pdgId()
                print "motherId: %d" %motherId
                if (abs(motherId) > 13):
                    break
                if (abs(motherId) < 7):
                    numAntiMuons+=1
                    antimuon = mu
                    histMuPlusGenPt.Fill(mu.pt())
                    histMuPlusGenEta.Fill(mu.eta())
                    histMuPlusGenPhi.Fill(mu.phi())
            print "The number of antimuons so far are: %d" %numAntiMuons
#        print "=================================================================================="





    print"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print "                                          pt:eta:phi:id:status"
    print "The muon     information in this event is: %f:%f:%f:%d:%d" %(muon.pt(),muon.eta(),muon.phi(),muon.pdgId(),muon.status())
    print "The antimuon information in this event is: %f:%f:%f:%d:%d" %(antimuon.pt(),antimuon.eta(),antimuon.phi(),antimuon.pdgId(),antimuon.status())
    print"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    histInvariantMass.Fill(math.sqrt(2*muon.pt()*antimuon.pt() * (math.cosh(muon.eta() - antimuon.eta()) - math.cos(muon.phi() - antimuon.phi() ) ) ) )
print "The number of muons:antimuons are: %d:%d" %(numMuons, numAntiMuons)
print "The number times through the loop muon:antimuon: %d:%d" %(numloopMuons, numloopAntiMuons)
of.cd()
histInvariantMass.Write()
histMuMinusGenPt.Write()
histMuMinusGenEta.Write()
histMuMinusGenPhi.Write()
histMuPlusGenPt.Write()
histMuPlusGenEta.Write()
histMuPlusGenPhi.Write()
of.Close()

#t = f.Get("Events")
#muHandle = Handle('std::vector<pat::Muon>')
#muLabel  = "slimmedMuons"

#histTkLayWMeas = rt.TH1D("histTkLayWMeas","TrackerLayersWithMeasurement",25,-0.5,24.5)
#histVdPxHits = rt.TH1D("histVdPxHits", "ValidPixelHits",25,-0.5,24.5)
#histVdMuHits = rt.TH1D("histVdMuHits", "ValidMuonHits", 25,-0.5,54.5);
#histPtRes = rt.TH1D("histPtRes", "Ptresolution", 25, -0.15, 0.80); 
#histDxy = rt.TH1D("histDxy", "Dxy", 25, -0.70, 0.70)
#histMatchMuStation = rt.TH1D("histMatchMuStation", "MatchedMuonStations", 20, -0.5, 19.5)
#histPt =rt.TH1D("histPt", "Pt", 200, 0.0, 2000)
#histEta = rt.TH1D("histEta", "#Eta", 50, -4.0, 4.0)
#histPhi = rt.TH1D("histPhi", "#Phi", 30, -3.2, 3.2)



#for i,ev in enumerate(events):
#    ev.getByLabel(muLabel, muHandle)
#    muons = muHandle.product()
#    if args.debug:
#        print("Event {0:8d} had {1:3d} muons".format(i,len(muons)))
#        pass
#    for mu in muons:
#        if mu.globalTrack().isNonnull():
#            if args.debug:
#                print("Found muon with {} tracker layers with measurement".format(mu.globalTrack().hitPattern().trackerLayersWithMeasurement()))
#                pass
#            histTkLayWMeas.Fill(mu.globalTrack().hitPattern().trackerLayersWithMeasurement())
#            histVdPxHits.Fill(mu.globalTrack().hitPattern().numberOfValidPixelHits())
#            histVdMuHits.Fill(mu.globalTrack().hitPattern().numberOfValidMuonHits())
#            histPtRes.Fill(mu.tunePMuonBestTrack().ptError()/mu.tunePMuonBestTrack().pt())
#            histDxy.Fill(mu.tunePMuonBestTrack().dxy())
#            histMatchMuStation.Fill(mu.numberOfMatchedStations())
#            histPt.Fill(mu.tunePMuonBestTrack().pt())
#            histEta.Fill(mu.tunePMuonBestTrack().eta())
#            histPhi.Fill(mu.tunePMuonBestTrack().phi())
#        else:
#            if args.debug:
#                print("Muon had no valid global track reference")
#                pass
#            pass
#        pass
#    if args.debug and (i > 10): break
#    pass
#t = f.Get("Runs")
#xsHandle = Handle('GenRunInfoProduct')
#xsLabel = "generator"

#for i,run in enumerate(runs):
#    runs.getByLabel(xsLabel, xsHandle)
#    gen = xsHandle.product()

#    print gen.internalXSec().value()
#    pass


#of.cd()
#histTkLayWMeas.Write()
#histVdPxHits.Write()
#histVdMuHits.Write()
#histPtRes.Write()
#histDxy.Write()
#histMatchMuStation.Write()
#histPt.Write()
#histEta.Write()
#histPhi.Write()
#of.Close()
