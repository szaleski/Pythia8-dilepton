import ROOT as rt
import sys,os
import argparse

from DataFormats.FWLite import Events, Handle, Runs

parser = argparse.ArgumentParser()
#parser.add_argument("infile", help="Make plots for cosmics",type=str)
#parser.add_argument("--local", help="local file",action="store_true")
parser.add_argument("-i", help="Input miniAOD File", type=str)
parser.add_argument("-d", "--debug", help="debugging information",action="store_true")

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
#lines = [line.rstrip('\n') for line in open(fileTest.out)]
of = rt.TFile("output.root","recreate")

histTkLayWMeas = rt.TH1D("histTkLayWMeas","TrackerLayersWithMeasurement",25,-0.5,24.5)
histVdPxHits = rt.TH1D("histVdPxHits", "ValidPixelHits",25,-0.5,24.5)
histVdMuHits = rt.TH1D("histVdMuHits", "ValidMuonHits", 25,-0.5,54.5);
histPtRes = rt.TH1D("histPtRes", "Ptresolution", 25, -0.15, 0.80); 
histDxy = rt.TH1D("histDxy", "Dxy", 25, -0.70, 0.70)
histMatchMuStation = rt.TH1D("histMatchMuStation", "MatchedMuonStations", 20, -0.5, 19.5)
histPt =rt.TH1D("histPt", "Pt", 200, 0.0, 2000)
histEta = rt.TH1D("histEta", "#Eta", 50, -4.0, 4.0)
histPhi = rt.TH1D("histPhi", "#Phi", 30, -3.2, 3.2)

for x in array:
    print "\n Now Reading file:"
    print x
    print "\n"
    f = rt.TNetXNGFile("%s/%s"%(eossrc,x),"read")
    events  = Events("%s/%s"%(eossrc,x))
    runs = Runs("%s/%s"%(eossrc, x))
                     
                     
    t = f.Get("Events")
    muHandle = Handle('std::vector<pat::Muon>')
    muLabel  = "slimmedMuons"



    rt.TH1.SetDefaultSumw2()
    for i,ev in enumerate(events):
        ev.getByLabel(muLabel, muHandle)
        print "Getting Event Label"
        muons = muHandle.product()
        if args.debug:
            print("Event {0:8d} had {1:3d} muons".format(i,len(muons)))
            pass
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
            else:
                if args.debug:
                    print("Muon had no valid global track reference")
                    pass
                pass
            pass
        if args.debug and (i > 10): break
        pass
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


#print array[2]
