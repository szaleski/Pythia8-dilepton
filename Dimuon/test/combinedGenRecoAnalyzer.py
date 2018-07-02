import ROOT as rt
import sys,os
import argparse
import math

from DataFormats.FWLite import Events, Handle, Runs

#function that makes dictionary of histograms to fill with various parameters
#place: 1 = outside loop, 2 = inside loop 
def makeDictionary(particle, place, objectP, isAntiParticle):
    hists = {}
    if particle == "electron":
        if not isAntiParticle:
            hists = {
                histEleMinusGenPt : objectP.pt,
                histEleMinusGenEta : objectP.eta,
                histEleMinusGenPhi : objectP.phi
                }
        else:
            hists = {
                histElePlusGenPt : objectP.pt, 
                histElePlusGenEta : objectP.eta, 
                histElePlusGenPhi : objectP.phi
                }
    if particle == "muon" and place == 1: 
        if not isAntiParticle:
            hists = {
                histMuMinusGenPt : objectP.pt,
                histMuMinusGenEta : objectP.eta,
                histMuMinusGenPhi : objectP.phi
                }
        else:
            hists = {
                histMuPlusGenPt : objectP.pt, 
                histMuPlusGenEta : objectP.eta, 
                histMuPlusGenPhi : objectP.phi
                } 
    if particle == "muon" and place == 2:
        hists = {
            histTkLayWMeas : objectP.globalTrack().hitPattern().trackerLayersWithMeasurement, 
            histVdPxHits : objectP.globalTrack().hitPattern().numberOfValidPixelHits, 
            histVdMuHits : objectP.globalTrack().hitPattern().numberOfValidMuonHits, 
            
            histDxy : objectP.tunePMuonBestTrack().dxy, 
            histMatchMuStation : objectP.numberOfMatchedStations,    
            histPt:objectP.tunePMuonBestTrack().pt, 
            histEta:objectP.tunePMuonBestTrack().eta,
            histPhi:objectP.tunePMuonBestTrack().phi
       }
    return hists

def findInvariantMass(antiPt, pt, antiEta, eta, antiPhi, phi):
    deltaEta = abs(eta - antiEta)
    deltaPhi = abs(phi - antiPhi)
    x = math.cosh(deltaEta) 
    y = math.cos(deltaPhi)
    InvariantMass = 2*pt*antiPt * abs(x-y) 
    if InvariantMass > 0:
        return math.sqrt(InvariantMass)
    return 0 

def changePt(pt, kappa):
    a = kappa*pt
    b = 1-a
    y = pt/b 
    return y 
def changePt1(pt, kappa):
    a = kappa*pt
    b = 1-a
    y = pt/b 
    return y
# Function to get (anti)particles that we are interested in from the MC bank
# returns 0 is particle has no mother
def getParticle(p, debug):
    particle = p
    if (debug):
        print p.pdgId()
        print p.status()
        print "Number of particle mothers: %d" %p.numberOfMothers()
        print "Particle mother is: %d \t with status: %d" %(p.mother().pdgId(), p.mother().status())
    nu = p
    motherId = nu.pdgId()
       # print(motherId)
    #print "nuId: %d" %nu.pdgId()
    while (abs(motherId) > 6): # not a quark
        if(nu.mother() and nu):
            nu = nu.mother()
            motherId = nu.pdgId()
            #print("motherId", motherId)
            if (args.debug):
                print "motherId: %d" %motherId
            if (abs(motherId) > 13): #not a particle
                particle = 0
            if (abs(motherId) < 7): #is a quark 
                particle = p
        else:
            print("error")
            return 0 
                #print(particle) 
    return particle


#Function to print Muon Event information
def printGenInfo(numLoopTot, numParticles, numAntiParticles):
    pairFound = 0
    #print "The number of events is: %d" %numLoopTot
    if (numLoopTot != numAntiParticles):
        print "Look HERE!"
        numLoopTot+=1
        #print "The number of events is: %d" %numLoopTot
    if (numParticles != numAntiParticles):
        print "number of Muons different from antimuons!!!"
    else:
        pairFound = 1
    return pairFound
    print"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print "                                          pt:eta:phi:id:status"
    print "The muon     information in this event is: %f:%f:%f:%d:%d" %(particle.pt(),particle.eta(),particle.ephi(),particle.pdgId(),particle.status())
    print "The antimuon information in this event is: %f:%f:%f:%d:%d" %(antiparticle.pt(),antiparticle.eta(),antiparticle.phi(),antiparticle.pdgId(),antiparticle.status())
    print"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#            print "The number of muons so far are: %d" %numMuons                    


# Begin main method

parser = argparse.ArgumentParser()
#parser.add_argument("infile", help="Make plots for cosmics",type=str)
#parser.add_argument("--local", help="local file",action="store_true")
parser.add_argument("-d", "--debug", help="debugging information",action="store_true")
parser.add_argument("-i", help="Input miniAOD File", type=str)
parser.add_argument("-p", help= "Input Particle", type=str) 

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

particle1 = args.p 
print(particle1) 

#Instantiate Histograms

histTkLayWMeas = rt.TH1D("histTkLayWMeas","TrackerLayersWithMeasurement",25,-0.5,24.5) #counts number of layers with hits
histVdPxHits = rt.TH1D("histVdPxHits", "ValidPixelHits",25,-0.5,24.5) #number of hits in the pixel layers that are valid
histVdMuHits = rt.TH1D("histVdMuHits", "ValidMuonHits", 25,-0.5,54.5);#number of valid hits in the muon chambers
histPtRes = rt.TH1D("histPtRes", "Ptresolution", 25, -0.15, 0.80); 
histDxy = rt.TH1D("histDxy", "Dxy", 25, -0.70, 0.70) 
histMatchMuStation = rt.TH1D("histMatchMuStation", "MatchedMuonStations", 20, -0.5, 19.5) #number of muon stations that match to the muon tracks 

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

histEleMinusGenPt = rt.TH1D("histEleMinusGenPt","Pt",100,0.,500.)
histEleMinusGenEta = rt.TH1D("histEleMinusGenEta", "Eta",50,-3.5,3.5)
histEleMinusGenPhi = rt.TH1D("histEleMinusGentPhi", "Phi", 50,-3.15,3.15);
histElePlusGenPt = rt.TH1D("histElePlusGenPt","Pt",100,0.,500.)
histElePlusGenEta = rt.TH1D("histElePlusGenEta", "Eta",50,-3.5,3.5)
histElePlusGenPhi = rt.TH1D("histElePlusGentPhi", "Phi", 50,-3.15,3.15);

histMass300 = rt.TH2D("histMass300", "histMass300", 100, 2000, 3000, 100, 0, 0.5) 
histMass = rt.TH1D("histMass", "Mass", 100, 0, 1000) 
histPerDiff = rt.TH1D("histPerDiff", "PercentDiff", 100, 0, 0.05) 
#array of hists to write for muons
histsMu = [histMuMinusGenPt, histMuMinusGenEta, histMuMinusGenPhi, histMuPlusGenPt, histMuPlusGenEta, histMuPlusGenPhi, histInvariantMass, histTkLayWMeas, histVdPxHits, histVdMuHits, histPtRes, histDxy, histMatchMuStation, histPt, histEta, histPhi] 
#array of hists to write for electrons
histsEle = [histEleMinusGenPt, histEleMinusGenEta, histEleMinusGenPhi, histElePlusGenPt, histElePlusGenEta, histElePlusGenPhi, histInvariantMass]

kappa = 0.00005


#Counters to track number of muons
totalPairs = 0
totalParticles = 0
totalAntiParticles = 0
count = 0
count1 =0 


#Loop over ROOT file list

for x in array:
    print "\n Now Reading file:"
    print x
    print "\n"
    f = rt.TNetXNGFile("%s/%s"%(eossrc,x),"read")
    events  = Events("%s/%s"%(eossrc,x))
    runs = Runs("%s/%s"%(eossrc, x))
                     
               
#Get Events Tree and create handles for GEN and RECO
    t = f.Get("Events")
    genHandle = Handle('std::vector<reco::GenParticle>')
    genLabel  = "prunedGenParticles"
    if particle1 == "muon":
        handle = Handle('std::vector<pat::Muon>')
        label = "slimmedMuons"
        #print(handle) 
    else:
        handle = Handle('std::vector<pat::Electron>')
        label = "slimmedElectrons"

    rt.TH1.SetDefaultSumw2()

#Counters to track number of muons in event and store particle object
    numloopParticles = 0
    numParticles = 0
    numloopAntiParticles = 0
    numAntiParticles = 0
    particle = 0
    antiparticle = 0
    numLoopTot = 0

#Event Loop
    for i,ev in enumerate(events):
        ev.getByLabel(genLabel, genHandle)
        genParticles = genHandle.product()
        ev.getByLabel(label, handle)
        particles = handle.product()
#Initialize muon objects
        particle  = 0
        antiparticle = 0
#Begin GEN looping
#Loop through Particle list
        for p in genParticles:
            #Check for muon 
            if ((p.pdgId() == 13 or p.pdgId() == 11)): 
                if (p.status() == 1 and particle == 0 and (p.mother().status() != 1)):
                    particle = getParticle(p, args.debug)
                    if particle!=0:
                        numParticles+=1
                        totalParticles+=1
                    #print "The number of muons so far are: %d" %numParticles

#Check for antimuon
            if ((p.pdgId() == -13 or p.pdgId() == -11)):
                if (p.status() == 1 and antiparticle == 0 and (p.mother().status() != 1)):
                    antiparticle = getParticle(p, args.debug)
                    if antiparticle != 0:
                        numAntiParticles+=1
                        totalAntiParticles+=1
                    #print "The number of antimuons so far are: %d" %numAntiParticles

                
#        print "=================================================================================="

#Fill histograms
        if particle!=0 and antiparticle!=0:
            numLoopTot+=1
            pairFound = printGenInfo(numLoopTot, numParticles, numAntiParticles)
            totalPairs += pairFound
        
#loop to fill histograms for electrons and positrons
        if particle1 == "electron":
            hists1 = makeDictionary(particle1, 1, particle, False)
            hists2 = makeDictionary(particle1, 1, antiparticle, True)
            for (key, val),(key2, val2) in zip(hists1.items(), hists2.items()):
                key.Fill(val()) 
                key2.Fill(val2()) 
            histInvariantMass.Fill(math.sqrt(2*particle.pt()*antiparticle.pt() * (math.cosh(particle.eta() - antiparticle.eta()) - math.cos(particle.phi() - antiparticle.phi() ) ) ) )
#loops to fill histograms for muons and antimuons 
        if particle1 == "muon":
            hists3 = makeDictionary(particle1, 1, particle, False) 
            hists4 = makeDictionary(particle1, 1, antiparticle, True)
            for (key, val),(key2, val2) in zip(hists3.items(), hists4.items()):
                key.Fill(val()) 
                key2.Fill(val2())
            mass1 = math.sqrt(2*particle.pt()*antiparticle.pt() * (math.cosh(particle.eta() - antiparticle.eta()) - math.cos(particle.phi() - antiparticle.phi() ) ) )
            print("starting mass", mass1)
        
            #print("starting mass2", mass2)
            
            antiPt = changePt(antiparticle.pt(), kappa)
            print("antipt", antiPt)
            
            pt = changePt(particle.pt(), kappa)
            print("pt", pt)
            mass = findInvariantMass(antiPt, pt, antiparticle.eta(), particle.eta(), antiparticle.phi(), particle.phi())
            print("mass", mass)
            #masses = [2100, 2200, 2300, 2400, 2500]
            #for m in masses: 
               # print(m)
            if mass1 > 2100:
                count = count + 1
                print("count", count) 
                
            if mass > 2100:
                #print(m)
                count1 = count1 + 1
                #print("mass", mass)   
                print("count1", count1) 
                if not count == 0:
                    percentDiff = (count1-count)
                    percentDiv = float(percentDiff)/count
                    print("percentDiff", percentDiv) 
                    histMass300.Fill(mass1, percentDiv)
                #print(percentDiff, mass1, mass) 
             
            histInvariantMass.Fill(math.sqrt(2*particle.pt()*antiparticle.pt() * (math.cosh(particle.eta() - antiparticle.eta()) - math.cos(particle.phi() - antiparticle.phi() ) ) ) )
                    #

                    #print("mass:", mass) 

        #if particle=0, then loop to next particle
        else:
            continue
#loop through muon hists
        if particle1 == "muon": 
            for p in particles:
                if p.globalTrack().isNonnull():
                    if args.debug:
                        print("Found muon with {} tracker layers with measurement".format(mu.globalTrack().hitPattern().trackerLayersWithMeasurement()))
                        pass
                    hists = makeDictionary(particle1, 2, p, "s")
                    histPtRes.Fill(p.tunePMuonBestTrack().ptError()/p.tunePMuonBestTrack().pt())
                    for key, val in hists.items():
                        key.Fill(val())
                    
#loop through electron hists
        else:
            for p in particles:         
                histPt.Fill(p.pt())
#            print "The pt and track pt are: %f:%f" %(mu.pt(),mu.p4().pt())
                print "The kinematic properties are: %d:%f:%f:%f" %(p.gsfTrack().charge(),p.pt(), p.superCluster().eta(), p.superCluster().phi())
                histEta.Fill(p.superCluster().eta())
                histPhi.Fill(p.superCluster().phi())
        

#print "The number of muons:antimuons are: %d:%d" %(numParticles, numAntiParticles)
print "The number times through the loop muon:antimuon: %d:%d" %(numloopParticles, numloopAntiParticles)
print "The total number of muons:antimuons:pairs are: %d:%d:%d" %(totalParticles, totalAntiParticles, totalPairs)
of.cd()



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
#write hists for muons
print("countmass1", count)
print("count mass", count1)
if particle1 == "muon":
    for i in histsMu:
        i.Write()
#write hists for electrons 
else:
    for i in histsEle:
        i.Write()
histMass300.Fit("expo") 
histMass300.Write()
#histMass300Down.DrawClone("Same") 
of.Close()


'''
                    
'''
