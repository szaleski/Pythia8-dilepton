from optparse import OptionParser
import ROOT as r


parser = OptionParser()
parser.add_option("-f", "--file1", dest='filename1',
	help="",type="string")

parser.add_option("-g", "--file2", dest='filename2',
	help="",type="string")
(options,args)=parser.parse_args()

filedir="/afs/cern.ch/user/s/szaleski/CMSSW_744_MCGen/src/GenStudy/Dimuon/test/"
file1 = r.TFile(filedir+options.filename1,"read")
file2 = r.TFile(filedir+options.filename2,"read")
print options.filename1,file1
print options.filename2,file2
r.gROOT.ProcessLine(".x CutComparisonNew.cc++(\"%s\",\"%s\")"%(options.filename1,options.filename2))
raw_input("press enter to quit")
