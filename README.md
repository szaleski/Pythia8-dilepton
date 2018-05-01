# Pythia8-dilepton

This packages is designed to perform GEN-only private MC productions using the CMSSW Pythia8 release that is included in the working CMSSW release. It also includes basic TTree-making files and simplistic histogramming that is intended primarly as a starting point to springboard to do whatever specific plotting is desired.

# Soon-to-be Outdated CMSSW_8_0_21 Instuctions (Valid for 2016 running)

First create a clean CMSSW_8_0_21 release area:
```
cmsrel CMSSW_8_0_21
cd CMSSW_8_0_21/src
cmsenv
```

Next clone the Pythia8 repository here

```
git clone -b cmssw8021 https://github.com/szaleski/Pythia8-dilepton.git GenStudy
cd GenStudy
```

Next use scram to build and compile the project
```
scram b -j 8

```

If there are compiler errors, rsolve these (should be minor, if any).
Upon successful build change to test directory and try running the genrator executable using cmsRun.
This creates a lot screen output and can be time consuming to run so we redirect the output into a .txt file using '>' and run both functions in the background using '&'
It does take command line input. Every option has a default value except for maxEvents...which defaults to infity!
YOU MUST ASSIGN A VALUE TO THIS VARIABLE EVERY TIME!!!
If you forget, simply kill the process with "Ctrl+C"

```
cd Dimuon/test

cmsRun ciCrossSecCalc_edm_cfg.py maxEvents=50000 minMass=300 maxMass=800 Lambda=10000 helicityLL=-1 ciGen=1 ULE=off pdgId=11 &> outFile.txt &
```

This creates an output ROOT file that is labeled by most of the input parameters. This file has all of the RAW event MC info. It is meant to be analysed using Dilepton_cfg.py
which is our TTree maker (or EDAnalyzer in CMS language). It takes an input .root file from GEN production with the .root extension and also takes in a user created output file name WITHOUT the .root extension.
This also displays alot of screen output so we redirect and run in the background again.

```
cmsRun Dilepton_cfg.py inFile=<name of file created by generator> filename=<some output filename> &> AnalyzerOut.txt &
```
 
The output file has branches and histograms that can be used to quickly check to see if the values make sense.

These can be opened in a TBrowser in ROOT by doing the following:
```
root -l <some output filename>.root
[0]TBrowser b

```

If everything looks OK, it is likely that you would like to plot these histograms. 
We offer a basic histogramming and plotting script, makeHistograms.cc that provides simple access to the TTree variables and some basic ROOT plotting functionality.
It takes three input variables, inputfile name, particle type, and output file name.
To run this do


```
#for BASH shell
root -x -b -q makeHistograms.cc++(\"outputFile.root\",\"muon\",\"histogramOut\")

#for tcsh shell
root -x -b -q makeHistograms.cc++\(\"outputFile.root\",\"electron\",\"histogramOut\"\)

```

Once this has run, it will have made an output ROOT file that has stored the histograms that can be either:
a) modified in the TBrowser
b) saved as a .C ROOT macro that can be modified that way

The .root file will also have all TCanvas objects saved inside it.
The TCanvas objects also get written out to .png file...if a different format is desired, 
this can easily be modified inside makeHistogram.cc in the SaveAs() function by changing  .png to .pdf, or some other desired extension.

# CRAB submission steps for Workflow

When running local simulations do the following:

GENSIM:

in test/genSimCrabConfig folder:

One can run locally doing the mc16GenSim_cfg.py, using the same command line input as with ciCrossSecCalc_edm_cfg.py. E.g.

```
cmsRun mc16GenSim_cfg.py maxEvents=10000 ciGen=1 helicityLL=-1 minMass=300 Lambda=10000 >& outputFile.log &

```
To submit via CRAB, the crabConfig file needs to be modied. It is highly recomended that the CRAB configuration file name has it's own unique name that clearly identifies the sample that you want to submit. It is also strongly recommended that all naming information used within the configuration file is unique to the sample, especially the request name. 

Make the following changes:
config.General.requestName = '<some unique name that identifies the sample>'
config.JopType.psetName = 'mc16GenSim_cfg.py'
config.JobType.pyCfgParams =[<list of command line input that will be used for the sample just as one would for ciCrossSecCalc_edm_cfg.py>]

config.outputPrimaryDataSet = '<PrimaryDataset name for DAS that will identify the lepton flavor and Lambda value (N.B. all daughter samples e.g. RECO, AOD, miniAOD will inherit from this name)>'
config.Data.outputDatasetTag = '<secondary DAS name that will be used to uniquely identify a sample within the primary dataset folder>'

For first time setup only you will need to create the crab_projects folder:
```
mkdir crab_projects
```

This is the directory that CRAB will store local CRAB output into.

To submit to CRAB do the following (once only per cd):
```
cmsenv
voms-proxy-init -voms cms --valid 168:00
source cvmfs/cms.cern.ch/crab3/crab.sh
```
Then for each CRAB job that you would like to submit do the following:
```
crab submit -c <crabConfigFile>
```
You can check the status of the jobs using:
```
crab status --dir crab_projects/<CRAB task request name>
```

If any jobs need to be resubmitted, do the following (N.B. only jobs in FAILED state can be resubmitted):

```
crab resubmit --dir crab_projects/<CRAB task request name>  --jobids=<comma separated list of jobs e.g. 1,6,10-70>
```

Jobs can also be followed using the task monitoring dashboard.

Once jobs are completed, the RECO process can be started.


RECO:
To start change to the RECO CRAB directory, test/recoCrabConfig.
To run Locally, make sure that the input file path in the file matches the location of the GENSIM output ROOT file. Then simply do:
```
cmsRun mc16RECO_cfg.py
```

To submit to CRAB, the process is similar to submitting the GENSIM step. The primary difference comes from modifying the CRAB config file. To do this, open the CRAB config file in a text editor.
The following will need changes:
config.General.requestName = '<CRAB task name that identifies the unique sample>'
config.Data.inputDatset =' <DAS datset name for output GENSIM sample>'
To find the ouput dataset name you will need to query DAS using the following:

"dataset = /<GENSIM outputPrimaryDataset variable from GENSIM CRABCONFIG>/<username>-<outputDatasetTag variable from GENSIM CRABCONFIG>-*/USER" using the prod/phys03 Database (DBS). Simply copy and paste the appropriate sample name from the DAS page and paste it into the inputDataset variable in the RECO CRAB config file.

Lastly:
config.Data.outputDatasetTag = '<output Tag that uniquely identifies the RECO (DIGIRAW) sample>'

Once these changes are made, CRAB submission is the same as GENSIM, simply use the RECO crabconfig file instead of the GENSIM in the command.


AOD:


For AOD, change to the test/aodCrabConfig directory and make the changes to the same fields specified in the RECO step in the crab config file, but update to use RECO as input dataset and change the request name and output tag to reflect the AOD.




miniAOD:

For miniAOD, change to the test/miniAODCrabConfig directory and make the changes to the same fields specified in the RECO step in the crab config file, but update to use the AOD as input dataset and change the request name and output tag to reflect miniAOD.


# Outdated CMSSW_7_4_4 Instructions (Valid for 2015 running)
In order to check this out you should have checked out a CMSSW release.
It has been tested on CMSSW_7_4_4.
To checkout package, please do the following in the src directory:
``` 
git clone https://github.com/szaleski/Pythia8-dilepton.git GenStudy
```
or
    
```
git clone git@github.com:szaleski/Pythia8-dilepton.git GenStudy
```

After the package has been checked out you can run the Generator  dimuon_Pythia8_gen.py that resides in the test directory.
Running this generator creates a .root output file. Also, it is possible to specify optional parameters at the command line, specifically, the number of events to produce (maxEvents), minimum mass cut (minMass), maximum mass cut (maxMass), model type (model), and output file name (outName). Any of the specified parameters may be omitted from the command, and will revert to their default values found in the mcCommandLineOptions_cfi.py file in the python dirctory. These can be hardcoded into the generator file though if one wishes.

To run please do the following:
1) change to test directory

```   
cmsRun dimuon_Pythia8_gen.py maxEvents=<number of events> minMass=<low mass cut> maxMass=<high mass cut> model=<model name> outName="<file name>"
```
This creates a .root output file that can be analyzed using Dimuon_cfg.py. This cfg file analyzes the TTree generated by dimuon_Pythia8_gen.py.

To run this file please do the following:
```
cmsRun Dimuon_cfg.py
```

It creates histograms for the following quantities of Boson, and dimuon pairs: pT, mass, energy, pseudorapidity, phi, theta, and charge. It also compares values of theta and phi for the muon pairs.

Do not use Plots.cc, this was an earlier attempt of Plots2.cc
Plots2.cc can be used to generate different canvases to be analyzed in the .root output file created by Dimuon_cfg.py
