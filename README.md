# fly_project nanoaodrdframework 

Analysis of CMS experiment Nano AOD format using ROOT RDataFrame

The package shows how to analyze using the data frame concept.
In a data frame concept, the operations on data should be defined
in terms of functions, which means that result of applying
a function should produce outputs and no internal states are allowed.
In ROOT, RDataFrame allows us to borrow the concept for analysis.

This package code implements a hadronic channel analysis using NanoAOD format.
However, it can be adapted for any flat ROOT trees.

## I. Introduction

- The advantage of using RDataFrame is that you don't have to
worry about looping, filling of histograms at the right place,
writing out trees. In some sense you can concentrate more on
thinking about data analysis.

- Draw back is that you have to think of a new way of doing old things.
This could getting used to. Alos complicated algorithms may not
be so easy to implement with RDataFrame

- Purely data frame concept is not ideal for HEP data since
we have to treat real data and simulated data differently.
We must apply different operations depending on the input data,
   such as scale factors and corrections.
Therefore, some form of global state should be stored. 
In this package, the data frame concept is used together
with object oriented concept for this purpose.


## II. Code

- The code consists of a main class NanoAODAnalyzerrdframe. 
There is one.h header file and one .cpp source file for it.
The class has several methods:
    - Object definitions (selectElectrons, selectMuons, ...)
    - Additional derived variables definition (defineMoreVars)
    - Histogram definitions (bookHists)
    - Selections (defineCuts)
    - Plus some utility methods/functions (generate_4vec,calculateDR,reconstruction,..)

- Users should modify: object definitions, define additional variables, histogram definitions, selections.

## III. Clone Repository 

if using standalone:
* ROOT >= 6.24
* correctionlib >=2.0.0 (https://github.com/cms-nanoAOD/correctionlib)
* rapidjson (https://github.com/Tencent/rapidjson)
* nlohman/json (https://github.com/nlohmann/json)

if CMSSW is available: the package doesn't depend on any CMSSW libraries, but some of the dependent libraries are available, so this makes it easier to compile.

* Setup CMSSW

      > cmsrel CMSSW_12_3_X

      > cd CMSSW_12_3_X/src

CMSSW 12.3 has ROOT 6.22 but there are run time errors when using this. So you must setup ROOT 6.24.

      > source /cvmfs/cms.cern.ch/slc7_amd64_gcc900/lcg/root/6.24.07-db9b7135424812ed5b8723a7a77f4016/bin/thisroot.(c)sh  

  * Clone with SSH :
      > git clone git@github.com:cdozen/test_nanofly.git

  * Clone with HTTPS:
      > git clone https://github.com/cdozen/test_nanofly.git
  

## IV. Compiling

        > make -j(24)

  this will compile and create a libnanoaodrdframe.so shared library that can be loaded in a ROOT session or macro:  gSystem->Load("libnanoadrdframe.so");


## V. Running over large dataset

      > ./submitanalysisjob.py jobconfiganalysis.py

  * "submitanalysisjob.py" script loop over the directories and submit jobs over the configuation python module (e.g: processnanoad.py)
      > os.system('./processnanoaod.py '+ indir + ' ' + outdir + ' ' + jobconfmod  + ' > ' + outfile + ' 2>&1 &')

  * "processnanoaod.py" script can automatically run over all ROOT files in an input directory and run over the subclass (e.g: BaseAnalyser.cpp)
      > aproc = ROOT.BaseAnalyser(t, outputroot)

  * "jobconfiganalysis.py" script contatins all the options are to be set. It contains 3 python directories:
    * one for configurations for datatype, data year, runtype, various POG corrections, golden JSON, input and output root tree names, etc. ( The correction is an optional argument. As a default no corrections are applied on this step. User should modify processnanoaod.py if one wishes to do so (search for skipcorrections in that file).)
       
            # options for Nanoaodrdframe
            config = {
            'intreename': "Events",       # tree name of input file(s)
            'outtreename': "outputTree2", # tree name of output file(s) it cannot be the same as thinput tree name or it'll crash
            'year': 2018,                 #data year (2016,2017,2018)
            'runtype': 'UL',              # is ReReco or Ultra Legacy
            'datatype': -1,               # -1: checking "gen" branch if it's exist run over MC. 1: Data0:MC
            }

    * The second dictionary contains processing options, to split processing into multiple jobs, produce one output file per input file, skip files already processed, whether to dive into subdirectories recursively and process the files there. 
       
            #processing options
            procflags = {
            'split': 1,                   # how many jobs?
            'allinone': True,             # if "False" one output file per input file, if True then onoutput file for everything
            'skipold': True,              # if "True" then skip existing analyzed files
            'recursive': True,            # travel through the subdirectories and their subdirecties wheprocessing.
                                      # be careful not to mix   MC and real DATA in them.
            'saveallbranches': False,     # if False then only selected branches which is done in thecpp file will be saved
            'nrootfiles': 1,              #How many input files?
            }
    
    * The third dictionary contains input directory where your input files are, output directory, and text output file for any print out or error messages.

            #input directory where your input root tree resides && output directory wheere your output should go && dump of stderr/stdout to file
            nanoaod_inputdir_outputdir_pairs = [
                    ['testinputdata/MC/2018','analyzed/test_nano_v1.root', 'stderrout.out' ],
            ]

Ref.:https://github.com/suyong-choi/nanoaodrdframe
