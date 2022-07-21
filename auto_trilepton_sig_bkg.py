#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv, exit
from os import getcwd, chdir, path, system, popen
from subprocess import Popen
import numpy as np
from matplotlib import pyplot as plt
from termcolor import colored
import re
import numpy as np
import json

################################################################################

# Move up one directory
def gotoParentDir():
    chdir("..")

# Run a given terminal command
def run(cmd):
    print(colored(cmd, "cyan"))
    system(cmd)

# Set current working directory to this script's directory
def setWorkingDir():
    scriptPath = path.realpath(__file__)
    scriptDir = path.dirname(scriptPath)
    chdir(scriptDir)
    return

# Compile and run the analysis program
def compileRun():
    run("make -j24 > make.log 2>&1") # Redirect stdout & stderr to make.log
    run("./submitanalysisjob.py jobconfiganalysis.py")
    run("cat stderrout.log")
    return

def fetch(directory, command="ls"):
    chdir(directory)
    filename = popen(command).read().strip()
    gotoParentDir()
    return filename

# List files in a given directory in a condensed format
def ls(directory):
    command = f"ls {directory}"
    util = lambda cmd: popen(cmd).read().strip().replace("\n", " ") or ""
    print(colored(command, "cyan"), util(command))
    return

################################################################################

if __name__ == "__main__":
    files = {
        "sig": ["TT2L2Nu_2_Muons_signal.root"],
               #"TT2L2Nu_1_Muon_1_Electron_signal.root",
               #"TT2L2Nu_2_Electrons_signal.root"],
        "bkg": {
            "WZ":    ["TT2L2Nu_2_Muons_WZ.root"],
                     #"TT2L2Nu_1_Muon_1_Electron_WZ.root",
                     #"TT2L2Nu_2_Electrons_WZ.root"],
            "ttW":   ["TT2L2Nu_2_Muons_ttW.root"],
                     #"TT2L2Nu_1_Muon_1_Electron_ttW.root",
                     #"TT2L2Nu_2_Electrons_ttW.root"],
            "WWW":   ["TT2L2Nu_2_Muons_WWW.root"],
                     #"TT2L2Nu_1_Muon_1_Electron_WWW.root",
                     #"TT2L2Nu_2_Electrons_WWW.root"],
            "WH":    ["TT2L2Nu_2_Muons_HWminusJets.root",
                      "TT2L2Nu_2_Muons_HWplusJets.root"],
                     #"TT2L2Nu_1_Muon_1_Electron_HWminusJets.root",
                     #"TT2L2Nu_1_Muon_1_Electron_HWplusJets.root",
                     #"TT2L2Nu_2_Electrons_HWminusJets.root",
                     #"TT2L2Nu_2_Electrons_HWplusJets.root"],
            "ttbar": ["TT2L2Nu_2_Muons_ttbar.root"],
                     #"TT2L2Nu_1_Muon_1_Electron_1_ttbar.root",
                     #"TT2L2Nu_1_Muon_1_Electron_2_ttbar.root",
                     #"TT2L2Nu_1_Muon_1_Electron_3_ttbar.root",
                     #"TT2L2Nu_2_Electrons_ttbar.root"],
        }
    }
    files = {
        "sig": ["TT2L2Nu_2_Muons_signal.root",
                "TT2L2Nu_1_Muon_1_Electron_signal.root",
                "TT2L2Nu_2_Electrons_signal.root"],
        "bkg": {
            "WZ":    ["TT2L2Nu_2_Muons_WZ.root",
                      "TT2L2Nu_1_Muon_1_Electron_WZ.root",
                      "TT2L2Nu_2_Electrons_WZ.root"],
            "ttW":   ["TT2L2Nu_2_Muons_ttW.root",
                      "TT2L2Nu_1_Muon_1_Electron_ttW.root",
                      "TT2L2Nu_2_Electrons_ttW.root"],
            "WWW":   ["TT2L2Nu_2_Muons_WWW.root",
                      "TT2L2Nu_1_Muon_1_Electron_WWW.root",
                      "TT2L2Nu_2_Electrons_WWW.root"],
            "WH":    ["TT2L2Nu_2_Muons_HWminusJets.root",
                      "TT2L2Nu_2_Muons_HWplusJets.root",
                      "TT2L2Nu_1_Muon_1_Electron_HWminusJets.root",
                      "TT2L2Nu_1_Muon_1_Electron_HWplusJets.root",
                      "TT2L2Nu_2_Electrons_HWminusJets.root",
                      "TT2L2Nu_2_Electrons_HWplusJets.root"],
            "ttbar": ["TT2L2Nu_2_Muons_ttbar.root",
                      "TT2L2Nu_1_Muon_1_Electron_1_ttbar.root",
                      "TT2L2Nu_1_Muon_1_Electron_2_ttbar.root",
                      "TT2L2Nu_1_Muon_1_Electron_3_ttbar.root",
                      "TT2L2Nu_2_Electrons_ttbar.root"],
        }
    }
    storage   = "rootfilesstorage"
    usage     = "rootfiles"
    outputDir = "analysed/"

    for color in "blue", "cyan", "green", "grey", "magenta", "red", "white", "yellow":
        print(colored(color, color), end="")
    print()
    
    analyser = "src/TopHiggsTrileptonAnalyser.cpp"
    print(f"Analyser: {analyser}")
    with open(analyser, "r") as src:
        lines = list(filter(lambda line: "//(" in line and ")//" in line,
                            map(lambda line: line.strip(), src.readlines())))

    run(f"ls -lh {usage}")
    run(f"ls -lh {storage}")

    steps = {}    
    def stepNcut(line):
        search = re.search("//\(([a-zA-Z0-9_]+)\)//(.*)", line)
        step, instr, full = search.group(1), search.group(2), search.group(0)
        if step in steps:
            steps[step] += [[full, instr]]
        else:
            steps[step] = [[full, instr]]
        return step, instr, full
    lines = list(map(stepNcut, lines))

    for step in sorted(steps, key=steps.get):
        print()
        print()
        print()
        print()
        print("#" * 80)
        print()
        print()
        print()
        print()
        cuts = steps[step]
        ##################################################################################
        # Enable step ``step``
        with open(analyser, "r") as src:
            curr = src.read()
            for full, instr in cuts:
                print(colored(f"{full}\n/*()*/{instr}", "yellow"))
                curr = curr.replace(full, f"/*()*/{instr}") # comment => instruction
        with open(analyser, "w") as src:
            src.write(curr)
        ##################################################################################

        run("make clean")
        print()
        print()
        filenamesStr = fetch(usage)
        filenamesArr = filenamesStr.split("\n")
       #result = "_and_".join(map(lambda fname: fname[8:].replace(".root", ""), filenamesArr)) + f"_{step}.root"
        result = f"signal_{step}.root"
        print("Processing:")
        print(colored(filenamesStr, "green"))
        print("Saving result as:")
        print(colored(result, "green"))
        print("#"*80)
        compileRun()
        print("#"*80)
        run(f"mv analysed/analysed.root {outputDir}{result}")
        [run(f"mv {usage}/{filename} {storage}/") for filename in filenamesArr]
        run(f"ls -lh {usage}")
       #run(f"ls -lh {storage}")

        for (bkg, filenamesArrStorage) in files["bkg"].items():
            print()
            print()
            print()
            print()
            filenamesStrStorage = "\n".join(filenamesArrStorage)
           #result = "_and_".join(map(lambda fname: fname[8:].replace(".root", ""), filenamesArrStorage)) + f"_{step}.root"
            result = f"{bkg}_{step}.root"
            print("Processing:")
            print(colored(filenamesStrStorage, "green"))
            print("Saving result as:")
            print(colored(result, "green"))
            [run(f"mv {storage}/{filename} {usage}/") for filename in filenamesArrStorage]
            run(f"ls -lh {usage}")
            #run(f"ls -lh {storage}")
            print("#"*80)
            compileRun()
            print("#"*80)
            run(f"mv analysed/analysed.root {outputDir}{result}")
            [run(f"mv {usage}/{filename} {storage}/") for filename in filenamesArrStorage]
            run(f"ls -lh {usage}")
            #run(f"ls -lh {storage}")

        [run(f"mv {storage}/{filename} {usage}/") for filename in filenamesArr]
        run(f"ls -lh {usage}")
        #run(f"ls -lh {storage}")

        ##################################################################################
        # Disable step ``step``
        with open(analyser, "r") as src:
            curr = src.read()
            for full, instr in cuts:
                print(colored(f"/*()*/{instr}\n{full}", "yellow"))
                curr = curr.replace(f"/*()*/{instr}", full) # instruction => comment
        with open(analyser, "w") as src:
            src.write(curr)
        ##################################################################################
        content = ""
        with open("jobconfiganalysis.py", "r") as f:
            content = f.read()
        with open("jobconfiganalysis.py", "w") as f:
            f.write(content.replace("'isSignal': False", "'isSignal': True"))
