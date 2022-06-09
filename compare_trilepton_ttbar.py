#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv, exit
from os import getcwd, chdir, path, system, popen
from subprocess import Popen
import numpy as np
from matplotlib import pyplot as plt
from termcolor import colored

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
def compileRun(filename):
    run("make -j24 > make.log 2>&1") # Redirect stdout & stderr to make.log
    if "signal" in filename:
        content = ""
        with open("jobconfiganalysis.py", "r") as f:
            content = f.read()
        with open("jobconfiganalysis.py", "w") as f:
            f.write(content.replace("'isSignal': False", "'isSignal': True"))
    else:
        content = ""
        with open("jobconfiganalysis.py", "r") as f:
            content = f.read()
        with open("jobconfiganalysis.py", "w") as f:
            f.write(content.replace("'isSignal': True", "'isSignal': False"))
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

# Swap data sets
def swapDirectories(dir1, file1, dir2, file2):
    ls(dir1)
    ls(dir2)
    system(f"mv {dir1}/{file1} {dir2}/{file1}")
    ls(dir1)
    ls(dir2)
    system(f"mv {dir2}/{file2} {dir1}/{file2}")
    ls(dir1)
    ls(dir2)
    return

# Open a ROOT file
def analyse(rootfilename):
    print(colored("Starting analysis", "green"))
    run(f"root -l {rootfilename}")
    print(colored("Done", "green"))
    return

################################################################################

if __name__ == "__main__":
    storage = "rootfilesstorage"
    usage = "rootfiles"
    rootfilename = "analysed/analysed_#"

    for color in "blue", "cyan", "green", "grey", "magenta", "red", "white", "yellow":
        print(colored(color, color), end="")
    print()

    # Fetch the filename of the current data set
    currFile = fetch(usage)
    print("Processing", colored(currFile, "yellow"))

    # Prepare the compilation & analysis
    run("source /cvmfs/cms.cern.ch/slc7_amd64_gcc900/lcg/root/6.24.07-db9b7135424812ed5b8723a7a77f4016/bin/thisroot.sh")
    run("make clean")

    compileRun(currFile)
    currRootfilename = rootfilename.replace("#", currFile)
    run(f"mv analysed/analysed.root {currRootfilename}")
    #analyse(currRootfilename)

    # Fetch the filename of the other data set
    otherFile = fetch(storage)

    # Swap data sets
    swapDirectories(usage, currFile, storage, otherFile)

    print("Processing", colored(fetch(usage, "ls"), "yellow"))
    compileRun(otherFile)
    otherRootfilename = rootfilename.replace("#", otherFile)
    run(f"mv analysed/analysed.root {otherRootfilename}")
    #analyse(otherRootfilename)

    # TODO: joint analysis/comparison/stack plots

    # Swap data sets a second time to restore the original state
    swapDirectories(usage, otherFile, storage, currFile)

    content = ""
    with open("jobconfiganalysis.py", "r") as f:
        content = f.read()
    with open("jobconfiganalysis.py", "w") as f:
        f.write(content.replace("'isSignal': False", "'isSignal': True"))
