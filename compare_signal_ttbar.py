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

# Set current working directory to this script's directory
def setWorkingDir():
    scriptPath = path.realpath(__file__)
    scriptDir = path.dirname(scriptPath)
    chdir(scriptDir)
    return

# Compile and run the analysis program
def compileRun():
    system("make -j24")
    system("./submitanalysisjob.py jobconfiganalysis.py")
    return

def fetch(directory, command):
    chdir(directory)
    filename = popen(command).read().strip()
    gotoParentDir()
    return filename

def ls(directory):
    command = f"ls {directory}"
    util = lambda cmd: popen(cmd).read().strip().replace("\n", " ") or "-"
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

def analyse(rootfilename):
    print(colored("Starting analysis", "green"))
    system(f"root -l {rootfilename}")
    print(colored("Done", "green"))
    return

################################################################################

if __name__ == "__main__":
    storage = "rootfilesstorage"
    usage = "rootfiles"
    rootfilename = "analyzed/analyzed.root"

    for color in "blue", "cyan", "green", "grey", "magenta", "red", "white", "yellow":
        print(colored(color, color), end="")
    print()

    #setWorkingDir()
    #gotoParentDir() # Move up one directory (from /src/fly/src to /src/fly)

    currFile = fetch(usage, "ls")
    print("Processing", colored(currFile, "yellow"))
    system("cms_env")
    system("cmsenv")
    system("source /cvmfs/cms.cern.ch/slc7_amd64_gcc900/lcg/root/6.24.07-db9b7135424812ed5b8723a7a77f4016/bin/thisroot.sh")
    system("make clean")
    compileRun()
    analyse(rootfilename)

    """
    # Fetch the filename of the other data set
    otherFile = fetch(storage, "ls -d TT2L2Nu*")

    # Swap data sets
    swapDirectories(usage, currFile, storage, otherFile)

    print("Processing", colored(fetch(usage, "ls"), "yellow"))
    compileRun()
    analyse(rootfilename)

    # TODO: joint analysis/comparison/stack plots

    # Swap data sets a second time to restore original state
    swapDirectories(usage, otherFile, storage, currFile)
    """
