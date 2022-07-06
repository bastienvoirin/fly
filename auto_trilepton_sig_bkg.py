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

################################################################################

if __name__ == "__main__":
    storage = "rootfilesstorage"
    usage = "rootfiles"
    rootfilename = "analysed/analysed_#"

    for color in "blue", "cyan", "green", "grey", "magenta", "red", "white", "yellow":
        print(colored(color, color), end="")
    print()
    
    analyser = "src/TopHiggsTrileptonAnalyser.cpp"
    print(f"analyser: {analyser}")
    with open(analyser, "r") as src:
        lines = list(filter(lambda line: ")//" in line,
                            map(lambda line: line.strip(),
                                src.readlines())))
    steps = {}    
    def stepNcut(line):
        search = re.search("//\((\d+)\)//(.*)", line) # find //(N)//addCut("", ""); lines
        print(search)
        step, instr, full = search.group(1), search.group(2), search.group(0)
        if step in steps:
            steps[step] += [[full, instr]]
        else:
            steps[step] = [[full, instr]]
        return step, instr, full
    lines = list(map(stepNcut, lines))
    print(*lines, sep="\n")
    #print(json.dumps(steps, indent=4))

    for step in sorted(steps, key=steps.get):
        cuts = steps[step]
        ##################################################################################
        # Enable step ``step``
        with open(analyser, "r") as src:
            curr = src.read()
            for full, instr in cuts:
                print(colored(f"{full} => /*()*/{instr}", "yellow"))
                curr = curr.replace(full, f"/*()*/{instr}") # comment => instruction
        with open(analyser, "w") as src:
            src.write(curr)
        ##################################################################################
        # Fetch the filename of the current data set
        currFile = fetch(usage)
        run("make clean")
        print("Processing", colored(currFile, "yellow"))

        compileRun(currFile)
        currRootfilename = rootfilename.replace("#", currFile)
        run(f"mv analysed/analysed.root {currRootfilename}_{step}")

        # Fetch the filename of the other data sets
        otherFiles = fetch(storage).split("\n")

        for otherFile in otherFiles:
            # Swap data sets
            swapDirectories(usage, currFile, storage, otherFile)

            print("Processing", colored(fetch(usage, "ls"), "yellow"))
            compileRun(otherFile)
            otherRootfilename = rootfilename.replace("#", otherFile)
            run(f"mv analysed/analysed.root {otherRootfilename}_{step}")
            # TODO: joint analysis/comparison/stack plots

            # Swap data sets a second time to restore the original state
            swapDirectories(usage, otherFile, storage, currFile)
        ##################################################################################
        # Disable step ``step``
        with open(analyser, "r") as src:
            curr = src.read()
            for full, instr in cuts:
                print(colored(f"/*()*/{instr} => {full}", "yellow"))
                curr = curr.replace(f"/*()*/{instr}", full) # instruction => comment
        with open(analyser, "w") as src:
            src.write(curr)
        ##################################################################################
        content = ""
        with open("jobconfiganalysis.py", "r") as f:
            content = f.read()
        with open("jobconfiganalysis.py", "w") as f:
            f.write(content.replace("'isSignal': False", "'isSignal': True"))
