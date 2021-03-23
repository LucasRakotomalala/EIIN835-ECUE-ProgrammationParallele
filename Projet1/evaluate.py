from subprocess import STDOUT, check_output, TimeoutExpired, CalledProcessError
from os import listdir
from os.path import isfile, join
import subprocess
import sys
import random

dataFolder="data/"
binFolder="bin/"
srcFolder="src/"

datasets=[1,2,3]

def runDataset(name, i) :
    output=""
    try:
        binFile=binFolder+name
        if not isfile(binFile) :
            return "0"
        output = subprocess.check_output([binFile, dataFolder+"/test"+str(i)], stderr=STDOUT, timeout=15,universal_newlines=True).strip().replace('\n', '')
        with open(dataFolder+"/result"+str(i)) as f:
            for line in f:
                if line.strip().replace('\n', '')==output:
                    return "2"
                else:
                    return "1"
    except  CalledProcessError as e :
        #we still have to check the output because
        #some student call returns non 0 value in case of success !
        output=e.output.strip().replace('\n', '')
        if (e.returncode < 0) :
            return "0"
        with open(dataFolder+"/result"+str(i)) as f:
            for line in f:
                if line.strip().replace('\n', '')==output:
                    return "2"
                else:
                    return "1"
    except (TimeoutExpired) as e :
        return "0"



def runProject(name) :
    result=""
    for d in datasets :
        result = result + ";" + runDataset(name, d)
    return result

def compileAndRunProject(f):
    resultCompilation=compile(f)
    print(f + ";" + str(resultCompilation) + runProject(f))


def compile(f) :
    name=binFolder + f
    try:
        output = subprocess.check_output(["gcc","-std=c99", "-o" , name,  srcFolder+ f+".c" ,"-lm", "-fopenmp"], stderr=STDOUT, universal_newlines=True).strip().replace('\n', '')
        resultCompilation=0
        if (output==""):
            resultCompilation=1
    except  CalledProcessError as e :
        resultCompilation=0
    return resultCompilation

def compileAndRunProjects() :
    onlyfiles = [f for f in sorted(listdir(srcFolder), key=lambda s: s.lower()) if isfile(join(srcFolder, f))]
    #we shuffle the projects for execution to avoid side effects
    random.shuffle(onlyfiles)
    for f in onlyfiles:
        resultCompilation=compile(f[:-2])
        print(f[:-2] + ";" + str(resultCompilation)  + runProject(f[:-2]))

if len(sys.argv) > 1:
    compileAndRunProject(sys.argv[1])
else :
    compileAndRunProjects()
