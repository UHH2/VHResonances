from glob import glob
import sys
sys.path.append("../macros")
from ModuleRunnerBase import *
sys.path.append(GenericPath().PersonalCode)
from parallelise import *
from fileManipulation import *

import ROOT

def cont_event(original_dir):
    count = 0
    for sample in glob(original_dir+"*/*/*/*"):
        if not ".xml" in sample:
            continue
        count += 1
    return count

@timeit
def condor_control(SubmitDir, Module, Samples, Collections, Channels, Systematics, controls_ = [""], option="", forPlotting=False, doParallel=True):
    internal_option= {"List": "-l", "Submit": "-s", "Resubmit": "-r", "Add": "-a", "Merge": "-a", "ForceAdd": "-f", "ForceMerge": "-f"}
    list_processes  = []
    count = 0
    commonpath = SubmitDir+"/"+Module
    for collection in Collections:
        for channel in Channels:
            for syst in Systematics:
                middlePath = commonpath+"/"+collection+"/"+channel+"channel/"+syst+"/"
                for sample in Samples:
                    fname = Module+"_"+sample+".xml"
                    if not os.path.isfile(middlePath+fname): continue
                    if all(not control in fname for control in controls_):
                        continue
                    count += 1
                    if not doParallel: os.chdir(middlePath)
                    if option in internal_option:
                        command = ['sframe_batch.py', internal_option[option], fname]
                        if forPlotting:
                            command = ['sframe_batch.py', internal_option[option], "-T", fname]
                    else:
                        command = ['sframe_batch.py', fname]
                    if doParallel:
                        command = [middlePath]+command
                        list_processes.append(command)
                    else:
                        process = subprocess.Popen(command)
                        process.wait()
                    if not doParallel:
                        if "Submit" in option and all_events<100:
                            time.sleep(3)
                        os.chdir(commonpath)
    if doParallel:
        print len(list_processes)
        # for i in list_processes:
        #     print i
        parallelise(list_processes, 48, cwd=True)

def local_run(SubmitDir, Module, Samples, Collections, Channels, Systematics, controls_ = [""], option="", isNice=True, skip=True,nProcess=20):
    list_processes  = []
    i = 0
    commonpath = SubmitDir+"/"+Module
    for collection in Collections:
        for channel in Channels:
            for syst in Systematics:
                middlePath = commonpath+"/"+collection+"/"+channel+"channel/"+syst+"/"
                for sample in Samples:
                    fname = middlePath+"workdir_"+Module+"_"+sample+""
                    if all(not control in fname for control in controls_):
                        continue
                    controls = []
                    if not os.path.isfile(fname+"/missing_files.txt"): continue
                    with open(fname+"/missing_files.txt", "r") as missing_files:
                        lines = missing_files.readlines()
                        if len(lines)==0:
                            continue
                        for line in lines:
                            rootfile = line.split()[0]
                            process = line.split()[1]
                            file = fname+"/"+line.split()[2]
                            with open(file, "r") as xml_file:
                                xml_lines = xml_file.readlines()
                                for xml_line in xml_lines:
                                    if "ENTITY" in xml_line and "OUTDIR" in xml_line:
                                        outdir = xml_line.split()[-1][1:-2]
                                        break
                            if os.path.isfile(outdir+"/"+rootfile):
                                controls.append([rootfile, process])
                                continue
                            if isNice:
                                list_processes.append( ["nice", process, file] )
                            else:
                                list_processes.append( [process, file] )
                            i += 1
                    comment_lines(fname, "/missing_files.txt", controls, remove=True)
    if option == "Check":
        print len(list_processes)
        # for i in list_processes: print i
        return len(list_processes)
    if option == "Local":
        print len(list_processes)
        for i in list_processes: print i
        if(not skip): parallelise(list_processes, nProcess)


def local_run2(original_dir,option, controls_ = [""], isNice=True, skip=True,nProcess=20):
    if option == "NoSplit" :
        list_processes  = []
        for el in glob(original_dir+"*/*/*/*"):
            if not ".xml" in el or "workdir" in el:
                continue
            if all(not control in el for control in controls_):
                continue
            if isNice:
                list_processes.append( ["nice","sframe_main", el] )
            else:
                list_processes.append( ["sframe_main", el] )
        for i in list_processes:
            print i
        parallelise(list_processes, nProcess)
    list_processes  = []
    i = 0
    for el in glob(original_dir+"*/*/*/*"):
        if all( not control in el for control in controls_):
            continue
        if os.path.isdir(el) and "workdir" in el:
            controls = []
            with open(el+"/missing_files.txt", "r") as missing_files:
                lines = missing_files.readlines()
                if len(lines)==0:
                    continue
                for line in lines:
                    rootfile = line.split()[0]
                    process = line.split()[1]
                    file = el+"/"+line.split()[2]
                    with open(file, "r") as xml_file:
                        xml_lines = xml_file.readlines()
                        for xml_line in xml_lines:
                            if "ENTITY" in xml_line and "OUTDIR" in xml_line:
                                outdir = xml_line.split()[-1][1:-2]
                    if os.path.isfile(outdir+"/"+rootfile):
                        controls.append([rootfile, process])
                        continue
                    if isNice:
                        list_processes.append( ["nice", process, file] )
                    else:
                        list_processes.append( [process, file] )
                    i += 1
            comment_lines(el, "/missing_files.txt", controls, remove=True)
    if option == "Check":
        print len(list_processes)
        # for i in list_processes: print i
        return len(list_processes)
    if option == "Local":
        print len(list_processes)
        for i in list_processes: print i
        if(not skip): parallelise(list_processes, nProcess)




def FirstNonEmptyFile(list_):
    for el in list_:
        MyFile = ROOT.TFile.Open(el)
        MyTree = MyFile.Get("AnalysisTree")
        if MyTree.GetEntriesFast()!=0:
            # print MyTree.GetEntriesFast()
            return el
    return ""
