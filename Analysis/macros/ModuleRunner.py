from glob import glob
import numpy as np

from ModuleRunnerBase import *
sys.path.append("../python")
from Utils import *
from parallelise import *
from CreateConfigFiles import *
from functions import *
import ROOT

#This file contains all the classes and functions to clean up the steer.py code


class ModuleRunner(ModuleRunnerBase):
    def __init__(self,year="2017",controls=[""]):
        ModuleRunnerBase.__init__(self,year)
        self.Module = ""
        self.ConfigFile = ""
        self.defineModules()
        self.CompileModules()
        self.controls = controls
    def RunCommand(self,command="", isPython=False, **kwargs):
        cwd = os.getcwd()
        os.chdir(self.Path_ANALYSIS+"Analysis")
        list_processes = []
        list_logfiles = []
        print "RunCommand:", command
        if isPython:
            process = subprocess.Popen("python python/"+command+".py", shell=True)
            process.wait()
        else:
            if len(kwargs)>0:
                i = 0
                def LoopOver(arg, defaultList):
                    if not arg in kwargs: return [""]
                    else:
                        if len(kwargs[arg])>1: return kwargs[arg]
                        else: return defaultList
                for year in LoopOver("years", self.Samples_dict.keys()+["RunII"]):
                    for collection in LoopOver("Collections", self.Collections):
                        for channel in LoopOver("Channels", self.Channels):
                            for histFolder in LoopOver("histFolders", kwargs["histFolders"]):
                                list_processes.append(["./"+command,histFolder,channel+"channel" if channel!="" else channel,collection,year])
                                list_logfiles.append("log_"+str(i)+".txt")
                                i +=1
                for i in list_processes:
                    print i
                # print "Number of processes", len(list_processes)
                # parallelise(list_processes, 20,list_logfiles)
                parallelise(list_processes, 20)
            else:
                process = subprocess.Popen("./"+command, shell=True)
                process.wait()
        os.chdir(cwd)
    def defineModules(self):
        self.ModuleSamples_dict  = {"Test"               : {"sample": self.Samples_original[self.year], "all": self.Samples_NamesAll},
                                    "GenericCleaning"    : {"sample": self.Samples_original[self.year], "all": self.Samples_NamesAll},
                                    "ProbeNN"            : {"sample": self.Samples_Category[self.year], "all": self.Samples_CategoryAll},
                                    "NeuralNetwork"      : {"sample": self.Samples_Category[self.year], "all": self.Samples_CategoryAll},
                                    "Preselection"       : {"sample": self.Samples_original[self.year], "all": self.Samples_NamesAll},
                                    "Selection"          : {"sample": self.Samples_Category[self.year], "all": self.Samples_CategoryAll},
                                    "SignalRegion"       : {"sample": self.Samples_Category[self.year], "all": self.Samples_CategoryAll},
                                    "SF"                 : {"sample": self.Samples_original[self.year], "all": self.Samples_NamesAll},
                                    "VariableRStudies"   : {"sample": [x for x in self.Samples_original[self.year] if not "DATA".lower() in x.lower()], "all": self.Samples_NamesAll},
                                    }

    def CompileModules(self):
        cwd = os.getcwd()
        os.chdir(self.Path_ANALYSIS)
        process = subprocess.Popen("make -j 20", shell=True)
        process.wait()
        os.chdir(self.Path_ANALYSIS+"Analysis")
        os.system("mkdir -p "+self.Path_ANALYSIS+"Analysis/obj")
        os.system("mkdir -p "+self.Path_ANALYSIS+"Analysis/OtherPlots")
        process = subprocess.Popen("make -j 20", shell=True)
        process.wait()

    def SetModule(self,module, Collections=[], Channels=[], Systematics=[]):
        self.Module     = module
        self.ModuleFile = self.Module + "Module.cxx"
        self.ConfigFile = self.Module + "Config.xml"
        self.ModuleStorage = self.Path_STORAGE+self.year+"/"+self.Module
        self.Samples = self.ModuleSamples_dict[self.Module]["sample"]
        self.AllSamples = self.ModuleSamples_dict[self.Module]["all"]

        self.Collections = Collections if Collections else self.Collections
        self.Channels = Channels if Channels else self.Channels
        self.Systematics = Systematics if Systematics else self.Systematics
        # if "SignalRegion" in self.Module:
        #     self.Samples = np.array(self.Samples)
        #     mask = np.array([self.signal in x for x in self.Samples])
        #     temp = self.Samples[mask]
        #     self.Samples = list(self.Samples[~mask]) + [self.signal+"Tobb"+x[len(self.signal):] for x in temp] + [self.signal+"ToWW"+x[len(self.signal):] for x in temp] + [self.signal+"_extra"+x[len(self.signal):] for x in temp]
        # else:
        #     self.Samples = self.Samples
        if self.Module == "SF":
            # self.Samples = [sample for sample in self.Samples if "DATA" in sample or "MC_TT" in sample]
            self.Samples = [sample for sample in self.Samples if "MC_TT" in sample]
        if "Test" in self.Module or "GenericCleaning" in self.Module or "VariableRStudies" in self.Module:
            self.Collections = ["All"]
        if ( "GenericCleaning" in self.Module or "Test" in self.Module or "NeuralNetwork" in self.Module or "VariableRStudies" in self.Module ):
            self.Channels = ["lepton"]
        print "\n****************************************\t","\n \tModule Name: \t", self.Module, "\n****************************************\t", "\nRunning \t", self.ConfigFile, "\nUsing \t \t", self.ModuleFile, "\nSamples:\t", len(self.Samples), self.Samples, "\n", self.Channels, "\n", self.Systematics, "\n", self.Collections

    def DeleteWorkdirs(self):
        for path in [self.Path_STORAGE+self.year+"/",self.SubmitDir,self.Path_SFRAME]:
            for collection in self.Collections:
                for channel in self.Channels:
                    for syst in self.Systematics:
                        for sample in self.Samples:
                            if all(not control in collection+channel+syst+sample for control in self.controls):
                                continue
                            cmd = "rm -fr "+path+self.Module+"/"+collection+"/"+channel+"channel/"+syst+"/*"+sample+"*"
                            # print cmd
                            a = os.system(cmd)

    def CreateConfigFiles(self):
        for collection in self.Collections:
                for channel in self.Channels:
                    for syst in self.Systematics:
                        if all(not control in self.year+collection+channel+syst for control in self.controls):
                            continue
        		a = os.system("rm -fr " +self.SubmitDir+self.Module+"/"+collection+"/"+channel+"channel/"+syst+"/")
        CreateConfigFiles(self.year, self.Samples, self.AllSamples, self.Collections, self.Channels, self.Systematics, self.controls, self.Path_ANALYSIS, self.SubmitDir+self.Module+"/", self.ConfigFile, self.Path_SFRAME, self.lumi_pb)

    def CondorControl(self,option="", forPlotting=False):
        condor_control(self.SubmitDir, self.Module, self.Samples, self.Collections, self.Channels, self.Systematics, self.controls,option=option, forPlotting=forPlotting, doParallel=(False if "GenericCleaning" in self.Module else True))

    def RunLocal(self,option, skip="",nProcess=20, isNice=True):
        local_run(self.SubmitDir, self.Module, self.Samples, self.Collections, self.Channels, self.Systematics, self.controls, option=option, isNice=isNice, skip=(skip==""),nProcess=nProcess)

    @timeit
    def SecureMerge(self, mergeCategory=False, forPlotting=True):
        extraText = "_noTree" if forPlotting else ""
        list_processes = []
        Samples = self.Samples_Category[self.year] if mergeCategory else self.Samples
        for collection in self.Collections:
            for channel in self.Channels:
                for syst in self.Systematics:
                    middlePath = collection+"/"+channel+"channel/"+syst+"/"
                    for sample in Samples:
                        if all(not control in collection+channel+syst+sample for control in self.controls):
                            continue
                        mode = "MC" if "MC" in sample else "DATA"
                        commonpath = self.ModuleStorage+"/"+middlePath
                        filespath  = commonpath+"workdir_"+self.Module+"_"+sample+"/" if not mergeCategory or self.signal in sample else [commonpath+self.PrefixrootFile+mode+"."+name_+extraText+".root" for name_ in self.Samples_dict[self.year][sample]]
                        newFile    = commonpath+self.PrefixrootFile+mode+"."+sample+extraText
                        newFile   += "_merge.root" if mergeCategory else ".root"
                        if mergeCategory:
                            if len(filespath)<=1 or self.signal in filespath : continue
                            if forPlotting:
                                list_processes.append(["hadd", "-f", "-T", newFile]+filespath)
                            else:
                                list_processes.append(["hadd", "-f", newFile]+ filespath)
                        else:
                            list_ = glob(filespath+"/*root")
                            if len(list_)==1 and forPlotting :
                                list_processes.append(["hadd", "-f", "-T", newFile, list_[0]])
                            else:
                                # list_processes.append([self.Path_ANALYSIS+"Analysis/python/MergeLargeRootFiles.py", str(forPlotting), filespath, newFile])
                                list_processes.append(["hadd", "-f", "-T", newFile]+list_)
        for i in list_processes:
            print i[0:5], "+ others" if len(i)>5 else ""
        print "Number of processes", len(list_processes)
        parallelise(list_processes, 20)

    @timeit
    def StoreModuleOutput(self, process="Move"):
        list_processes = []
        for collection in self.Collections:
            for channel in self.Channels:
                for syst in self.Systematics:
                    middlePath = collection+"/"+channel+"channel/"+syst+"/"
                    commonpath = self.Path_SFRAME+self.Module+"/"+middlePath
                    newPath = self.ModuleStorage+"/"+middlePath
                    a = os.system("mkdir -p "+newPath)
                    for x in glob(commonpath+"*"):
                        list_processes.append(["mv", x, newPath] if process=="Move" else ["cp", "-r", x, newPath])
        print len(list_processes)
        parallelise(list_processes, 20)

    @timeit
    def CreateXml(self):
        for collection in self.Collections:
            for channel in self.Channels:
                for syst in self.Systematics:
                    middlePath = collection+"/"+channel+"channel/"+syst+"/"
                    for sample_ in self.Samples_Category[self.year]:
                        if all(not control in sample_+channel+collection for control in self.controls):
                            continue
                        decays = [""]
                        # if self.Module == "SignalRegion" and self.signal in sample_:
                        #     decays = ["ToWW","Tobb","_extra"]
                        for decay in decays:
                            sample = sample_.replace(self.signal,self.signal+decay)
                            # print sample
                            mode = "MC" if "MC" in sample else "DATA"
                            filePrefix = self.PrefixrootFile+mode+"."
                            commonpath = self.ModuleStorage+"/"+middlePath
                            out = open(commonpath+sample+".xml", 'w')
                            # print commonpath+sample+".xml"
                            nComment = 0
                            for dir in self.Samples_dict[self.year][sample] if not self.signal in sample and self.Samples != self.Samples_Category[self.year] else [sample]:
                                print "Search in", commonpath+"workdir_"+self.Module+"_"+dir+"/*root"
                                for f_ in glob(commonpath+"workdir_"+self.Module+"_"+dir+"/*root"):
                                    try:
                                        ntuple = ROOT.TFile(str(f_))
                                        AnalysisTree = ntuple.Get("AnalysisTree")
                                        isToWrite =  AnalysisTree.GetEntriesFast() > 0
                                        ntuple.Close()
                                    except Exception as e:
                                        isToWrite = False
                                    if not isToWrite: nComment += 1
                                    extraText1 = "" if isToWrite else "<!-- "
                                    extraText2 = "" if isToWrite else " -->"
                                    out.write(extraText1+'<In FileName="'+f_+'" Lumi="0.0"/>'+extraText2+'\n')
                        out.write('<!-- File Commented: '+str(nComment)+' -->\n')
                        out.close()

    @timeit
    def MakePlots(self):
        cwd = os.getcwd()
        os.chdir(self.Path_SPlotter)
        for collection in self.Collections:
            for channel in self.Channels:
                for syst in self.Systematics:
                    if all(not control in self.year+collection+channel+syst for control in self.controls):
                        continue
                    a = os.system("mkdir -p "+self.Path_STORAGE+self.year+"/"+self.Module+"/"+collection+"/"+channel+"channel/"+syst+"/Plots")
                    process = subprocess.Popen("Plots -f Analysis/"+self.Module+"Plotter_"+channel+"channel_"+collection+"_"+syst+"_"+self.year+".steer", shell=True)
                    process.wait()
        os.chdir(cwd)

    def MakeRunII(self, Collections=[], Channels=[], Systematics=[]):
        year = "RunII"
        list_processes = []
        list_processes_plots = []
        for module in ["Preselection","Selection","SignalRegion"]:
            self.SetModule(module, Collections=Collections, Channels=Channels, Systematics=Systematics)
            for collection in self.Collections:
                for channel in self.Channels:
                    for syst in self.Systematics:
                        if all(not control in module+collection+channel+syst for control in self.controls):
                            continue
                        path_RunII = self.Path_STORAGE+year+"/"+self.Module+"/"+collection+"/"+channel+"channel/"+syst+"/"
                        a = os.system("mkdir -p "+path_RunII+"Plots")
                        if self.Module=="SignalRegion":
                            list_processes_plots.append(["Plots", "-f", "Analysis/"+self.Module+"Plotter_"+channel+"channel_"+collection+"_"+syst+"_"+year+".steer" ])
                        for sample in self.Samples_Category[self.year]:
                            mode = "MC" if "MC" in sample else "DATA"
                            filespath = path_RunII+self.PrefixrootFile+mode+"."+sample+"_noTree.root"
                            command = ["hadd", "-f", "-T", filespath.replace(self.year,year)]
                            for histo in glob(filespath.replace("RunII","201*").replace(self.year,"201*")):
                                command.append(histo)
                            list_processes.append(command)
        # for i in list_processes:
        #     print i
        print "Number of processes", len(list_processes)
        parallelise(list_processes, 20)
        # print "Number of processes", len(list_processes_plots)
        # cwd = os.getcwd()
        # os.chdir(self.Path_SPlotter)
        # for proc in list_processes_plots:
        #     process = subprocess.Popen(" ".join(proc), shell=True)
        #     process.wait()
        # os.chdir(cwd)

    @timeit
    def DoChecks(self):
        check = False
        check = True
        list_toCheck = []
        orig_stderr = sys.stderr
        errname = os.getcwd()+"/err_"+self.year+".txt"
        if check:
            with open(errname, 'w') as f_err:
                sys.stderr = f_err
                self.CreateXml()
                sys.stderr = orig_stderr
        for collection in self.Collections:
            for channel in self.Channels:
                for syst in self.Systematics:
                    middlePath = collection+"/"+channel+"channel/"+syst+"/"
                    for sample_ in self.Samples_Category[self.year]:
                        if all(not control in sample_+channel+collection for control in self.controls):
                            continue
                        sample = sample_.replace(self.signal,self.signal)
                        mode = "MC" if "MC" in sample else "DATA"
                        filePrefix = self.PrefixrootFile+mode+"."
                        commonpath = self.ModuleStorage+"/"+middlePath
                        xml = commonpath+sample+".xml"
                        if "muonchannel" in xml and "SingleElectron" in xml: continue
                        if "electronchannel" in xml and "SingleMuon" in xml: continue
                        with open(xml) as out:
                            lines = out.readlines()
                            if len(lines)==0:
                                print "Empty xml:", xml
                            elif len(lines)==1:
                                if "File Commented" in lines[0]:
                                    print "No files found in:", xml
                                else:
                                    print "Unexpected error:", xml
                            else:
                                ncomment = -1
                                if "File Commented" in lines[-1]: ncomment = int(lines[-1].split()[3])
                                else: print "Unexpected error:", xml
                                if ncomment>0 and ncomment==(len(lines)-1):
                                    print "Found files ",ncomment,"out of ",len(lines)-1," in:", xml
                    for sample_ in self.Samples_original[self.year]:
                        if all(not control in sample_+channel+collection for control in self.controls):
                            continue
                        sample = sample_.replace(self.signal,self.signal)
                        mode = "MC" if "MC" in sample else "DATA"
                        filePrefix = self.PrefixrootFile+mode+"."
                        path = self.Path_ANALYSIS+"/config/SubmittedJobs/"+self.year+"/"+self.Module+"/"+middlePath+"workdir_"+self.Module+"_"+sample+"/Stream_"+sample+"/"
                        for err in glob(path+sample+"_*.e*"):
                            num = err.replace(path,"").split(".")[0].replace(sample+"_","")
                            list_toCheck.append(sample)
                            # print err
                            # for x in glob(path+sample+".e*"+num):
                            #     print "\t", x
        print set(list_toCheck)
        mylist = []
        with open(errname, 'r') as f_err:
            for l in f_err.readlines():
                if "probably not closed, trying to recover" in l: mylist.append(l.split()[3])
        print len(mylist)
        # if not check: os.system("rm -fr "+errname)
        # self.ReRunList(mylist)

    @timeit
    def ReRunList(self, mylist=[]):
        Path_SFRAME = self.Path_SFRAME.replace(self.year+"/","") # TODO not nice!!!
        list_processes = []
        list_copy = []
        for x in mylist:
            x = x.replace("//", "/")
            root = x.split("/")[-1]
            workdir = x.replace(root,"").replace(self.Path_STORAGE,"")
            number = root.split("_")[-1].replace(".root","")
            newnumber = str(int(number)+1)
            xml = self.Path_ANALYSIS+"config/SubmittedJobs/"+workdir+root.replace("uhh2.AnalysisModuleRunner.DATA.", "").replace("uhh2.AnalysisModuleRunner.MC.","").replace(number+".root", newnumber+".xml")
            cmd = "mkdir -p "+Path_SFRAME+"/"+workdir
            os.system(cmd)
            list_processes.append(["sframe_main",xml])
            list_copy.append(["mv", Path_SFRAME+"/"+workdir+"/"+root, x])

        for x in list_processes:
            print x
        print len(list_processes)
        parallelise(list_processes, 10)

        for x in list_copy:
            print x
        print len(list_copy)

        parallelise(list_copy, 10)



    @timeit
    def HowToSpeedCondor(self):
        Samples = self.Samples
        timeList = {}
        nFileList = {}
        checks = ["*.o*","*.l*","NF"]
        checks = ["*.o*","NF"]
        for check in checks:
            rt = "Real time"
            if ".l" in check: rt = "Memory (MB)"
            # if ".l*" in check: rt = "Disk (KB)"
            if "NF" in check: rt = "FileName"
            allmax = 0
            for collection in self.Collections:
                for channel in self.Channels:
                    for syst in self.Systematics:
                        middlePath = collection+"/"+channel+"channel/"+syst+"/"
                        for sample in Samples:
                            if all(not control in sample+channel+collection for control in self.controls):
                                continue
                            path = self.SubmitDir+self.Module+"/"+middlePath+"workdir_"+self.Module+"_"+sample+"/Stream_"+sample+"/"+sample+check
                            if "NF" in check: path = self.SubmitDir+self.Module+"/"+middlePath+"/workdir_"+self.Module+"_"+sample+"/"+sample+"*"
                            val = []
                            for el in glob(path):
                                with open(el, "U") as file:
                                    lines = file.readlines()
                                    if "NF" in check:
                                        sec =0
                                    for line in lines:
                                        if not rt in line: continue
                                        if ".l" in check:
                                            # print line.split()[3], line
                                            sec = float(line.split()[3])
                                        elif "NF" in check:
                                            sec +=1
                                        else:
                                            sec = float(line.split()[10])
                                        if not "NF" in check:
                                            val.append(sec)
                                if "NF" in check:
                                    val.append(sec)
                            if len(val)==0 :
                                continue
                            if len(glob(path))==0: continue
                            max_ = np.amax(np.array(val))
                            min_ = np.amin(np.array(val))
                            std_ = np.std(np.array(val))
                            if (max_>(0*3600) and max_<(200*3600) and ".o" in check):
                                print check, "\t", collection, "\t", channel[:4], "\t", syst, "\t", sample, " "*(30-len(sample)), round(max_/3600,2), "\t", round(min_/3600,2), "\t", round(std_*100/max_,2)
                                timeList.setdefault(sample,[]).append(max_)
                            elif (".l" in check and max_>(1*1024)):
                                print check, "\t", collection, "\t", channel[:4], "\t", syst, "\t", sample, " "*(30-len(sample)), round(max_/1024,2), "\t", round(min_/1024,2), "\t", round(std_*100/max_,2)
                            elif ("NF" in check):
                                print check, "\t", collection, "\t", channel[:4], "\t", syst, "\t", sample, " "*(30-len(sample)), round(max_,2), "\t", round(min_,2), "\t", round(std_*100/max_,2), "\t", len(val), "\t", np.sum(val)
                                nFileList.setdefault(sample,[]).append(max_)
                            allmax = np.amax(np.array([allmax,max_]))
                            allstd = np.std(np.array([allmax,max_]))
            print check, allmax/3600, allstd*100/allmax
        for x in timeList:
            print x, " "*(30-len(x)), np.amax(np.array(timeList[x]))/3600, np.amin(np.array(timeList[x]))/3600, np.std(np.array(timeList[x]))/3600, 3.*np.amax(np.array(nFileList[x]))/(np.amax(np.array(timeList[x]))/3600)
