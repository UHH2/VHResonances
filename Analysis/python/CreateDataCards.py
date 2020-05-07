import sys,os, time
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fileManipulation import *
from parallelise import *

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/macros")
from ModuleRunnerBase import *


nSpace_text   = 45
nSpace_number = 25
def nSpace(word=""):
    return word+" "*(nSpace_number -len(word))
def nSpace_Long(word=""):
    return word+" "*(nSpace_text -len(word))

def MergeListString(myList):
    return " ".join([nSpace_Long(str(x)) if len(myList)>1 and "param" in str(myList[1]).lower() else nSpace(str(x)) for x in myList[:1]]+[nSpace(str(x)) for x in myList[1:]])

def CreateDataCard(path, filename, mode = "", signalName="MC_ZprimeToZH", ChannelNames=["ZprimeZH", "bbchannel"], bkgNames=["MC_DY", "MC_TTbar", "MC_Diboson"], Parameters=[], Observed=[-1,-1], Rates=[[0,0,0],[0,0,0]], year="2016"):
    nChannel = len(ChannelNames)
    nBkg = len(bkgNames)
    nSys = 2 +len(filter(lambda x: "param" in x, Parameters)) # 2 = Lumi and shape for signal
    # nSys = 1 +len(filter(lambda x: "param" in x, Parameters)) # 2 = Lumi
    if nChannel!=len(Observed):
        raise Exception("nChannel!=len(Observed): "+str(nChannel)+" "+str(len(Observed)))
    if nChannel*(nBkg+1)!=np.array(Rates).flatten().shape[0]:
        raise Exception("nChannel*(nBkg+1)!=len(rates): "+str(nChannel*(nBkg+1))+" "+str(np.array(Rates).flatten().shape[0]))
    separator = "\n"+"-"*(2*nSpace_number)
    lines = []
    if isHbb:
        lines.append("# Version of the "+str(ModuleRunnerBase(year).lumi_fb)+"/fb Zprime->Z(ll)H(bb) analysis") #TODO Make it more general for other channels
    else:
        lines.append("# Version of the "+str(ModuleRunnerBase(year).lumi_fb)+"/fb Zprime->Z(ll)H(WW) analysis") #TODO Make it more general for other channels
    lines.append(separator)
    lines.append("imax  "+str(nChannel)+" number of channels")
    lines.append("jmax  "+str(nBkg)+" number of backgrounds")
    lines.append("kmax  "+str(nSys)+" number of nuisance parameters (sources of systematical uncertainties)")
    lines.append(separator)
    for cn in ChannelNames:
        lines.append("shapes * "+cn+" ws_"+mode+".root "+cn+":$PROCESS")
    lines.append(separator)
    lines.append("# name of channels, and number of observed events (total number of event in Data)\n")
    lines.append(nSpace("bin")+MergeListString(ChannelNames))
    lines.append(nSpace("observation")+MergeListString(Observed))
    lines.append(separator)
    lines.append("# name of the channel you are considering, name of the process (signal,bkg,...),")
    lines.append("# process unique ID (positive number for backgrounds, and zero or negative for signal)")
    lines.append("# expected events for each process (total number of events in MC)\n")
    lines.append(nSpace("bin")+MergeListString(sorted(ChannelNames*(nBkg+1))))
    lines.append(nSpace("process")+MergeListString([signalName]+bkgNames)*nChannel)
    lines.append(nSpace("process")+MergeListString(range(0,nBkg+1))*nChannel)
    lines.append(nSpace("rate")+MergeListString(Rates if Rates else [""]*(nBkg+1)*nChannel))
    lines.append(separator)
    lines.append("# list of independent sources of uncertainties, and give their effect (syst. error)\n")
    lines.append(nSpace("lumi     lnN")+MergeListString(([str(ModuleRunnerBase(year).lumi_sys/100+1.)]+["-"]*nBkg)*nChannel))
    lines.append(nSpace("sigma    lnU")+MergeListString((["1.150"]+["-"]*nBkg)*nChannel))
    lines.append(separator)
    for param in Parameters:
        lines.append(param)
    lines.append(separator)
    with open(path+filename, "w") as outputfile:
        for line in lines:
            outputfile.write(line+"\n")


class CreateDataCards:
    def __init__(self,studies = "nominal", RunCombine = True, Method="AsymptoticLimits", extraOptions = ["-t", "-1"], extraOptionsText = "Asimov" , isHbb = False):
        self.studies = studies
        self.RunCombine = RunCombine
        self.isHbb = isHbb
        self.Method = Method
        self.extraOptions = extraOptions
        self.extraOptionsText = extraOptionsText
        # self.extraOptions = ["-t", "-1"] #TODO move to main
        # self.extraOptions+= ["--freezeParameters","sg_p0,sg_p1"] #TODO move to main

        self.AnalysisDir = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/Limits/"+studies+"/"
        if (self.isHbb): self.AnalysisDir += "/Hbb/"
        # self.histFolders = ["btag_DeepBoosted_H4qvsQCD", "NN", "NN_1","NN_2", "CNN", "tau42"]
        # self.histFolders = ["btag_DeepBoosted_H4qvsQCD"]
        # self.histFolders = ["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02"]
        self.histFolders = ["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02"]
        self.years = ["2016", "2017", "2018", "RunII"]
        self.channels = ["muonchannel", "electronchannel"]
        self.collections = ["Puppi"]
        # self.years = ["2016"]
        # self.channels = ["muonchannel"]
        # self.mode = "CB"
        # self.mode = "Exp_3"
        self.mode = "Exp_2"
        self.ResetLists()
    def ResetLists(self):
        self.list_processes = []
        self.list_logfiles = []
        self.Observed = {}
        self.ChannelNames = {}
        self.bkgNames = {}
        self.RatesBkg = {}
        self.RatesSignal = {}
        self.Parameters = {}
        self.ParametersSignal = {}
    def ExtractParameters(self,uniqueName,workingDir,year,collection,channel,histFolder):
        ChannelName = channel+"_"+year
        AnalysisOutput = "OutputFit_"+histFolder+".txt"
        self.ChannelNames.setdefault(uniqueName,{"single":[], "leptonchannel": [], "fullRunII": []})
        self.ChannelNames[uniqueName]["leptonchannel"].append(ChannelName) # TODO not yet implemented. Do we need it?
        self.ChannelNames[uniqueName]["fullRunII"].append(ChannelName) # TODO not yet implemented. Do we need it?
        self.ChannelNames[uniqueName]["single"] = [ChannelName]
        self.Observed.setdefault(uniqueName,[])
        self.Observed[uniqueName] = [-1]
        self.bkgNames.setdefault(uniqueName,[])
        self.RatesBkg.setdefault(uniqueName,{})
        self.RatesSignal.setdefault(uniqueName,{})
        self.Parameters.setdefault(uniqueName,[])
        self.ParametersSignal.setdefault(uniqueName,{})
        with open(workingDir+AnalysisOutput, "U") as file:
            lines = file.readlines()
            for line in lines:
                if "integral" in line and self.mode in line and "bkg_pred" in line:
                    self.bkgNames[uniqueName].append(line.split()[0]  if "h_" in self.mode else line.split()[0])
                    self.RatesBkg[uniqueName].setdefault(ChannelName,[]).append(line.split()[-1] if "h_" in self.mode else line.split()[-4]) #TODO fix for all cases. Before was -3
                if "signal" in line:
                    self.RatesSignal[uniqueName][int(line.split()[0][1:])] = line.split()[-1]
                if "param" in line.lower():
                    # if "rateParam" in line:
                    #     Parameters.append(MergeListString(["TF","rateParam",ChannelName]+[self.mode if "h_" in self.mode else "default_value"]+[line.split()[-1]]))
                    # elif "sg_" in line:
                    if "sg_" in line:
                        self.ParametersSignal[uniqueName].setdefault(line.split()[0],[]).append(MergeListString(line.split()[1:]))
                    elif self.mode in line and not self.mode == line.split()[0] and "bkg_pred" in line:
                        self.Parameters[uniqueName].append(MergeListString(line.split()))
        for index, par in enumerate(self.Parameters[uniqueName]):
            self.Parameters[uniqueName][index] = par.replace("default_value",self.bkgNames[uniqueName][-1])
    def DataCardName(self,uniqueName,signalName):
        return "DataCard_"+uniqueName+"_"+signalName+"_"+self.mode+".txt"
    def DataCardOutPutName(self,workingDir,filename):
        return workingDir+"out/"+filename.replace(".txt","")+"_"+self.Method+"_"+self.extraOptionsText+".out"
    def GetWorkdirName(self,uniqueName):
        workingDir = self.AnalysisDir+"/"+uniqueName+"/datacards/"
        os.system("mkdir -p "+workingDir+"out")
        uniqueName = uniqueName.replace("/","_")
        return uniqueName,workingDir
    def AppendDataCard(self,workingDir,outname,filename,signalName=""):
        # cmd = ["combine", "-M", "AsymptoticLimits", "-d", workingDir+filename, "-t", "-1", "-n", signalName]
        cmd = [workingDir, "combine", "-M", self.Method, "-d", filename, "-n", signalName]
        cmd+= self.extraOptions
        self.list_processes.append(cmd)
        self.list_logfiles.append(outname)
        if self.RunCombine and os.path.exists(outname): os.remove(outname)
    def CombineCards(self,workingDir,filename,signalName,inputCards):
        a = os.system("combineCards.py "+inputCards+" > "+workingDir+filename)
        outname = self.DataCardOutPutName(workingDir,filename)
        self.AppendDataCard(workingDir=workingDir,outname=outname,filename=filename,signalName=signalName)
    @timeit
    def WriteDataCards(self):
        for histFolder in self.histFolders:
            for collection in self.collections:
                for year in self.years:
                    for channel in self.channels:
                        ChannelName = channel+"_"+year
                        uniqueName = year+"/"+collection+"/"+channel+"/"+histFolder
                        uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                        self.ExtractParameters(uniqueName,workingDir,year,collection,channel,histFolder)
                        for sn in self.RatesSignal[uniqueName]:
                            signalName = "M"+str(sn)
                            filename = self.DataCardName(uniqueName,signalName)
                            Rates = [self.RatesSignal[uniqueName][sn]]+self.RatesBkg[uniqueName][ChannelName]
                            CreateDataCard(workingDir, filename, mode=histFolder, signalName= signalName, ChannelNames=self.ChannelNames[uniqueName]["single"], Observed=self.Observed[uniqueName], bkgNames=self.bkgNames[uniqueName],Rates=Rates, Parameters=self.Parameters[uniqueName]+self.ParametersSignal[uniqueName][signalName], year=year)
                            outname = self.DataCardOutPutName(workingDir,filename)
                            self.AppendDataCard(workingDir=workingDir,outname=outname,filename=filename,signalName=signalName)
                            # TODO output tree with different workingDir/name
                            # TODO Parallelize/ Optimize for multi toys (now takes 10 minutes)
    @timeit
    def CombineChannels(self):
        for histFolder in self.histFolders:
            for collection in self.collections:
                for year in self.years:
                    uniqueName = year+"/"+collection+"/leptonchannel/"+histFolder
                    uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                    for sn in self.RatesSignal[uniqueName.replace("leptonchannel","muonchannel")]:
                        signalName = "M"+str(sn)
                        filename = self.DataCardName(uniqueName,signalName)
                        self.CombineCards(workingDir,filename,signalName, inputCards=" ".join([(channel+"_"+year+"="+workingDir+filename).replace("leptonchannel",channel) for channel in self.channels]))
    @timeit
    def CombineYear(self):
        for histFolder in self.histFolders:
            for collection in self.collections:
                for channel in self.channels+["leptonchannel"]:
                    uniqueName = "fullRunII/"+collection+"/"+channel+"/"+histFolder
                    uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                    for sn in self.RatesSignal[uniqueName.replace("fullRunII","2016").replace("leptonchannel","muonchannel")]:
                        signalName = "M"+str(sn)
                        filename = self.DataCardName(uniqueName,signalName)
                        self.CombineCards(workingDir,filename,signalName,inputCards=" ".join([(lambda x,y: x if y!="leptonchannel" else x.replace("lepton","muon")+" "+x.replace("lepton","electron"))((channel+"_"+year+"="+workingDir+filename).replace("fullRunII",year),channel) for year in self.years if year!="RunII"]))
    @timeit
    def RunDataCards(self):
        # TODO if one wants to filter the running cards.
        # list_logfiles   = list(map(lambda x: list_logfiles[x] if ("M30" in list_processes[x][-1]) else "" , range(len(list_processes))))
        # list_processes = list(map(lambda x: list_processes[x] if ("M30" in list_processes[x][-1]) else "" , range(len(list_processes))))
        # list_logfiles  = list(filter(lambda x: x!="", list_logfiles))
        # list_processes = list(filter(lambda x: x!="", list_processes))
        # for i in range(len(list_processes)):
        #     print list_processes[i], list_logfiles[i]
        # print len(list_processes)
        # print " ".join(list_processes[0])
        # print list_logfiles[0]
        if self.RunCombine: parallelise(self.list_processes, 20, self.list_logfiles, cwd=True)


if __name__ == '__main__':
    Method="AsymptoticLimits"
    extraOptions = ["-t", "-1"]
    extraOptionsText = "Asimov"
    isHbb = False
    # TODO Don't create datacards if not needed.
    RunCombine = True
    # RunCombine = False
    studies = "nominal"
    # studies = "noSignalFlatUncertainty"
    DataCards = CreateDataCards(studies=studies, RunCombine=RunCombine, Method=Method, extraOptions=extraOptions, extraOptionsText=extraOptionsText, isHbb=isHbb)
    DataCards.WriteDataCards()
    DataCards.CombineChannels()
    DataCards.CombineYear()
    DataCards.RunDataCards()
