from Utils import *
from parallelise import *

nSpace_text   = 50
nSpace_number = 25
nSpace_small  = 10
def nSpace_short(word=""):
    return word+" "*(nSpace_small -len(word))
def nSpace(word=""):
    return word+" "*(nSpace_number -len(word))
def nSpace_Long(word=""):
    return word+" "*(nSpace_text -len(word))

def MergeListString(myList):
    return " ".join([nSpace_Long(str(x)) if len(myList)>1 and "param" in str(myList[1]).lower() else nSpace(str(x)) for x in myList[:1]]+[nSpace(str(x)) for x in myList[1:]])

# def lnN(nominal, varUp, varDown):
#     if varUp==0: raise Exception("varUp==0")
#     if varDown==0: raise Exception("varDown==0")
#     return (1+abs(nominal-varUp)/nominal,1+abs(nominal-varDown)/nominal)

def lnN(nominal, variation):
    if variation==0: raise Exception("var==0")
    return 1+abs(nominal-variation)/nominal


def DataCardTemplate(path, filename, mode = "", signalName="MC_ZprimeToZH", ChannelNames=["ZprimeZH", "bbchannel"], bkgNames=["MC_DY", "MC_TTbar", "MC_Diboson"], Parameters=[], Observed=[-1,-1], Rates=[[0,0,0],[0,0,0]], lumi=(100,1)):
    nChannel = len(ChannelNames)
    nBkg = len(bkgNames)
    # nSys = 2 +len(filter(lambda x: "param" in x, Parameters)) +len(filter(lambda x: "lnN" in x, Parameters)) # 2 = Lumi and shape for signal
    nSys = 1 +len(filter(lambda x: "param" in x, Parameters)) +len(filter(lambda x: "lnN" in x, Parameters)) # 1 = Lumi and shape for signal
    # nSys = 1 +len(filter(lambda x: "param" in x, Parameters)) # 2 = Lumi
    if nChannel!=len(Observed):
        raise Exception("nChannel!=len(Observed): "+str(nChannel)+" "+str(len(Observed)))
    if nChannel*(nBkg+1)!=np.array(Rates).flatten().shape[0]:
        raise Exception("nChannel*(nBkg+1)!=len(rates): "+str(nChannel*(nBkg+1))+" "+str(np.array(Rates).flatten().shape[0]))
    separator = "\n"+"-"*(2*nSpace_number)
    lines = []
    if isHbb:
        lines.append("# Version of the "+str(round(lumi[0],1))+"/fb Zprime->Z(ll)H(bb) analysis") #TODO Make it more general for other channels
    else:
        lines.append("# Version of the "+str(round(lumi[0],1))+"/fb Zprime->Z(ll)H(WW) analysis") #TODO Make it more general for other channels
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
    lines.append(nSpace(nSpace_short("lumi")+"lnN")+MergeListString(([str(lumi[1]/100+1.)]+["-"]*nBkg)*nChannel))
    # lines.append(nSpace(nSpace_short("sigma")+"lnN")+MergeListString((["1.150"]+["-"]*nBkg)*nChannel))
    # lines.append(nSpace(nSpace_short("sigma")+"lnN")+MergeListString((["1.0001"]+["-"]*nBkg)*nChannel))
    lines.append(separator)
    for param in Parameters:
        lines.append(param)
    lines.append(separator)
    with open(path+filename, "w") as outputfile:
        for line in lines:
            outputfile.write(line+"\n")


class CreateDataCards(VariablesBase):
    def __init__(self, histFolders, years, channels, collections, studies = "nominal", RunCombine = True, RunSystematics=True, Method="AsymptoticLimits", extraOptions = ["-t", "-1"], extraOptionsText = "Asimov" , isHbb = False):
        VariablesBase.__init__(self)
        self.studies = studies
        self.RunCombine = RunCombine
        self.RunSystematics = RunSystematics
        self.isHbb = isHbb
        self.Method = Method
        self.extraOptions = extraOptions
        self.extraOptionsText = extraOptionsText
        self.histFolders = histFolders
        self.years = years
        self.channels = channels
        self.collections = collections
        print "\n"+("*"*40)+"\t","\n \tModule Name: \t", type(self).__name__, "\n"+("*"*40)+"\t", "\nUsing \t \t", self.years, "\n", self.histFolders, "\n", self.channels, "\n"
        self.AnalysisDir = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/Limits/"+self.studies+"/"
        if (self.isHbb): self.AnalysisDir += "/Hbb/"

        # self.fitFunction = "CB"
        self.fitFunction = "Exp_3"
        self.fitFunction = "NO"
        # self.mode = "bkg_pred"
        self.fitFunction = "Exp_2"
        # self.mode = "DY_SR"
        self.mode = "bkg_pred"
        self.mode = "data"
        self.BRs= {"invisible": 0.2, "muon": 0.1, "electron": 0.1, "chargedlepton":0.1}
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

    def ExtractParameters(self,uniqueName,workingDir,year,channel,histFolder):
        AnalysisOutput = "OutputFit_"+histFolder+".txt"
        ChannelName = self.ChannelName(channel,year)
        # self.ChannelNames.setdefault(uniqueName,{"single":[], "leptonchannel": [], "fullRunII": []})
        self.ChannelNames.setdefault(uniqueName,{"single":[]})
        # self.ChannelNames[uniqueName]["leptonchannel"].append(ChannelName) # TODO not yet implemented. Do we need it?
        # self.ChannelNames[uniqueName]["fullRunII"].append(ChannelName) # TODO not yet implemented. Do we need it?
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
                if "data_CR" in line: continue
                if "integral" in line and self.fitFunction in line and self.mode in line:
                    self.bkgNames[uniqueName].append(line.split()[0]  if "h_" in self.fitFunction else line.split()[0])
                    self.RatesBkg[uniqueName].setdefault(ChannelName,[]).append(line.split()[-1] if "h_" in self.fitFunction else line.split()[-2]) #TODO fix for all cases. Before was -3
                if "signal" in line:
                    BR = self.BRs["invisible"] if "invisible" in uniqueName else self.BRs["chargedlepton"]
                    if not any(syst in line for syst in self.Systematics+self.Systematics_Scale):
                        self.RatesSignal[uniqueName][line.split()[0]] = str(float(line.split()[-1])*BR)
                    else:
                        self.ParametersSignal[uniqueName].setdefault(line.split()[0],[]).append(MergeListString(["rate",str(float(line.split()[-1])*BR)]))
                if "param" in line.lower():
                    # if "rateParam" in line:
                    #     Parameters.append(MergeListString(["TF","rateParam",ChannelName]+[self.fitFunction if "h_" in self.fitFunction else "default_value"]+[line.split()[-1]]))
                    # elif "sg_" in line:
                    if "sg_" in line:
                        self.ParametersSignal[uniqueName].setdefault(line.split()[0],[]).append(MergeListString(line.split()[1:]))

                    elif self.fitFunction in line and not self.fitFunction == line.split()[0] and self.mode in line:
                        self.Parameters[uniqueName].append(MergeListString(line.split()))
        for index, par in enumerate(self.Parameters[uniqueName]):
            self.Parameters[uniqueName][index] = par.replace("default_value",self.bkgNames[uniqueName][-1])

    def ChannelName(self,channel,year):
        return channel.replace("channel","")+"_"+year

    def DataCardName(self,uniqueName,signalName):
        return "DataCard_"+uniqueName+"_"+signalName+"_"+self.fitFunction+"_"+self.extraOptionsText+".txt"

    def DataCardOutPutName(self,workingDir,filename):
        return workingDir+"out/"+filename.replace(".txt","")+"_"+self.Method+".out"

    def GetWorkdirName(self,uniqueName):
        workingDir = self.AnalysisDir+"/"+uniqueName+"/datacards/"
        os.system("mkdir -p "+workingDir+"out")
        uniqueName = uniqueName.replace("/","_")
        return uniqueName,workingDir

    def AppendDataCard(self,workingDir,outname,filename,signalName=""):
        # cmd = ["combine", "-M", "AsymptoticLimits", "-d", workingDir+filename, "-t", "-1", "-n", signalName]
        cmd = [workingDir, "combine", "-M", self.Method, "-d", filename, "-n", signalName+"_"+self.fitFunction]
        cmd+= self.extraOptions
        self.list_processes.append(cmd)
        self.list_logfiles.append(outname)
        if self.RunCombine and os.path.exists(outname): os.remove(outname)

    def CombineCards(self,workingDir,filename,signalName,inputCards):
        if not RunCombine: a = os.system("combineCards.py "+inputCards+" > "+workingDir+filename)
        outname = self.DataCardOutPutName(workingDir,filename)
        self.AppendDataCard(workingDir=workingDir,outname=outname,filename=filename,signalName=signalName)

    def CreateDataCard(self, path, filename, mode = "", signalName="MC_ZprimeToZH", uniqueName ="",channel="muonchannel",year="2016"):
        ChannelNames = self.ChannelNames[uniqueName]["single"]
        Observed     = self.Observed[uniqueName]
        bkgNames     = self.bkgNames[uniqueName]
        Rates        = [self.RatesSignal[uniqueName][signalName]]+self.RatesBkg[uniqueName][self.ChannelName(channel,year)]
        lumi         = (float(self.lumi_map[year]["lumi_fb"]),float(self.lumi_map[year]["uncertainty"]))
        Parameters   = self.Parameters[uniqueName]+self.ParametersSignal[uniqueName][signalName]
        # if self.RunSystematics:
        #     for sys in self.Systematics + self.Systematics_Scale:
        #         for var in ["Up","Down"]:
        #             if signalName+sys+var in self.ParametersSignal[uniqueName]:
        #                 Parameters += self.ParametersSignal[uniqueName][signalName+sys+var]
        if self.RunSystematics:
            for sys in self.Systematics + self.Systematics_Scale:
                if sys=="nominal": continue
                nominal = float(self.RatesSignal[uniqueName][signalName])
                varUp = 0
                varDown = 0
                # print uniqueName, signalName+sys, nominal,
                for var in ["Up","Down"]:
                    if signalName+sys+var in self.ParametersSignal[uniqueName]:
                        for par in self.ParametersSignal[uniqueName][signalName+sys+var]:
                            if not "rate" in par: continue
                            # print var,
                            if var == "Up": varUp = float(par.split()[-1])
                            else: varDown = float(par.split()[-1])
                            # print (varUp, varDown),
                #     else: #TODO Fix me
                #         varUp = nominal
                #         varDown = nominal
                if not "muon"  in channel and "MuonScale" in sys: continue
                if not "muon"  in channel and "isolation" in sys: continue
                if not "muon"  in channel and "tracking"  in sys: continue
                if not "muon"  in channel and "reco"      in sys: continue
                if "invisible" in channel and "pu"        in sys: continue
                if "invisible" in channel and "btag"      in sys: continue
                if "invisible" in channel and "prefiring" in sys: continue
                if "invisible" in channel and "id"        in sys: continue
                if "invisible" in channel and "trigger"   in sys: continue
                if nominal==0:#TODO Fix me
                    print channel, year,sys,signalName, nominal,varUp,varDown
                    continue
                if "M600" == signalName and (varUp==0 or varDown==0):#TODO Fix me
                    print channel, year,sys,signalName, nominal,varUp,varDown
                    continue
                # if varUp==0: varUp = nominal
                # if varDown==0: varDown = nominal
                # print filename,channel, year,sys,signalName, nominal,varUp,varDown, lnN(nominal,varUp,varDown)
                # Parameters += [nSpace(nSpace_short(sys)+"lnN")+MergeListString(([str(lnN(nominal,varUp,varDown)[0])+"/"+str(lnN(nominal,varUp,varDown)[1])]+["-"]*len(bkgNames))*len(ChannelNames))]
                Parameters += [nSpace(nSpace_short(sys)+"lnN")+MergeListString(([str(lnN(nominal,varDown))+"/"+str(lnN(nominal,varUp))]+["-"]*len(bkgNames))*len(ChannelNames))]
                        # Parameters += self.ParametersSignal[uniqueName][signalName+sys+var]
        DataCardTemplate(path, filename, mode = mode , signalName = signalName, ChannelNames = ChannelNames, bkgNames = bkgNames, Parameters = Parameters, Observed = Observed, Rates = Rates, lumi = lumi)

    @timeit
    def WriteDataCards(self):
        for histFolder in self.histFolders:
            for collection in self.collections:
                for year in self.years:
                    for channel in self.channels:
                        uniqueName = year+"/"+collection+"/"+channel+"/"+histFolder
                        uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                        self.ExtractParameters(uniqueName,workingDir,year,channel,histFolder)
                        for signalName in self.RatesSignal[uniqueName]:
                            filename = self.DataCardName(uniqueName,signalName)
                            if not RunCombine: self.CreateDataCard(workingDir, filename, mode = histFolder, signalName = signalName, uniqueName = uniqueName, channel = channel, year = year)
                            outname = self.DataCardOutPutName(workingDir,filename)
                            self.AppendDataCard(workingDir = workingDir, outname = outname, filename = filename, signalName = signalName)
                            # TODO output tree with different workingDir/name
                            # TODO Parallelize/ Optimize for multi toys (now takes 10 minutes)

    @timeit
    def CombineChannels(self, isLepton="false"):
        for histFolder in self.histFolders:
            for collection in self.collections:
                for year in self.years:
                    uniqueName = year+"/"+collection+"/chargedleptonchannel/"+histFolder
                    uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                    for signalName in self.RatesSignal[uniqueName.replace("chargedleptonchannel","muonchannel")]:
                        if len(list(filter(lambda x: x!="invisiblechannel", self.channels)))>=2:
                            filename = self.DataCardName(uniqueName,signalName)
                            self.CombineCards(workingDir,filename,signalName, inputCards=" ".join([(channel+"_"+year+"="+workingDir+filename).replace("chargedleptonchannel",channel) for channel in self.channels if channel!="invisiblechannel"]))
                    uniqueName = year+"/"+collection+"/leptonchannel/"+histFolder
                    uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                    for signalName in self.RatesSignal[uniqueName.replace("leptonchannel","muonchannel")]:
                        if len(self.channels)==3:
                            filename = self.DataCardName(uniqueName,signalName)
                            self.CombineCards(workingDir,filename,signalName, inputCards=" ".join([(channel+"_"+year+"="+workingDir+filename).replace("leptonchannel",channel) for channel in self.channels]))

    @timeit
    def CombineYear(self):
        for histFolder in self.histFolders:
            for collection in self.collections:
                for channel in self.channels+(["leptonchannel"] if self.channels!=["invisiblechannel"] else []):
                    uniqueName = "fullRunII/"+collection+"/"+channel+"/"+histFolder
                    uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                    key = uniqueName.replace("fullRunII","2016").replace("leptonchannel","muonchannel")
                    if (channels==["invisiblechannel"]): key = key.replace("muonchannel", "invisiblechannel")
                    for signalName in self.RatesSignal[key]:
                        filename = self.DataCardName(uniqueName,signalName)
                        self.CombineCards(workingDir,filename,signalName,inputCards=" ".join([(lambda x,y: x if y!="leptonchannel" else x.replace("lepton","muon")+" "+x.replace("lepton","electron"))((channel+"_"+year+"="+workingDir+filename).replace("fullRunII",year),channel) for year in self.years if year!="RunII"]))

    @timeit
    def RunDataCards(self):
        # TODO if one wants to filter the running cards.
        self.list_logfiles  = list(map(lambda x: self.list_logfiles[x]  if ("RunII" in self.list_processes[x][0] and not "full" in self.list_processes[x][0]) else "" , range(len(self.list_processes))))
        self.list_processes = list(map(lambda x: self.list_processes[x] if ("RunII" in self.list_processes[x][0] and not "full" in self.list_processes[x][0]) else "" , range(len(self.list_processes))))
        # self.list_processes = list(map(lambda x: self.list_logfiles[x] if ("inv" in self.list_logfiles[x][0]) else "" , range(len(self.list_logfiles))))
        # self.list_processes = list(map(lambda x: self.list_processes[x] if ("inv" in self.list_processes[x][0]) else "" , range(len(self.list_processes))))
        # self.list_logfiles  = list(map(lambda x: self.list_logfiles[x]  if ("lepton" in self.list_processes[x][0] and "RunII" in self.list_processes[x][0] and not "full" in self.list_processes[x][0]) else "" , range(len(self.list_processes))))
        # self.list_processes = list(map(lambda x: self.list_processes[x] if ("lepton" in self.list_processes[x][0] and "RunII" in self.list_processes[x][0] and not "full" in self.list_processes[x][0]) else "" , range(len(self.list_processes))))
        # # self.list_logfiles  = list(map(lambda x: self.list_logfiles[x]  if (not "RunII" in self.list_processes[x][0] and not "full" in self.list_processes[x][0]) else "" , range(len(self.list_processes))))
        # self.list_processes = list(map(lambda x: self.list_processes[x] if (not "RunII" in self.list_processes[x][0] and not "full" in self.list_processes[x][0]) else "" , range(len(self.list_processes))))
        self.list_logfiles  = list(filter(lambda x: x!="", self.list_logfiles))
        self.list_processes = list(filter(lambda x: x!="", self.list_processes))
        # for i in range(len(self.list_processes)):
        #     print self.list_processes[i], self.list_logfiles[i]
        print len(self.list_processes)
        for x in self.list_processes: print x
        # print " ".join(self.list_processes[0])
        # print self.list_processes[0]
        if self.RunCombine: parallelise(self.list_processes, 20, self.list_logfiles, cwd=True)

    @timeit
    def RunCommandsPerDataCard(self):
        list_processes = {}
        list_logfiles = []
        index = 1
        for histFolder in self.histFolders:
            for collection in self.collections:
                for year in self.years:
                    for channel in self.channels+["leptonchannel"]:
                        uniqueName = year+"/"+collection+"/"+channel+"/"+histFolder
                        uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                        if not "lepton" in channel: self.ExtractParameters(uniqueName,workingDir,year,channel,histFolder)
                        for signalName in self.RatesSignal[uniqueName.replace("lepton","muon")]:
                            DCname = self.DataCardName(uniqueName,signalName).replace(".txt","")
                            mass = signalName.replace("M","")
                            command = [workingDir]
                            list_processes.setdefault("text2workspace", []).append([workingDir, "text2workspace.py", DCname+".txt", "-m", mass])
                            list_processes.setdefault("combineTool_In", []).append([workingDir, "combineTool.py","-M","Impacts", "-d", DCname+".root", "-m", mass, "--doInitialFit", "--robustFit", "1"])
                            list_processes.setdefault("combineTool_Do", []).append([workingDir, "combineTool.py","-M","Impacts", "-d", DCname+".root", "-m", mass, "--robustFit", "1", "--doFits"])
                            list_processes.setdefault("combineTool_Im", []).append([workingDir, "combineTool.py","-M","Impacts", "-d", DCname+".root", "-m", mass, "-o", "impacts_"+mass+".json"])
                            list_processes.setdefault("plotImpacts",    []).append([workingDir, "plotImpacts.py","-i","impacts_"+mass+".json","-o","impacts_"+mass])
                            list_logfiles.append("log_"+str(index)+".txt")
                            index += 1
        for command in ["text2workspace","combineTool_In","combineTool_Do","combineTool_Im","plotImpacts"]:
            print command, len(list_processes[command])
            for x in list_processes[command]: print x[1:]
            parallelise(list_processes[command], 20, list_logfiles, cwd=True)


if __name__ == '__main__':
    args = parse_arguments()

    RunCombine    = (True if args.RunCombine=="True" else False)
    print "RunCombine: " + str(RunCombine)

    # extraOptions = ["-t", "-1"]
    # extraOptionsText = "Asimov"
    extraOptions = ["--run", "both"]
    extraOptionsText = "Expected"
    # extraOptions = ["--run", "observed"]
    # extraOptionsText = "Observed"
    # Method="FitDiagnostics"
    # extraOptions = ["--saveWorkspace", "--saveShapes"]
    # extraOptionsText = "Diagnostics"
    isHbb = False
    RunSystematics = True
    # RunSystematics = False
    if not RunSystematics: extraOptionsText += "NoSys"
    # if not RunSystematics: extraOptionsText += "NoSys0"
    # if RunSystematics: extraOptionsText += "Sys0"
    studies = "nominal"
    # studies = "noSignalFlatUncertainty"

    # histFolders = ["DeepAk8_H4qvsQCD_massdep", "DeepAk8_ZHccvsQCD_MD", "DeepAk8_HccvsQCD_MD", "DeepAk8_H4qvsQCD_MD",
    #                "DeepAk8_H4qvsQCD_massdep_HccvsQCD_MD", "DeepAk8_H4qvsQCD", "DeepAk8_HccvsQCD", "DeepAk8_ZHccvsQCD",
    #                "DeepAk8_H4qvsQCD_massdep_HccvsQCD", "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD", "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD_MD"]

    # histFolders = ["DeepAk8_HccvsQCD", "DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD", "DeepAk8_ZHccvsQCD_MD2"]
    # histFolders = ["DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD2"]
    # histFolders = ["DeepAk8_ZHccvsQCD_MD2"]
    histFolders = ["DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD2"]
    # histFolders = ["DeepAk8_HccvsQCD2"]

    # histFolders = ["DeepAk8_ZHccvsQCD_MD"]

    # years = ["2016", "2017", "2018", "RunII"]
    channels = ["muonchannel", "electronchannel", "invisiblechannel"]
    # channels = ["muonchannel", "electronchannel"]
    # channels = ["invisiblechannel"]
    collections = ["Puppi"]
    years = ["RunII"]
    # channels = ["muonchannel"]


    DataCards = CreateDataCards(histFolders, years, channels, collections, studies=studies, RunCombine=RunCombine, RunSystematics=RunSystematics, Method=Method, extraOptions=extraOptions, extraOptionsText=extraOptionsText, isHbb=isHbb)
    DataCards.WriteDataCards()
    DataCards.CombineChannels()
    # DataCards.CombineYear()
    DataCards.RunDataCards()
    # DataCards.RunCommandsPerDataCard()




# name="DataCard_RunII_Puppi_muonchannel_btag_DeepBoosted_H4qvsQCDmassdep_x3_M1000_Exp_2_simple2"
# mpoint="2000"
# text2workspace.py ${name}.txt -m ${mpoint}
# combineTool.py -M Impacts -d ${name}.root -m ${mpoint} --doInitialFit --robustFit 1
# combineTool.py -M Impacts -d ${name}.root -m ${mpoint} --robustFit 1 --doFits
# combineTool.py -M Impacts -d ${name}.root -m ${mpoint} -o impacts.json
# plotImpacts.py -i impacts.json -o impacts_simple2

#
# mpoint="2000"
# name="DataCard_RunII_Puppi_leptonchannel_btag_DeepBoosted_H4qvsQCD_cc_M${mpoint}_Exp_2_Expected"
# text2workspace.py ${name}.txt -m ${mpoint}
# combineTool.py -M Impacts -d ${name}.root -m ${mpoint} --doInitialFit --robustFit 1 --rMin 0 --rMax 10
# combineTool.py -M Impacts -d ${name}.root -m ${mpoint} --robustFit 1 --doFits --rMin 0 --rMax 10
# combineTool.py -M Impacts -d ${name}.root -m ${mpoint} -o impacts_${mpoint}.json
# plotImpacts.py -i impacts_${mpoint}.json -o impacts_${mpoint}
