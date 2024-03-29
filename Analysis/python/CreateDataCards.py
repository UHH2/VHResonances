from Utils import *
from parallelise import *

nParallel     = 25
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
    groups = {}
    groups["Theory"]      = ["NNPDF", "murmuf"]
    groups["Detector"]    = ["lumi", "pu"]
    groups["DeepAk8"]     = ["taggerSF"]
    groups["Tagging"]     = ["btag"] + groups["DeepAk8"]
    groups["Jets"]        = ["JEC", "JER"] + groups["Tagging"]
    if not "inv" in ChannelNames[0]:
        groups["Leptons"] = ["id", "prefiring","trigger"]
    if "muon" in ChannelNames[0]:
        groups["Muon"]    = ["MuonScale", "isolation", "reco", "tracking"]
        groups["Leptons"]+= groups["Muon"]
    groups["Systematics"] = list(set([el for sys in groups.keys() for el in groups[sys]]))
    groups["Signal"]      = []
    groups["Background"]  = []
    for param in Parameters:
        lines.append(param)
        sys = param.split()[0]
        if "sg_" in sys: groups["Signal"].append(sys)
        if "data_" in sys: groups["Background"].append(sys)
    lines.append(separator)
    lines.append(separator)
    groups["Stastical"] = groups["Signal"] + groups["Background"]
    groups["All"]       = groups["Systematics"] + groups["Stastical"]
    for sys in groups.keys():
        line = sys+" group ="
        for el in groups[sys]:
            line += " "+el
        lines.append(line)
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
        # self.fitFunction = "Exp_3"
        # self.fitFunction = "NO"
        # self.mode = "bkg_pred"
        self.fitFunction = "Exp_2"
        # self.mode = "DY_SR"
        # self.mode = "bkg_pred"
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
		    print BR
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
                # if "NNPDF" in sys: continue
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
                if DoControl([""], year+channel+sys, channel, signalName): continue
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
        # print len(self.list_processes)
        # for x in self.list_processes: print x
        if self.RunCombine: parallelise(self.list_processes, nParallel, self.list_logfiles, cwd=True)

    @timeit
    def RunCommandsPerDataCard(self):
        list_processes = {}
        list_logfiles = []
        index = 1
        for histFolder in self.histFolders:
            for collection in self.collections:
                for year in self.years:
                    for channel in self.channels+["leptonchannel"]:
                    # for channel in ["leptonchannel"]:
                        uniqueName = year+"/"+collection+"/"+channel+"/"+histFolder
                        uniqueName, workingDir = self.GetWorkdirName(uniqueName)
                        if not "lepton" in channel: self.ExtractParameters(uniqueName,workingDir,year,channel,histFolder)
                        for signalName in self.RatesSignal[uniqueName.replace("lepton","muon")]:
                            DCname = self.DataCardName(uniqueName,signalName).replace(".txt","")
                            mass = signalName.replace("M","")
                            if mass == "1000": continue
                            # if mass != "3000": continue
                            # for tool,method,freeze in list(itertools.product(["Scan","Impacts","FitDiagnostics","GOF"],["_Asimov","_noAsimov"],["","_Freeze"])):
                            for tool,method,freeze in list(itertools.product(["GOF"],["_Asimov","_noAsimov"],["","_Freeze"])):
                                if "FitDiagnostics" in tool:
                                    if "_Freeze" in freeze: continue
                                if "GOF" in tool:
                                    if "_Freeze" in freeze: continue
                                    if "_Asimov" in method: continue
                                folder = workingDir+tool+method+freeze+"/M"+mass
                                os.system("mkdir -p "+folder)
                                os.system("cp "+workingDir+DCname+".txt "+folder)
                                if "muon" in channel or "ele" in channel or "inv" in channel:
                                    os.system("cp "+workingDir+"ws_*.root "+folder)
                                ExtraCommand  = ["--robustFit", "1"]
                                if "_Asimov" in method:
                                    ExtraCommand += ["-t", "-1"]
                                if "_Freeze" in freeze:
                                    if "Impacts" in tool:
                                        ExtraCommand += ["--freezeNuisanceGroups", "Background"]
                                    if "Scan" in tool:
                                        ExtraCommand += ["--freezeNuisanceGroups", "Systematics"]

                                if "FitDiagnostics" in tool:
                                    command = [folder, "combine","-M","FitDiagnostics", "-d", DCname+".txt", "-m", mass, "-n", signalName+"_"+self.fitFunction, "--saveWorkspace", "--saveShapes", "--saveOverallShapes"]
                                    list_processes.setdefault(tool, []).append(command)
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                elif "Impacts" in tool:
                                    command = [folder, "combineTool.py","-M","Impacts", "-d", DCname+".root", "-m", mass]
                                    list_processes.setdefault("text2workspace", []).append([folder, "text2workspace.py", DCname+".txt", "-m", mass])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    list_processes.setdefault(tool+"_In", []).append(command+ExtraCommand+["--doInitialFit"])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    list_processes.setdefault(tool+"_Do", []).append(command+ExtraCommand+["--doFits"])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    list_processes.setdefault(tool+"_Im", []).append(command+["-o", "Impacts_"+mass+".json"])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    list_processes.setdefault(tool+"_Plot", []).append([folder, "plotImpacts.py","-i","Impacts_"+mass+".json","-o","Impacts_"+mass])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                elif "Scan" in tool:
                                    command = [folder, "combine","-M","MultiDimFit", "-m", mass, "--rMin", "-0.5", "--rMax", "3",]
                                    list_processes.setdefault(tool+"_snapshot", []).append(command+ExtraCommand+["-d", DCname+".txt","-n", method+"_"+tool+"_snapshot", "--saveWorkspace"])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    list_processes.setdefault(tool, []).append(command+ExtraCommand+["higgsCombine"+method+"_"+tool+"_snapshot.MultiDimFit.mH"+mass+".root", "-n", method+"_"+tool+freeze, "--snapshotName", "MultiDimFit", "--algo", "grid", "--points", "30"])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    pytool = self.Path_ANALYSIS+"../../CombineHarvester/CombineTools/scripts/plot1DScan.py"
                                    command = [folder, pytool, "higgsCombine"+method+"_"+tool+freeze+".MultiDimFit.mH"+mass+".root", "-o","Scan"+method+tool+freeze]
                                    list_processes.setdefault(tool+"_Plot", []).append(command)
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    if "_Freeze" in freeze:
                                        command = [folder, pytool, workingDir+tool+method+"/M"+mass+"/higgsCombine"+method+"_"+tool+".MultiDimFit.mH"+mass+".root", "-o","Scan_Combine"+method+tool+freeze]
                                        list_processes.setdefault(tool+"_PlotCombined", []).append(command+["--others", "higgsCombine"+method+"_"+tool+freeze+".MultiDimFit.mH"+mass+".root:FreezeAll:632", "--breakdown", "Syst,Stat"])
                                        list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                elif "GOF" in tool:
                                    command = [folder, "combine","-M","GoodnessOfFit", "-d", DCname+".txt", "-m", mass, "-n", method+"_"+tool, "--algo", "saturated"]
                                    list_processes.setdefault(tool, []).append(command+["-t", "100", "-s", "123456"])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    command = [folder, "combine","-M","GoodnessOfFit", "-d", DCname+".txt", "-m", mass, "-n", method+"_"+tool, "--algo", "saturated"]
                                    list_processes.setdefault(tool, []).append(command)
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1
                                    command = [folder, "combine","-M","GoodnessOfFit", "-d", DCname+".txt", "-m", mass, "-n", method+"_"+tool, "--algo", "saturated"]
                                    list_processes.setdefault(tool+"_Plot", []).append(command+["higgsCombine_noAsimov_GOF.GoodnessOfFit.mH1400.123456.root"])
                                    list_logfiles.append("log_"+str(index)+".txt"); index += 1

        # for command in ["FitDiagnostics","text2workspace","Impacts_In","Impacts_Do","Impacts_Im","Impacts_Plot","Scan_snapshot", "Scan_snapshot", "Scan", "Scan_Plot", "Scan_PlotCombined"]:
        # for command in ["Scan_snapshot", "Scan_snapshot", "Scan", "Scan_Plot", "Scan_PlotCombined"]:
        for command in ["GOF", "GOF_Plot"]:
            for x in list_processes[command]: print x
            print command, len(list_processes[command])
            if self.RunCombine:
                print command
                # parallelise(list_processes[command], nParallel, list_logfiles, cwd=True)


if __name__ == '__main__':
    args = parse_arguments()

    RunCombine    = (True if args.RunCombine=="True" else False)
    print "RunCombine: " + str(RunCombine)
    Method="AsymptoticLimits"
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

    # histFolders = ["DeepAk8_ZHccvsQCD_MD", "DeepAk8_ZHccvsQCD_MD2", "DeepAk8_HccvsQCD", "DeepAk8_H4qvsQCD", "DeepAk8_H4qvsQCD_massdep", "DeepAk8_H4qvsQCD_massdep_HccvsQCD", "tau42"]
    histFolders = ["DeepAk8_ZHccvsQCD_MD", "DeepAk8_H4qvsQCD", "DeepAk8_H4qvsQCD_massdep"]
    histFolders = ["DeepAk8_ZHccvsQCD_MD"]

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
