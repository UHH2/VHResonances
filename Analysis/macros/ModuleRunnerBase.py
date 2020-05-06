import os
import sys

class GenericPath:
    def __init__(self):
        self.user = os.environ["USER"]
        self.cmssw_base = os.environ["CMSSW_BASE"]
        if "beegfs" in os.getcwd():
            self.PersonalCode = "/beegfs/desy/user/"+self.user+"/PersonalCode/"
            self.isMaxwell = True
        else:
            self.PersonalCode = self.cmssw_base+"/src/UHH2/PersonalCode/"
            self.isMaxwell = False
        self.Path_NFS       = "/nfs/dust/cms/user/"+self.user+"/"
        self.Path_Maxwell   = "/beegfs/desy/user/"+os.environ["USER"]+"/"
        self.Path_ANALYSIS  = self.cmssw_base+"/src/UHH2/HiggsToWWTagger/"
        self.Path_SFRAME    = self.Path_NFS+"sframe_all/"
        self.Path_STORAGE   = self.Path_NFS+"WorkingArea/File/Analysis/"
        self.Path_SPlotter  = self.Path_NFS+"WorkingArea/SFramePlotter/"
    def Get(self, source):
        return getattr(self,source)
    def ExtractVariableFromConstants(self):
        with open(self.Path_ANALYSIS+"include/constants.hpp") as f_:
            self.lumi_map = {}
            line = f_.readline()
            keeploop = False
            while line:
                while "lumi_map" in line or keeploop:
                    keeploop = True
                    line = f_.readline()
                    if '};' in line:
                        keeploop = False; continue
                    if "201" in line.split()[1] or "RunII" in line.split()[1]:
                        year = str(line.split()[1][1:-2])
                        line = f_.readline()
                        while not "}}," in line:
                            self.lumi_map.setdefault(year,{}).setdefault(str(line.split()[1][1:-2]),str(line.split()[-1][:-2]))
                            line = f_.readline()
                line = f_.readline()

class ModuleRunnerBase(GenericPath):
    def __init__(self,year="2017"):
        GenericPath.__init__(self)
        self.year = year
        self.PrefixrootFile = "uhh2.AnalysisModuleRunner."
        self.defineDirectories()
        self.lumi_fb  = round(float(self.lumi_map[year]["lumi_fb"]),1)
        self.lumi_pb  = int(self.lumi_map[self.year]["lumi_pb"])
        self.lumi_sys = round(float(self.lumi_map[year]["uncertainty"]),1)
        # prettydic(self.__dict__)
        self.defineSamples()
    def defineDirectories(self):
        self.Path_SFRAME        = self.Path_SFRAME+self.year+"/"
        self.ConfigDir          = self.Path_ANALYSIS+"config/"
        self.SubmitDir          = self.ConfigDir+"SubmittedJobs/"+self.year+"/"
        self.ExtractVariableFromConstants()
    def defineSamples(self):
        self.Collections = ["Puppi"]
        self.Channels = ["muon", "electron"]
        self.Systematics = ["nominal", "JER_up", "JER_down", "JEC_up", "JEC_down"]
        self.signal = "MC_ZprimeToZH"
        self.SignalSamples = ["MC_ZprimeToZH_M600", "MC_ZprimeToZH_M800", "MC_ZprimeToZH_M1000", "MC_ZprimeToZH_M1200", "MC_ZprimeToZH_M1400", "MC_ZprimeToZH_M1600", "MC_ZprimeToZH_M1800", "MC_ZprimeToZH_M2000", "MC_ZprimeToZH_M2500", "MC_ZprimeToZH_M3000", "MC_ZprimeToZH_M3500", "MC_ZprimeToZH_M4000", "MC_ZprimeToZH_M4500", "MC_ZprimeToZH_M5000", "MC_ZprimeToZH_M5500", "MC_ZprimeToZH_M6000", "MC_ZprimeToZH_M7000", "MC_ZprimeToZH_M8000" ]
        Samples_dict_ = {}
        Samples_dict_["2016"] = {
            #Matching Level
            "MC_QCD"                : ["MC_DY_HT70to100", "MC_DY_HT100to200", "MC_DY_HT200to400", "MC_DY_HT400to600", "MC_DY_HT600to800", "MC_DY_HT800to1200", "MC_DY_HT1200to2500", "MC_DY_HT2500toInf"],
            #Analysis Level
            "MC_DY"                 : ["MC_DY_HT70to100", "MC_DY_HT100to200", "MC_DY_HT200to400", "MC_DY_HT400to600", "MC_DY_HT600to800", "MC_DY_HT800to1200", "MC_DY_HT1200to2500", "MC_DY_HT2500toInf"],
            "MC_TTbar"              : ["MC_TTbar"],
            "MC_WW_incl"            : ["MC_WW"],
            "MC_WZ_incl"            : ["MC_WZ"],
            "MC_ZZ_incl"            : ["MC_ZZ"],
            "MC_WW"                 : ["MC_WWTo4Q", "MC_WWToLNuQQ", "MC_WWTo2L2Nu"],
            "MC_WZ"                 : ["MC_WZToLNu2Q", "MC_WZTo2Q2Nu", "MC_WZTo2L2Q", "MC_WZTo1L3Nu", "MC_WZTo3LNu"],
            "MC_ZZ"                 : ["MC_ZZTo4Q", "MC_ZZTo2Q2Nu", "MC_ZZTo2L2Q", "MC_ZZTo2L2Nu", "MC_ZZTo4L"],
            "DATA_SingleElectron"   : ["DATA_SingleElectron_RunB", "DATA_SingleElectron_RunC", "DATA_SingleElectron_RunD", "DATA_SingleElectron_RunE", "DATA_SingleElectron_RunF", "DATA_SingleElectron_RunG", "DATA_SingleElectron_RunH" ],
            "DATA_SingleMuon"       : ["DATA_SingleMuon_RunB", "DATA_SingleMuon_RunC", "DATA_SingleMuon_RunD", "DATA_SingleMuon_RunE", "DATA_SingleMuon_RunF", "DATA_SingleMuon_RunG", "DATA_SingleMuon_RunH" ],
            }
        Samples_dict_["2017"] = {
            #Matching Level
            "MC_HWW"                : ["MC_HZ_HiggsToWWZToLL"],
            "MC_QCD"                : ["MC_QCD_Pt170to300", "MC_QCD_Pt300to470", "MC_QCD_Pt470to600", "MC_QCD_Pt600to800", "MC_QCD_Pt800to1000", "MC_QCD_Pt1000to1400", "MC_QCD_Pt1400to1800"],
            # "MC_QCD"                : ["MC_QCD_Pt15to30", "MC_QCD_Pt30to50", "MC_QCD_Pt50to80", "MC_QCD_Pt80to120", "MC_QCD_Pt120to170", "MC_QCD_Pt170to300", "MC_QCD_Pt300to470", "MC_QCD_Pt470to600", "MC_QCD_Pt600to800", "MC_QCD_Pt800to1000", "MC_QCD_Pt1000to1400", "MC_QCD_Pt1400to1800", "MC_QCD_Pt1800to2400", "MC_QCD_Pt2400to3200", "MC_QCD_Pt3200toInf"],
            "MC_Top"                : ["MC_TTToHadronic", "MC_TTToSemiLeptonic"],
            "MC_W"                  : ["MC_WZ", "MC_WJetsToQQ_HT400to600", "MC_WJetsToQQ_HT600to800"],
            "MC_Z"                  : ["MC_ZZ", "MC_WZ_Zmatch", "MC_ZJetsToQQ_HT800toInf"],
            #Analysis Level
            "MC_ZH_HWW"             : ["MC_HZ_HiggsToWWZToLL"],
            "MC_TTbar"              : ["MC_TTToHadronic", "MC_TTToSemiLeptonic"],
            "MC_WZ"                 : ["MC_WZ"],
            "MC_WJetsToQQ"          : ["MC_WJetsToQQ_HT400to600", "MC_WJetsToQQ_HT600to800"],
            "MC_ZZ"                 : ["MC_ZZ"],
            "MC_ZJetsToQQ"          : ["MC_ZJetsToQQ_HT800toInf"],
            "MC_QCD_Pt"             : ["MC_QCD_Pt170to300", "MC_QCD_Pt300to470", "MC_QCD_Pt470to600", "MC_QCD_Pt600to800", "MC_QCD_Pt800to1000"],
            # "MC_QCD_Pt"             : ["MC_QCD_Pt15to30", "MC_QCD_Pt30to50", "MC_QCD_Pt50to80", "MC_QCD_Pt80to120", "MC_QCD_Pt120to170", "MC_QCD_Pt170to300", "MC_QCD_Pt300to470", "MC_QCD_Pt470to600", "MC_QCD_Pt600to800", "MC_QCD_Pt800to1000", "MC_QCD_Pt1000to1400", "MC_QCD_Pt1400to1800", "MC_QCD_Pt1800to2400", "MC_QCD_Pt2400to3200", "MC_QCD_Pt3200toInf"],
            "DATA_SingleElectron"   : ["DATA_SingleElectron_RunB", "DATA_SingleElectron_RunC", "DATA_SingleElectron_RunD","DATA_SingleElectron_RunE","DATA_SingleElectron_RunF"],
            "DATA_SingleMuon"       : ["DATA_SingleMuon_RunB", "DATA_SingleMuon_RunC", "DATA_SingleMuon_RunD","DATA_SingleMuon_RunE","DATA_SingleMuon_RunF"],
            }
        Samples_dict_["2018"] = {
            #Matching Level
            "MC_QCD"                : ["MC_QCD_HT50to100", "MC_QCD_HT100to200", "MC_QCD_HT200to300"],
            "MC_Top"                : ["MC_TTToSemiLeptonic"],
            "MC_W"                  : ["MC_WZ", "MC_WW"],
            "MC_Z"                  : ["MC_ZZ", "MC_WZ_Zmatch"],
            #Analysis Level
            # "MC_TTbar"              : ["MC_TTToSemiLeptonic"],
            # "MC_QCD_HT"             : ["MC_QCD_HT50to100", "MC_QCD_HT100to200", "MC_QCD_HT200to300", "MC_QCD_HT300to500", "MC_QCD_HT500to700", "MC_QCD_HT700to1000", "MC_QCD_HT1000to1500", "MC_QCD_HT1500to2000", "MC_QCD_HT2000toInf"],
            "MC_DY"                 : ["MC_DY_HT70to100", "MC_DY_HT100to200", "MC_DY_HT200to400", "MC_DY_HT400to600", "MC_DY_HT600to800", "MC_DY_HT800to1200", "MC_DY_HT1200to2500", "MC_DY_HT2500toInf"],
            "MC_TTbar"              : ["MC_TTbar"],
            "MC_WW_incl"            : ["MC_WW"],
            "MC_WZ_incl"            : ["MC_WZ"],
            "MC_ZZ_incl"            : ["MC_ZZ"],
            "MC_WW"                 : ["MC_WWTo4Q", "MC_WWToLNuQQ", "MC_WWTo2L2Nu"],
            "MC_WZ"                 : ["MC_WZToLNu2Q", "MC_WZTo2Q2Nu", "MC_WZTo2L2Q", "MC_WZTo1L3Nu", "MC_WZTo3LNu"],
            "MC_ZZ"                 : ["MC_ZZTo4Q", "MC_ZZTo2Q2Nu", "MC_ZZTo2L2Q", "MC_ZZTo2L2Nu", "MC_ZZTo4L"],
            "DATA_SingleElectron"   : ["DATA_SingleElectron_RunA","DATA_SingleElectron_RunB", "DATA_SingleElectron_RunC", "DATA_SingleElectron_RunD"],
            "DATA_SingleMuon"       : ["DATA_SingleMuon_RunA","DATA_SingleMuon_RunB", "DATA_SingleMuon_RunC", "DATA_SingleMuon_RunD"],
            }


        self.RunControls = {"2016": ["B", "C", "D", "E", "F", "G", "H"],
                            "2017": ["B", "C", "D", "E", "F"],
                            "2018": ["A", "B", "C", "D"],
                            }
        Samples_dict_ = {
            "MC_DY"                 : ["MC_DY_HT100to200", "MC_DY_HT200to400", "MC_DY_HT400to600", "MC_DY_HT600to800", "MC_DY_HT800to1200", "MC_DY_HT1200to2500", "MC_DY_HT2500toInf"],
            "MC_TTbar"              : ["MC_TTTo2L2Nu", "MC_TTToHadronic", "MC_TTToSemiLeptonic"],
            "MC_WW_incl"            : ["MC_WW"],
            "MC_WZ_incl"            : ["MC_WZ"],
            "MC_ZZ_incl"            : ["MC_ZZ"],
            "MC_WW"                 : ["MC_WWTo4Q", "MC_WWToLNuQQ", "MC_WWTo2L2Nu"],
            "MC_WZ"                 : ["MC_WZToLNu2Q", "MC_WZTo2Q2Nu", "MC_WZTo2L2Q", "MC_WZTo1L3Nu", "MC_WZTo3LNu"],
            "MC_ZZ"                 : ["MC_ZZTo4Q", "MC_ZZTo2Q2Nu", "MC_ZZTo2L2Q", "MC_ZZTo2L2Nu", "MC_ZZTo4L"],
            "DATA_SingleElectron"   : ["DATA_SingleElectron_RunA", "DATA_SingleElectron_RunB", "DATA_SingleElectron_RunC", "DATA_SingleElectron_RunD", "DATA_SingleElectron_RunE", "DATA_SingleElectron_RunF", "DATA_SingleElectron_RunG", "DATA_SingleElectron_RunH" ],
            "DATA_SingleMuon"       : ["DATA_SingleMuon_RunA", "DATA_SingleMuon_RunB", "DATA_SingleMuon_RunC", "DATA_SingleMuon_RunD", "DATA_SingleMuon_RunE", "DATA_SingleMuon_RunF", "DATA_SingleMuon_RunG", "DATA_SingleMuon_RunH" ],
            self.signal             : self.SignalSamples,
            }

        self.Samples_dict = {}
        self.Samples_original = {}
        self.Samples_Category = {}

        for year in self.RunControls:
            for x in Samples_dict_:
                loop_over = Samples_dict_[x]
                if (year=="2016" and x=="MC_TTbar"): loop_over = ["MC_TTbar"]
                if (year!="2016" and (x=="MC_WW" or x=="MC_WZ" or x=="MC_ZZ") ): loop_over = Samples_dict_[x+"_incl"]
                self.Samples_dict.setdefault(year, {}).setdefault(x+"_"+year, [el+"_"+year for el in loop_over] )
            for x in ["DATA_SingleElectron","DATA_SingleMuon"]:
                self.Samples_dict[year][x+"_"+year] = [el for el in self.Samples_dict[year][x+"_"+year] if any("Run"+crl in el for crl in self.RunControls[year])]
            self.Samples_original.setdefault(year, list(dict.fromkeys([el for list_ in self.Samples_dict[year].values() for el in list_])))
            self.Samples_Category.setdefault(year, self.Samples_dict[year].keys())
            self.Samples_Category[year].remove(self.signal+"_"+year)
            self.Samples_Category[year].extend(self.Samples_dict[year][self.signal+"_"+year])

        self.Samples_NamesAll = []
        for x in self.Samples_original:
            self.Samples_NamesAll.extend(self.Samples_original[x])

        self.Samples_CategoryAll = []
        for x in self.Samples_Category:
            self.Samples_CategoryAll.extend(self.Samples_Category[x])

        Samples_matching2016  = ["MC_QCD"]
        Samples_matching2017  = ["MC_HWW", "MC_QCD", "MC_Top", "MC_W", "MC_Z"]
        Samples_matching2018  = ["MC_QCD", "MC_Top", "MC_W", "MC_Z"]

        self.Samples_matching = Samples_matching2018 if self.year=="2018" else Samples_matching2017 if self.year=="2017" else Samples_matching2016
