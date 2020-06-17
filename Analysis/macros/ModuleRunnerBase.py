import os
import sys

class GenericPath:
    def __init__(self):
        self.user = os.environ["USER"]
        self.cmssw_base = os.environ["CMSSW_BASE"]
        self.isMaxwell = True if "beegfs" in os.getcwd() else False
        self.Path_NFS       = "/nfs/dust/cms/user/"+self.user+"/"
        self.Path_Maxwell   = "/beegfs/desy/user/"+os.environ["USER"]+"/"
        self.Path_ANALYSIS  = self.cmssw_base+"/src/UHH2/VHResonances/"
        self.Path_SFRAME    = self.Path_NFS+"sframe_all/"
        self.Path_STORAGE   = self.Path_NFS+"WorkingArea/File/Analysis/"
        self.Path_SPlotter  = self.Path_NFS+"WorkingArea/SFramePlotter/"
    def Get(self, source):
        return getattr(self,source)

class VariablesBase(GenericPath):
    def __init__(self):
        GenericPath.__init__(self)
        self.PrefixrootFile     = "uhh2.AnalysisModuleRunner."
        self.signal             = "MC_ZprimeToZH"
        self.Collections        = ["Puppi"]
        self.Channels           = ["muon", "electron", "invisible"]
        self.Systematics        = ["nominal", "JER_up", "JER_down", "JEC_up", "JEC_down"]
        self.Systematics_Scale  = ["PU_up", "PU_down"]
        self.years              = ["2016", "2017", "2018"]

        self.SignalSamples  = ["MC_ZprimeToZH_M600", "MC_ZprimeToZH_M800", "MC_ZprimeToZH_M1000", "MC_ZprimeToZH_M1200", "MC_ZprimeToZH_M1400", "MC_ZprimeToZH_M1600", "MC_ZprimeToZH_M1800", "MC_ZprimeToZH_M2000", "MC_ZprimeToZH_M2500", "MC_ZprimeToZH_M3000", "MC_ZprimeToZH_M3500", "MC_ZprimeToZH_M4000", "MC_ZprimeToZH_M4500", "MC_ZprimeToZH_M5000", "MC_ZprimeToZH_M5500", "MC_ZprimeToZH_M6000", "MC_ZprimeToZH_M7000", "MC_ZprimeToZH_M8000" ]

        self.RunControls    = { "2016": ["B", "C", "D", "E", "F", "G", "H"],
                                "2017": ["B", "C", "D", "E", "F"],
                                "2018": ["A", "B", "C", "D"],
                                }

        self.Samples_dict_ = {# TODO MC_DY_HT70to100
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

        self.ExtractVariableFromConstants()

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

class ModuleRunnerBase(VariablesBase):
    def __init__(self,year="2016"):
        VariablesBase.__init__(self)
        self.year = year
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
    def defineSamples(self):

        self.Samples_dict = {}
        self.Samples_original = {}
        self.Samples_Category = {}

        for year in self.RunControls:
            for x in self.Samples_dict_:
                loop_over = self.Samples_dict_[x]
                if (year=="2016" and x=="MC_TTbar"): loop_over = ["MC_TTbar"]
                if (year!="2016" and (x=="MC_WW" or x=="MC_WZ" or x=="MC_ZZ") ): loop_over = self.Samples_dict_[x+"_incl"]
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
