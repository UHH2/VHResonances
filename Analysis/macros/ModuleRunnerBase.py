import os
import sys
import itertools


class GenericPath:
    ''' Class container for paths '''
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
    ''' Class container for list of objects '''
    def __init__(self):
        GenericPath.__init__(self)
        self.PrefixrootFile     = "uhh2.AnalysisModuleRunner."
        self.Collections        = ["Puppi"]
        self.Channels           = ["muon", "electron", "invisible"]
        self.Systematics        = ["nominal", "JER_up", "JER_down", "JEC_up", "JEC_down", "MuonScale_up", "MuonScale_down"]
        self.Systematics_Scale  = ["PU_up", "PU_down"]
        self.Signal             = "MC_ZprimeToZH"
        self.MainBkg            = "MC_DY"
        self.MassPoints         = [600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 7000, 8000]
        self.SignalSamples      = [self.Signal+mode+"_M"+str(mass) for mass in self.MassPoints for mode in ["","_inv"]]
        self.RunPeriods_Dict    = {"2016": ["B", "C", "D", "E", "F", "G", "H"],
                                   "2017": ["B", "C", "D", "E", "F"],
                                   "2018": ["A", "B", "C", "D"]
                                   }
        self.years              = sorted(self.RunPeriods_Dict.keys())
        self.AllRunPeriods      = list(set(itertools.chain.from_iterable(self.RunPeriods_Dict.values())))

        # TODO MC_DY_HT70to100
        self.Generic_SubSamples_Dict_ = {
            "MC_DY"                 : [proc+subsample for proc in ["MC_DY_HT", "MC_DY_inv_HT"] for subsample in ["100to200", "200to400", "400to600", "600to800", "800to1200", "1200to2500", "2500toInf",]] + ["MC_DY_inv_PtZ_"+subsample for subsample in ["50to100", "100to250", "250to400", "400to650", "650toInf"]],
            "MC_TTbar"              : ["MC_TTTo2L2Nu", "MC_TTToHadronic", "MC_TTToSemiLeptonic"],
            "MC_WW_incl"            : ["MC_WW"],
            "MC_WZ_incl"            : ["MC_WZ"],
            "MC_ZZ_incl"            : ["MC_ZZ"],
            "MC_WW"                 : ["MC_WWTo4Q", "MC_WWToLNuQQ", "MC_WWTo2L2Nu"],
            "MC_WZ"                 : ["MC_WZToLNu2Q", "MC_WZTo2Q2Nu", "MC_WZTo2L2Q", "MC_WZTo1L3Nu", "MC_WZTo3LNu"],
            "MC_ZZ"                 : ["MC_ZZTo4Q", "MC_ZZTo2Q2Nu", "MC_ZZTo2L2Q", "MC_ZZTo2L2Nu", "MC_ZZTo4L"],
            "DATA_SingleElectron"   : ["DATA_SingleElectron_Run"+str(run) for run in self.AllRunPeriods],
            "DATA_SingleMuon"       : ["DATA_SingleMuon_Run"+str(run) for run in self.AllRunPeriods],
            "DATA_MET"              : ["DATA_MET_Run"+str(run) for run in self.AllRunPeriods],
            self.Signal             : self.SignalSamples,
            }

        self.ExtractVariableFromConstants()
        self.defineSamples()

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

    def defineSamples(self):
        '''
        Create Useful list of subsample and processes
        Example Process = MC_DY, MC_TTbar, DATA_MET, ...
        Example SubSample = MC_DY_HT100to200, MC_TTToHadronic, DATA_MET_RunB, ...
        '''

        self.Samples_Year_Dict = {}
        self.SubSamples_Year_Dict = {}
        self.Processes_Year_Dict = {}

        ''' This part is a bit arbitrary. The idea is to catch all the combinations for different years'''
        for year in self.years:
            for subsample in sorted(self.Generic_SubSamples_Dict_):
                loop_over = self.Generic_SubSamples_Dict_[subsample]
                if (year=="2016" and subsample=="MC_DY"): loop_over = filter(lambda subsample: not "inv" in subsample or "PtZ" in subsample ,loop_over)
                if (year!="2016" and subsample=="MC_DY"): loop_over = filter(lambda subsample: not "PtZ" in subsample ,loop_over)
                if (year=="2016" and subsample=="MC_TTbar"): loop_over = ["MC_TTbar"]
                if (year!="2016" and (subsample=="MC_WW" or subsample=="MC_WZ" or subsample=="MC_ZZ") ): loop_over = self.Generic_SubSamples_Dict_[subsample+"_incl"]
                if "DATA" in subsample: loop_over = [subsample+"_Run"+str(run) for run in self.RunPeriods_Dict[year]]
                self.Samples_Year_Dict.setdefault(year, {}).setdefault(subsample+"_"+year, [el+"_"+year for el in sorted(loop_over)] )
                self.Processes_Year_Dict.setdefault(year, []).append(subsample+"_"+year)
            self.Processes_Year_Dict[year].remove(self.Signal+"_"+year)
            self.Processes_Year_Dict[year].extend(self.Samples_Year_Dict[year][self.Signal+"_"+year])
            self.SubSamples_Year_Dict.setdefault(year, sorted(list(dict.fromkeys([el for list_ in self.Samples_Year_Dict[year].values() for el in list_]))))

        # List of all sub sample for all years and all processes
        self.AllSubSamples_List = sorted(list(set(itertools.chain.from_iterable(self.SubSamples_Year_Dict.values()))))

        # List of all processes for all years
        self.AllProcesses_List = sorted(list(set(itertools.chain.from_iterable(self.Processes_Year_Dict.values()))))



class ModuleRunnerBase(VariablesBase):
    ''' Class container for list of objects for particular year '''
    def __init__(self,year="2016"):
        VariablesBase.__init__(self)
        self.year = year
        self.defineDirectories()
        self.lumi_fb  = round(float(self.lumi_map[self.year]["lumi_fb"]),1)
        self.lumi_pb  = int(self.lumi_map[self.year]["lumi_pb"])
        self.lumi_sys = round(float(self.lumi_map[self.year]["uncertainty"]),1)
        self.Samples_Dict    = self.Samples_Year_Dict[self.year]
        self.SubSamples_Dict = self.SubSamples_Year_Dict[self.year]
        self.Processes_Dict  = self.Processes_Year_Dict[self.year]

    def defineDirectories(self):
        self.Path_SFRAME        = self.Path_SFRAME+self.year+"/"
        self.ConfigDir          = self.Path_ANALYSIS+"config/"
        self.SubmitDir          = self.ConfigDir+"SubmittedJobs/"+self.year+"/"
