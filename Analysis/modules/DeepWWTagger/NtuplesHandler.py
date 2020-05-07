import os
import os.path
import sys
from glob import glob
import time
import math
import copy

from sklearn.utils import shuffle
import numpy as np
from root_numpy import root2array, rec2array, fill_hist
import ROOT

analysispath = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v2/CMSSW_10_2_16/src/UHH2/VHResonances/" if "/nfs" in os.getcwd() else "/beegfs/desy/user/amalara//VHResonances/"
sys.path.append(analysispath+"Analysis/macros/")
from ModuleRunnerBase import *
sys.path.append(GenericPath().Path_ANALYSIS+"Analysis/python")

ROOT.gInterpreter.ProcessLine('#include "'+analysispath+'include/constants.hpp"')

from tdrstyle_all import *
from parallelise import *



# ##     ##    ###    ######## ########  #### ##     ##     ######   #######  ########  ########
# ###   ###   ## ##      ##    ##     ##  ##   ##   ##     ##    ## ##     ## ##     ##    ##
# #### ####  ##   ##     ##    ##     ##  ##    ## ##      ##       ##     ## ##     ##    ##
# ## ### ## ##     ##    ##    ########   ##     ###        ######  ##     ## ########     ##
# ##     ## #########    ##    ##   ##    ##    ## ##            ## ##     ## ##   ##      ##
# ##     ## ##     ##    ##    ##    ##   ##   ##   ##     ##    ## ##     ## ##    ##     ##
# ##     ## ##     ##    ##    ##     ## #### ##     ##     ######   #######  ##     ##    ##




def MergeSort(matrix,pt):
    RecursiveMergeSort(matrix, pt)
    return matrix

def RecursiveMergeSort(matrix, pt):
    if matrix.shape[1]>1:
        mid = matrix.shape[1]//2
        lefthalf  = copy.deepcopy(matrix[:, :mid])
        righthalf = copy.deepcopy(matrix[:, mid:])
        RecursiveMergeSort(lefthalf, pt)
        RecursiveMergeSort(righthalf, pt)
        i=0; j=0; k=0;
        while i < lefthalf.shape[1] and j < righthalf.shape[1]:
            if lefthalf[pt][i] > righthalf[pt][j]:
                matrix[:,k]=lefthalf[:,i]
                i=i+1
            else:
                matrix[:,k]=righthalf[:,j]
                j=j+1
            k=k+1
        while i < lefthalf.shape[1]:
            matrix[:,k]=lefthalf[:,i]
            i=i+1; k=k+1;
        while j < righthalf.shape[1]:
            matrix[:,k]=righthalf[:,j]
            j=j+1; k=k+1;


#  ######   #######  ##    ## ##     ## ######## ########  ########    ##    ## ######## ##     ## ########  ##       ########  ######
# ##    ## ##     ## ###   ## ##     ## ##       ##     ##    ##       ###   ##    ##    ##     ## ##     ## ##       ##       ##    ##
# ##       ##     ## ####  ## ##     ## ##       ##     ##    ##       ####  ##    ##    ##     ## ##     ## ##       ##       ##
# ##       ##     ## ## ## ## ##     ## ######   ########     ##       ## ## ##    ##    ##     ## ########  ##       ######    ######
# ##       ##     ## ##  ####  ##   ##  ##       ##   ##      ##       ##  ####    ##    ##     ## ##        ##       ##             ##
# ##    ## ##     ## ##   ###   ## ##   ##       ##    ##     ##       ##   ###    ##    ##     ## ##        ##       ##       ##    ##
#  ######   #######  ##    ##    ###    ######## ##     ##    ##       ##    ##    ##     #######  ##        ######## ########  ######



class NtuplesHandlerBase(ModuleRunnerBase):
    def __init__(self, Channel="", Collection="", Sample="", SubSample=None, year="2017", extraText="", isTest=False):
        ModuleRunnerBase.__init__(self,year)
        if "beegfs" in os.getcwd(): self.SetForMaxwell()
        self.TreeName = "AnalysisTree"
        self.Module = "GenericCleaning"
        self.Mode = "MC"
        self.isTest = isTest
        self.Channel = Channel
        self.Collection = Collection
        self.Sample = Sample
        self.SubSample = SubSample
        self.FileExt = "test_"+self.SubSample+"_index.npy" if self.isTest else "_"+self.SubSample+"_index.npy"
        self.extraText = extraText
        # self.FileStorageInput = self.Path_STORAGE+self.year+"/"+self.Module+"/"+self.Collection+"/"+self.Channel+"channel"+"/"
        self.FileStorageInput = self.Path_STORAGE+self.year+"/"+self.Module+"/All/"+self.Channel+"channel"+"/"
        self.FileStorageOutput = self.Path_STORAGE.replace("Analysis","NeuralNetwork")+"input_varariables/VHResonances/"+self.year+"/"+self.Collection+"/"+self.Channel+"channel"+"/"
        self.Workdir = self.FileStorageInput+"workdir_"+self.Module+"_"+self.SubSample+"/"
        self.outdir = self.FileStorageOutput+Sample+"/Vars/"
        self.n_jet = 3
        self.n_PF = 100
        self.FS = np.float32
        self.MatchingTag = 33 #extract manually from TopJet header
        self.MatchingStatusTag = 34 #extract manually from TopJet header
        self.filesize=100000
        self.MakeVars()
    def SetForMaxwell(self):
        self.Path_STORAGE = self.Path_Maxwell+"WorkingArea/File/"
    def MakeVars(self):
        self.xsec = { "MC_QCD_Pt170to300": 0.003917790623,    "MC_QCD_Pt300to470": 0.0001188208272, "MC_QCD_Pt470to600": 0.00001999487323,
                     "MC_QCD_Pt600to800": 0.000001217875878, "MC_QCD_Pt800to1000": 0.0000003362947067,}

        self.MatchingDict   = {"MC_HWW"                 : {"MatchingName":"SubDaughtersMatched", "Matching": ROOT.HWWMatch, "MatchingStatus": ROOT.Hadronic},
                               "MC_Hbb"                 : {"MatchingName":"DaughterMatched",     "Matching": ROOT.HbbMatch, "MatchingStatus": ROOT.DaughterMatched},
                               "MC_ZH_HToBB_ZToQQ"      : {"MatchingName":"DaughterMatched",     "Matching": ROOT.HbbMatch, "MatchingStatus": ROOT.DaughterMatched},
                               "MC_ZH_HToBB_ZToLL"      : {"MatchingName":"DaughterMatched",     "Matching": ROOT.HbbMatch, "MatchingStatus": ROOT.DaughterMatched},
                               "MC_QCD"                 : {"MatchingName":"MotherMatched",       "Matching": ROOT.qMatch,   "MatchingStatus": ROOT.MotherMatched},
                               "MC_DY1JetsToLL"         : {"MatchingName":"MotherMatched",       "Matching": ROOT.qMatch,   "MatchingStatus": ROOT.MotherMatched},
                               "MC_DY2JetsToLL"         : {"MatchingName":"MotherMatched",       "Matching": ROOT.qMatch,   "MatchingStatus": ROOT.MotherMatched},
                               "MC_DY3JetsToLL"         : {"MatchingName":"MotherMatched",       "Matching": ROOT.qMatch,   "MatchingStatus": ROOT.MotherMatched},
                               "MC_DY4JetsToLL"         : {"MatchingName":"MotherMatched",       "Matching": ROOT.qMatch,   "MatchingStatus": ROOT.MotherMatched},
                               "MC_Top"                 : {"MatchingName":"SubDaughtersMatched", "Matching": ROOT.tWbMatch, "MatchingStatus": ROOT.Hadronic3},
                               "MC_TTTo2L2Nu"           : {"MatchingName":"SubDaughtersMatched", "Matching": ROOT.tWbMatch, "MatchingStatus": ROOT.Hadronic3},
                               "MC_TTbarSemiLeptonic"   : {"MatchingName":"SubDaughtersMatched", "Matching": ROOT.tWbMatch, "MatchingStatus": ROOT.Hadronic3},
                               "MC_TTToHadronic"        : {"MatchingName":"SubDaughtersMatched", "Matching": ROOT.tWbMatch, "MatchingStatus": ROOT.Hadronic3},
                               "MC_W"                   : {"MatchingName":"DaughterMatched",     "Matching": ROOT.WqqMatch, "MatchingStatus": ROOT.DaughterMatched},
                               "MC_Z"                   : {"MatchingName":"DaughterMatched",     "Matching": ROOT.ZqqMatch, "MatchingStatus": ROOT.DaughterMatched},
                               }
        self.matchingNames = {"NotMatched":[ROOT.NotMatched],"MotherMatched":[ROOT.MotherMatched],"DaughterMatched":[ROOT.DaughterMatched],"SubDaughtersMatched":[ROOT.Hadronic,ROOT.Hadronic1,ROOT.Hadronic2,ROOT.Hadronic3]}
        self.RootObjects = {"Event"         : "",
                            "TopJet"        : "hotvrPuppi.m_.m_" if self.Collection=="HOTVR" else "jetsAk8CHSSubstructure_SoftDropCHS.m_" if self.Collection=="CHS" else "jetsAk8PuppiSubstructure_SoftDropPuppi.m_",
                            "PFParticles"   : "PFParticles.m_",
                            }
                            # "Electron" : "slimmedElectronsUSER.m_",
                            # "Muon"     : "slimmedMuonsUSER.m_",
                            # "Jet"      : "jetsAk4Puppi.m_" if self.Collection=="CHS" else "jetsAk4CHS.m_"}
        self.RootVarNameListsDict = {}
        self.RootVarNameListsDict["Event"] = ["event", "weight_GLP", "weight_lumi", "weight_pu", "weight_pu_down", "weight_pu_up"]
        self.RootVarNameListsDict["Jet"] = ["pt", "eta", "phi", "energy","numberOfDaughters", "jetArea"]
        # "btag_combinedSecondaryVertex", "btag_combinedSecondaryVertexMVA", "btag_DeepCSV_probb", "btag_DeepCSV_probbb",
        # "btag_DeepFlavour_probbb", "btag_DeepFlavour_probb", "btag_DeepFlavour_problepb", "btag_DeepFlavour_probuds", "btag_DeepFlavour_probc", "btag_DeepFlavour_probg",
        # "neutralEmEnergyFraction", "neutralHadronEnergyFraction", "chargedEmEnergyFraction", "chargedHadronEnergyFraction", "muonEnergyFraction", "photonEnergyFraction",
        # "chargedMultiplicity", "neutralMultiplicity", "muonMultiplicity", "electronMultiplicity", "photonMultiplicity",
        # "puppiMultiplicity", "neutralPuppiMultiplicity", "neutralHadronPuppiMultiplicity", "photonPuppiMultiplicity", "HFHadronPuppiMultiplicity", "HFEMPuppiMultiplicity",]
        self.RootVarNameListsDict["TopJet"] = self.RootVarNameListsDict["Jet"] + ["tau1", "tau2", "tau3", "tau4", "softdropmass", "tags.first","tags.second",]
            # "mvahiggsdiscr", "btag_BoostedDoubleSecondaryVertexAK8", "tau1_groomed", "tau2_groomed", "tau3_groomed", "tau4_groomed",
            # "btag_MassDecorrelatedDeepBoosted_probHbb", "btag_MassDecorrelatedDeepBoosted_probHcc", "btag_MassDecorrelatedDeepBoosted_probHqqqq",
            # "btag_MassDecorrelatedDeepBoosted_probQCDc", "btag_MassDecorrelatedDeepBoosted_probQCDbb", "btag_MassDecorrelatedDeepBoosted_probQCDothers", "btag_MassDecorrelatedDeepBoosted_probQCDb", "btag_MassDecorrelatedDeepBoosted_probQCDcc",
            # "btag_MassDecorrelatedDeepBoosted_probTbqq", "btag_MassDecorrelatedDeepBoosted_probTbcq", "btag_MassDecorrelatedDeepBoosted_probTbq","btag_MassDecorrelatedDeepBoosted_probTbc",
            # "btag_MassDecorrelatedDeepBoosted_probWqq", "btag_MassDecorrelatedDeepBoosted_probWcq","btag_MassDecorrelatedDeepBoosted_probZcc", "btag_MassDecorrelatedDeepBoosted_probZqq","btag_MassDecorrelatedDeepBoosted_probZbb",]
        self.RootVarNameListsDict["Electron"] = ["pt", "eta", "phi", "energy", "charge"]
        self.RootVarNameListsDict["Muon"] = self.RootVarNameListsDict["Electron"]
        self.RootVarNameListsDict["PFParticles"] = self.RootVarNameListsDict["Electron"]+["particleID", "puppi_weight"]
        self.RootVarNamesDict = {}
        for Object in self.RootObjects:
            self.RootVarNamesDict[Object] = [self.RootObjects[Object]+var for var in self.RootVarNameListsDict[Object]]
            self.RootVarNamesDict[Object] = list(map(lambda x: x.replace(".m_", ".") if "tags" in x else x, self.RootVarNamesDict[Object]))
            self.RootVarNamesDict[Object] = list(map(lambda x: x+"_groomed" if ("tau" in x and self.Collection=="HOTVR") else x, self.RootVarNamesDict[Object]))
    @timeit
    def ConvertRoot2Numpy(self,Sample,firstfile=None,lastfile=None):
        list_ = sorted(glob(self.Workdir+"/*root"))
        if self.isTest: list_= list_[:10]
        list_= list_[firstfile:lastfile]
        for y in list_:
            print y
            try:
                self.Root2Numpy(filename=y, Sample=Sample)
            except Exception as e:
                    print "ERROR "+str(e) if not "none of the input files contain the requested tree" in str(e) else ""
    @timeit
    def Root2Numpy(self, filename, Sample=None):
        '''
        The Variables loaded have a shape (#events,#variables), i.e. len(self.RootVarNamesDict["TopJet"]).
        Each variable can be a number or an array. Usually event variables are numbers.
        There is a check in case this is not true to take only the first element.
        Other variables are store as array and need a bit of manipulation.
        Each array correspond to the info of all jets.
        '''
        start = 0 if self.isTest else None
        stop = 100 if self.isTest else None
        a = os.system("mkdir -p "+self.outdir)
        index = filename[filename.rfind('_')+1:-len(".root")]
        ext = self.FileExt.replace("index",str(index))
        if not Sample: Sample = self.Sample

        ''' Load Variables '''
        tag1      = self.RootVarNamesDict["TopJet"].index(self.RootObjects["TopJet"].replace(".m_", ".")+"tags.first")
        tag2      = self.RootVarNamesDict["TopJet"].index(self.RootObjects["TopJet"].replace(".m_", ".")+"tags.second")
        pt_index  = self.RootVarNameListsDict["TopJet"].index("pt")

        #TODO This step is by far the slowest one. Try to reduce the time.
        VarsEvent     = rec2array(root2array(filenames=filename, treename=self.TreeName, branches=self.RootVarNamesDict["Event"], start=start, stop=stop))
        VarsTopJet    = rec2array(root2array(filenames=filename, treename=self.TreeName, branches=self.RootVarNamesDict["TopJet"], start=start, stop=stop))
        VarsPF        = rec2array(root2array(filenames=filename, treename=self.TreeName, branches=self.RootVarNamesDict["PFParticles"], start=start, stop=stop))
        VarsPF_indexs = rec2array(root2array(filenames=filename,treename=self.TreeName,branches=[self.RootObjects["TopJet"]+x for x in ["numberOfDaughters","pfcand_indexs"]], start=start, stop=stop))

        ''' Check if no events are loaded '''
        if len(VarsEvent)==0: return None
        ''' Check if no jets are loaded '''
        if len(VarsTopJet[0][0])==0: return None

        ''' Reshape Event variables '''
        for col, var in enumerate(VarsEvent[0]):
            if (isinstance(var,np.ndarray)):
                VarsEvent[:,col] = np.array(map(lambda x: x[0], VarsEvent[:,col]))

        ''' Loop over events to construct the variable shape '''
        VarTopJetList=[]
        VarPFList=[]
        EventToRemoveList=[]
        for ev in range(len(VarsTopJet)):
            NumOfDaug = [0]+list(VarsPF_indexs[ev][0])
            PF_index  =     list(VarsPF_indexs[ev][1])
            if len(PF_index)==0:
                EventToRemoveList.append(ev)
                continue

            '''
            Now each element of the list correspond to 1 variable.
            Going to reduce unused tag variables. Selecting only matching and matching status.
            JET will be an array (#variables,#jets).
            '''
            JET = list(map(lambda x: x.reshape(1,x.shape[0]).astype(self.FS), VarsTopJet[ev]))
            MatchingMask = JET[tag1][0]== self.MatchingTag
            MatchingStatusMask = JET[tag1][0]== self.MatchingStatusTag
            JET[tag1] = JET[tag2][:,MatchingMask]
            JET[tag2] = JET[tag2][:,MatchingStatusMask]
            JET = np.concatenate(JET)

            """
            Select up to n_jet. If less, set all var to 0
            Select up to n_PF
            """
            NumOfDaug = list(map(lambda x: np.array(NumOfDaug[:x+1]).sum(), range(len(NumOfDaug))))
            PF_index  = list(map(lambda x: PF_index[NumOfDaug[x]:NumOfDaug[x+1]], range(len(NumOfDaug)-1)))
            PF = np.swapaxes(np.concatenate(list(map(lambda x: x.reshape(1,x.shape[0]).astype(self.FS), VarsPF[ev]))), 0, 1)
            PF = np.concatenate(list(map(lambda x:  np.expand_dims(np.swapaxes(MergeSort(np.swapaxes(PF[x], 0, 1),pt_index), 0, 1)[:self.n_PF],axis=0) if PF[x].shape[0] >= self.n_PF else
                                         np.expand_dims(np.vstack((np.swapaxes(MergeSort(np.swapaxes(PF[x], 0, 1),pt_index), 0, 1), np.zeros((self.n_PF-PF[x].shape[0], PF[x].shape[1])))),axis=0), PF_index)))
            PF = np.swapaxes(np.swapaxes(PF,0,1),1,2)
            VarTopJetList += [np.expand_dims(JET[:, :self.n_jet],axis=0)] if JET.shape[1] >= self.n_jet else [np.expand_dims(np.hstack((JET, np.zeros((JET.shape[0],self.n_jet-JET.shape[1])))),axis=0)]
            VarPFList     += [np.expand_dims(PF[:,:,:self.n_jet],axis=0)] if JET.shape[1] >= self.n_jet else [np.expand_dims(np.concatenate((PF,  np.zeros((PF.shape[0],PF.shape[1],self.n_jet-JET.shape[1]))),axis=2),axis=0)]
        VarsTopJet    = np.concatenate(VarTopJetList,axis=0)
        VarsPF  = np.concatenate(VarPFList,axis=0)
        """ Remove event if something is wrong with jets """
        for ev in reversed(EventToRemoveList):
            VarsEvent = np.delete(VarsEvent, ev, 0)
        """
        Concatenate all Variables in 1 3D array.
        Each row correspond to 1 jet. Event variables are added to jet info.
        """
        VarsTopJet_list = []
        VarsPF_list = []
        for i in range(self.n_jet):
            VarsTopJet_list.append(np.concatenate((VarsTopJet[:,:,i], VarsEvent[:,:]),axis=1))
            VarsPF_list.append(VarsPF[:,:,:,i])
        VarsTopJet = np.concatenate(VarsTopJet_list)
        VarsPF = np.concatenate(VarsPF_list)
        VarsTopJet, VarsPF = shuffle(VarsTopJet, VarsPF, random_state=0)
        print "Shape TopJet:", VarsTopJet.shape
        print "Shape PFCand:", VarsPF.shape
        if "QCD" in filename:
            for x in self.xsec:
                if x in filename: maxjets = max(1,int(len(VarsTopJet)*self.xsec[x]/self.xsec["MC_QCD_Pt300to470"]))
            VarsTopJet = VarsTopJet[:maxjets]
            VarsPF = VarsPF[:maxjets]
            print "It's QCD. Taken only", maxjets, self.xsec[x]/self.xsec["MC_QCD_Pt300to470"]*100,"%"
            print "Shape TopJet:", VarsTopJet.shape
            print "Shape PFCand:", VarsPF.shape
        mask  = VarsTopJet[:,tag1] == self.MatchingDict[self.Sample]["Matching"]
        mask *= VarsTopJet[:,tag2] == self.MatchingDict[self.Sample]["MatchingStatus"]
        print "Shape TopJet_others:", VarsTopJet[~mask].shape
        print "Shape PFCand_others:", VarsPF[~mask].shape
        print "Shape TopJet_match :", VarsTopJet[mask].shape
        print "Shape PFCand_match :", VarsPF[mask].shape
        for others in ["","_others"]:
            VarsTopJet_masked = VarsTopJet[mask] if others=="" else VarsTopJet[~mask]
            VarsPF_masked     = VarsPF[mask] if others=="" else VarsPF[~mask]
            if (len(VarsTopJet_masked)!=len(VarsPF_masked)): raise RuntimeError("Number of events don't match.")
            nElements = len(VarsTopJet_masked)
            subindex = 0
            print "save as", self.outdir+"/TopJet"+others+ext.replace(".npy","_"+str(subindex)+".npy")
            while nElements >= self.filesize:
                np.save(self.outdir+"/TopJet"+others+ext.replace(".npy","_"+str(subindex)+".npy"), VarsTopJet_masked[:self.filesize].astype(self.FS))
                np.save(self.outdir+"/PFParticles"+others+ext.replace(".npy","_"+str(subindex)+".npy"), VarsPF_masked[:self.filesize].astype(self.FS))
                VarsTopJet_masked = VarsTopJet_masked[self.filesize:]
                VarsPF_masked = VarsPF_masked[self.filesize:]
                nElements = len(VarsTopJet_masked)
                subindex+=1
            if nElements!=0:
                np.save(self.outdir+"/TopJet"+others+ext.replace(".npy","_"+str(subindex)+".npy"), VarsTopJet_masked.astype(self.FS))
                np.save(self.outdir+"/PFParticles"+others+ext.replace(".npy","_"+str(subindex)+".npy"), VarsPF_masked.astype(self.FS))

    def LoadRootObjects(self):
        self.Vars = {}
        ext = self.FileExt.replace("index","*")
        for Object in self.RootObjects:
            self.Vars[Object] = []
            for x in glob(self.outdir+Object+ext):
                self.Vars[Object].append(np.load(x).astype(self.FS))
            self.Vars[Object] = np.concatenate(self.Vars[Object])
            print "\tLoad Vars", Object, self.Vars[Object].shape





class NtuplesHandler(NtuplesHandlerBase):
    def __init__(self, Channel="lepton", Collection="CHS", Samples=[], year="2017", extraText="", isTest=False):
        # self.SamplesNames   = {"MC_HWW":0, "MC_Hbb":1, "MC_Top":2, "MC_QCD":3, "MC_W":4, "MC_Z":5}
        # self.SamplesNames   = {"MC_HWW":0, "MC_QCD":1}
        self.SamplesNames   = {"MC_HWW":0, "MC_QCD":1, "MC_Top":2, "MC_W":3, "MC_Z":4}
        self.colors         = {"MC_HWW":                ROOT.kRed+1,
                               "MC_Hbb":                ROOT.kOrange,
                               "MC_ZH_HToBB_ZToQQ":     ROOT.kYellow,
                               "MC_ZH_HToBB_ZToLL":     ROOT.kOrange+1,
                               "MC_QCD":                ROOT.kGreen-2,
                               "MC_DY1JetsToLL":        ROOT.kGreen-4,
                               "MC_DY2JetsToLL":        ROOT.kGreen+3,
                               "MC_DY3JetsToLL":        ROOT.kGreen+4,
                               "MC_DY4JetsToLL":        ROOT.kSpring+10,
                               "MC_Top":                ROOT.kBlue,
                               "MC_TTTo2L2Nu":          ROOT.kBlue+3,
                               "MC_TTbarSemiLeptonic":  ROOT.kBlue-7,
                               "MC_TTToHadronic":       ROOT.kBlue-9,
                               "MC_W":                  ROOT.kAzure+10,
                               "MC_Z":                  ROOT.kViolet-3,}
        self.SampleRef = "MC_HWW"
        NtuplesHandlerBase.__init__(self,Sample=self.SampleRef,SubSample="MC_HZ_HiggsToWWZToLL", Channel=Channel,Collection=Collection,extraText=extraText,isTest=isTest,year=year)
        self.pathPlots = self.FileStorageOutput+"plots/"
        self.VarNames = ["TopJet"+x for x in self.RootVarNameListsDict["TopJet"]] + self.RootVarNameListsDict["Event"]
        # self.SubSample = SubSample if SubSample else self.Samples_dict[self.Sample]
        if not os.path.exists(self.pathPlots):
            os.makedirs(self.pathPlots)
        if len(Samples)!=0:
            print "Original Sample", self.SamplesNames
            #TODO is it really working? Does what you want?
            self.ReduceSamples(Samples)
            print "ReducedSamples", self.SamplesNames
        self.NHBDict = {}
        for Sample in self.SamplesNames:
            self.NHBDict[Sample] = {}
            for SubSample in self.Samples_dict[Sample]:
                self.NHBDict[Sample][SubSample] = NtuplesHandlerBase(Sample=Sample, SubSample=SubSample, Channel=Channel,Collection=Collection,extraText=extraText,isTest=isTest,year=year)
    def ReduceSamples(self,Samples):
        temp_dict = {}
        for Sample in Samples:
            temp_dict[Sample] = self.SamplesNames[Sample]
        self.SamplesNames = temp_dict
    def LoadVars(self):
        for Sample in self.SamplesNames:
            print "LoadVars", Sample
            for Object in self.RootObjects:
                self.NHBDict[Sample][Object] = []
            for SubSample in self.Samples_dict[Sample]:
                print "\t", SubSample
                self.NHBDict[Sample][SubSample].LoadRootObjects()
                for Object in self.RootObjects:
                    self.NHBDict[Sample][Object].append(self.NHBDict[Sample][SubSample].Vars[Object])
            for Object in self.RootObjects:
                self.NHBDict[Sample][Object] = np.concatenate(self.NHBDict[Sample][Object])
                print Object, "Vars:", self.NHBDict[Sample][Object].shape
    def MergeVars(self, reprocess=False, others=""):
        debug = False
        # debug = True
        list_processes = []
        for Sample in self.SamplesNames:
            SubSample_ = self.Samples_dict[Sample][0]
            savename = self.NHBDict[Sample][SubSample_].outdir+"TopJet"+others+self.NHBDict[Sample][SubSample_].FileExt.replace("index","*")
            print savename
            savename = savename.replace("TopJet","/temp/TopJet").replace(SubSample_,Sample)
            print savename
            os.system("mkdir -p "+savename[:savename.rfind("/")])
            if debug: print Sample,others
            if debug: print "Save as", savename[savename.find("temp"):]
            index =0
            nElements = 0
            TJs = []
            PFs = []
            for SubSample in self.Samples_dict[Sample] if not reprocess else [SubSample_]:
                filenameTemplate = self.NHBDict[Sample][SubSample].outdir+"TopJet"+others+self.NHBDict[Sample][SubSample].FileExt.replace("index","*")
                if debug: print filenameTemplate
                if reprocess: filenameTemplate = filenameTemplate.replace(SubSample,Sample)
                if debug: print "Loop in", filenameTemplate[filenameTemplate.find("TopJet"):],
                filenameTemplate = sorted(glob(filenameTemplate))
                if debug: print len(filenameTemplate)
                findex = 0
                for fn in filenameTemplate:
                    findex +=1
                    print "Opening file ", findex, "out of", len(filenameTemplate),
                    print "\t", fn[fn.find("TopJet"):]
                    TopJet = np.load(fn)
                    if len(TopJet)>0:
                        PF = np.load(fn.replace("TopJet","PFParticles"))
                        if debug: print "\t pre ",  TopJet.shape, PF.shape
                        nElements += len(TopJet)
                        TJs.append(TopJet)
                        PFs.append(PF)
                        isFirst = True
                        while nElements >= self.filesize:
                            TJs = np.concatenate(TJs)
                            PFs = np.concatenate(PFs)
                            if isFirst:
                                TJs, PFs = shuffle(TJs, PFs, random_state=0)
                                isFirst = False
                            if not debug: np.save(savename.replace("*",str(index)),TJs[:self.filesize].astype(self.FS))
                            if not debug: np.save(savename.replace("*",str(index)).replace("TopJet","PFParticles"),PFs[:self.filesize].astype(self.FS))
                            TJs = [TJs[self.filesize:]]
                            PFs = [PFs[self.filesize:]]
                            nElements = len(TJs[0])
                            if debug: print "\tSaving index:", index
                            index += 1
                    if not debug:
                        if os.path.exists(fn): os.remove(fn)
                        if os.path.exists(fn.replace("TopJet","PFParticles")): os.remove(fn.replace("TopJet","PFParticles"))
            if nElements!=0:
                TJs = np.concatenate(TJs)
                PFs = np.concatenate(PFs)
                if not debug: np.save(savename.replace("*",str(index)),TJs.astype(self.FS))
                if not debug: np.save(savename.replace("*",str(index)).replace("TopJet","PFParticles"),PFs.astype(self.FS))
                if debug: print "Saving index:", index
            if debug: print "MERGING", others
            for Object in ["TopJet","PFParticles"]:
                for el in glob(savename.replace("TopJet",Object)):
                    if not debug: os.system("mv "+el+" "+el.replace("/temp/","/"))
    def SetRanges(self,var):
        min, max, bin = (0,100,100)
        if "tags" in var:                                           min, max, bin = (-0.5,5.5,6)
        if "charge" in var:                                         min, max, bin = (-2,2,5)
        if "Area" in var:                                           min, max, bin = (0,10,10)
        if "Fraction" in var or "tau" in var or "weight" in var:    min, max, bin = (0,1,100)
        if "btag" in var:                                           min, max, bin = (-1,1,100)
        if "eta" in var or "phi" in var:                            min, max, bin = (-math.pi,math.pi,100)
        if "pt" in var or "energy" in var:                          min, max, bin = (0,1000,100)
        if "PFpt" in var or "PFenergy" in var:                      min, max, bin = (0,1000,100)
        if "mass" in var:                                           min, max, bin = (0,200,100)
        if "particleID" in var:                                     min, max, bin = (-0.5,9.5,10)
        return min, max, bin
    @timeit
    def PlotVars(self,extraText="",Samples=[], Mode="Vars"):
        ROOT.gROOT.SetBatch(ROOT.kFALSE)
        ROOT.gROOT.SetBatch(ROOT.kTRUE)
        ROOT.gStyle.SetOptStat(1)
        ROOT.gStyle.SetOptStat(0)
        if len(Samples)==0:
            Samples = self.SamplesNames.keys()
        canvases = []
        Objects = self.RootObjects if Mode=="Vars" else self.matchingNames if Mode=="ReshapedVars" else None
        for Object in Objects:
            Vars = self.RootVarNameListsDict[Object] if Mode=="Vars" else self.VarNames
            for indexVarName, VarName in enumerate(Vars):
                min, max, bin = self.SetRanges(VarName)
                c_ = tdrCanvas(Object+VarName, min, max, 1.e-4, 1.e02, VarName, "A.U.")
                c_.SetLogy(1)
                leg = tdrLeg(0.50,0.60,0.9,0.9, 0.03)
                histos = []
                for Sample in Samples:
                    h_ = ROOT.TH1F(Object+VarName+Sample, "; "+VarName+"; A.U.", bin, min, max)
                    histos.append(h_)
                    if Mode=="Vars":
                        if Object=="TopJet":
                            array_ = np.concatenate([self.NHBDict[Sample].Vars[Object][:,:,i] for i in range(self.NHBDict[Sample].Vars[Object].shape[2])])
                        if Object=="Event":
                            array_ = self.NHBDict[Sample].Vars[Object]
                    else:
                        array_ = np.sum(self.NHBDict[Sample].ReshapeVars[Object],axis=1)/self.n_PF
                    fill_hist(h_, array_[:,indexVarName])
                    if (h_.Integral()!=0):
                        h_.Scale(1./h_.Integral())
                    tdrDraw(h_, "same", ROOT.kFullCircle, self.colors[Sample], 1, self.colors[Sample], 0, self.colors[Sample])
                    leg.AddEntry(h_, Sample ,"lep")
                leg.SetLineColorAlpha(1, 0.7)
                leg.Draw("same")
                c_name = self.pathPlots+Object+VarName+extraText+".pdf" if Mode=="Vars" else self.pathPlots+VarName+extraText+"_"+Object+".pdf"
                c_.SaveAs(c_name, "pdf");
                canvases.append(c_)
                time.sleep(2)
        return canvases
