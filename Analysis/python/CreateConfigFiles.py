
from Utils import *
from fileManipulation import *
from xml.dom.minidom import parseString


YearVars = {}
YearVars["JEC_Version"]         = {"2016": "Summer16_07Aug2017_V11",
                                   "2017": "Fall17_17Nov2017_V32",
                                   "2018": "Autumn18_V19",
                                   }
YearVars["JER_Version"]         = {"2016": "Summer16_25nsV1",
                                   "2017": "Fall17_V3",
                                   "2018": "Autumn18_V7b",
                                   }
YearVars["lumi_file"]           = {"2016": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2016/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root",
                                   "2017": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.root",
                                   "2018": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2018/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.root",
                                   }

YearVars["MCBtagEfficiencies"]  = {"2016": os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/OtherPlots/BTag/SF_2016.root",
                                   "2017": os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/OtherPlots/BTag/SF_2017.root",
                                   "2018": os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/OtherPlots/BTag/SF_2018.root",
                                   }

YearVars["BTagCalibration"]     = {"2016": os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/OtherPlots/BTag/subjet_DeepCSV_102XSF_V1.csv",
                                   "2017": os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/OtherPlots/BTag/subjet_DeepCSV_2016LegacySF_V1.csv",
                                   "2018": os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/OtherPlots/BTag/subjet_DeepCSV_94XSF_V4_B_F.csv",
                                   }

# YearVars["BTagCalibration"]     = {"2016": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2016/DeepCSV_2016LegacySF_WP_V1.csv",
#                                    "2017": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2017/DeepCSV_94XSF_WP_V4_B_F.csv",
#                                    "2018": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2018/DeepCSV_102XSF_WP_V1.csv",
#                                    }



def newNumber(year,sample,ConfigFile,syst):
    newNumber = 20
    if "DATA" in sample:
        if year=="2016":
            newNumber = 350
            if any(x in sample for x in ["DATA_SingleMuon_RunF", "DATA_SingleMuon_RunG", "DATA_SingleMuon_RunH"]):
                newNumber = 300
            if any(x in sample for x in ["DATA_SingleElectron_RunC", "DATA_SingleElectron_RunD", "DATA_SingleElectron_RunF"]):
                newNumber = 250
            if any(x in sample for x in ["DATA_SingleElectron_RunE", "DATA_SingleElectron_RunG", "DATA_SingleElectron_RunH"]):
                newNumber = 200
            if any(x in sample for x in ["DATA_SingleElectron_RunB"]):
                newNumber = 180
        if year=="2017":
            newNumber = 340
            if any(x in sample for x in ["DATA_SingleMuon_RunF"]):
                newNumber = 280
            if any(x in sample for x in ["DATA_SingleElectron_RunF", "DATA_SingleMuon_RunE"]):
                newNumber = 300
            if any(x in sample for x in ["DATA_SingleElectron_RunD"]):
                newNumber = 390
            if any(x in sample for x in ["DATA_SingleMuon_RunB"]):
                newNumber = 430
            if any(x in sample for x in ["DATA_SingleElectron_RunB"]):
                newNumber = 520
        if year=="2018":
            newNumber = 80
            if any(x in sample for x in ["DATA_SingleMuon_RunA", "DATA_SingleMuon_RunB", "DATA_SingleElectron_RunC", "DATA_SingleMuon_RunC", "DATA_SingleMuon_RunD"]):
                newNumber = 150
    if "MC_DY" in sample:
        newNumber = 60
        if any(x in sample for x in ["MC_DY_HT600to800_2016", "MC_DY_HT800to1200_2016", "MC_DY_HT800to1200_2017", "MC_DY_HT600to800_2017"]):
            newNumber = 70
        if any(x in sample for x in ["MC_DY_HT400to600_2017"]):
            newNumber = 80
        if any(x in sample for x in ["MC_DY_HT400to600_2016", "MC_DY_HT400to600_2018"]):
            newNumber = 100
        if any(x in sample for x in ["MC_DY_HT200to400_2017"]):
            newNumber = 140
        if any(x in sample for x in ["MC_DY_HT200to400_2018"]):
            newNumber = 155
        if any(x in sample for x in ["MC_DY_HT200to400_2016"]):
            newNumber = 200
        if any(x in sample for x in ["MC_DY_HT100to200_2017"]):
            newNumber = 240
        if any(x in sample for x in ["MC_DY_HT100to200_2018"]):
            newNumber = 330
        if any(x in sample for x in ["MC_DY_HT100to200_2016"]):
            newNumber = 300
        if any(x in sample for x in ["MC_DY_inv_PtZ_50To100", "MC_DY_inv_PtZ_100To250", "MC_DY_inv_PtZ_250To400", "MC_DY_inv_PtZ_400To650", "MC_DY_inv_PtZ_650ToInf", "MC_DY_inv_HT_100To200", "MC_DY_inv_HT_200To400", "MC_DY_inv_HT_400To600", "MC_DY_inv_HT_600To800", "MC_DY_inv_HT_800To1200", "MC_DY_inv_HT_1200To2500", "MC_DY_inv_HT_2500ToInf"]):
            newNumber = 100
    if "MC_TT" in sample:
        newNumber = 40 if year=="2016" else 200 if year=="2017" else 60
    if "MC_W" in sample:
        newNumber = 230
        if any(x in sample for x in ["MC_WWTo2L2Nu_2016", "MC_WZTo1L3Nu_2016", "MC_WZTo3LNu_2016"]):
            newNumber = 300
        if any(x in sample for x in ["MC_WZ_2017"]):
            newNumber = 700
        if any(x in sample for x in ["MC_WZ_2018"]):
            newNumber = 800
        if any(x in sample for x in ["MC_WW_2018"]):
            newNumber = 900
        if any(x in sample for x in ["MC_WZTo2L2Q_2016"]):
            newNumber = 40
    if "MC_ZZ" in sample:
        newNumber = 250
        if any(x in sample for x in ["MC_ZZTo4Q_2016"]):
            newNumber = 150
        if any(x in sample for x in ["MC_ZZ_2017"]):
            newNumber = 850
        if any(x in sample for x in ["MC_ZZ_2018"]):
            newNumber = 1000
    if "MC_ZprimeToZH" in sample:
        newNumber = 100
    if not "Preselection" in ConfigFile and not "SF" in ConfigFile and not "LeptonIDStudies" in ConfigFile:
        newNumber = 1 if "MC_DY" in sample else 1000
    if syst!="nominal":
        newNumber = int(0.9*newNumber)
    if "LeptonIDStudies" in ConfigFile:
        if not "MC_ZprimeToZH" in sample: newNumber = int(newNumber/3)
    # if "2017" in ConfigFile:
    #     isFast = False
    #     isFast = True
    #     isToReduce = any(x in sample for x in ["MC_ZJetsToQQ", "MC_WJetsToQQ", "MC_QCD"])
    #     newNumber = 10 if isToReduce else 50
    #     if isFast: newNumber = 5 if isToReduce else 10
    #     changes.append(NFileLine(newNumber))
    return str(max(1,int(newNumber/(defaulTimePerJob/TimePerJob))))

@timeit
def CreateConfigFiles(year, samples, all_samples, collections, channels, systematics, controls, original_dir, SubmitDir, ConfigFile, Path_SFRAME, lumi):
    with open(original_dir+"config/"+ConfigFile, "r") as search:
        for line in search:
            if "<ConfigSGE" in line:
                ConfigSGE = parseString(line).getElementsByTagName('ConfigSGE')[0]
                global TimePerJob; global defaulTimePerJob; defaulTimePerJob = 3.
                TimePerJob = int(ConfigSGE.attributes['TIME'].value) if ConfigSGE.hasAttribute('TIME') else defaulTimePerJob
    outdir = ConfigFile[0:ConfigFile.find("Config")]
    for collection in collections:
        for channel in channels:
            for syst in systematics:
                folders = collection+"/"+channel+"channel/"+syst+"/"
                a = os.system("mkdir -p "+Path_SFRAME+outdir+"/"+folders)
                path = SubmitDir+folders
                if not os.path.exists(path):
                    os.makedirs(path)
                for sample in samples:
                    if ("Electron" in sample and not "electron" in channel) : continue
                    if ("Muon" in sample and not "muon" in channel) : continue
                    if ("MET" in sample and not "invisible" in channel) : continue

                    if all(not control in collection+channel+syst+sample for control in controls):
                        continue
                    filename = outdir+"_"+sample+".xml"
                    a = os.system("cp "+original_dir+"config/"+ConfigFile+" "+path+filename)
                    a = os.system("cp "+original_dir+"JobConfig.dtd "+path)
                    comments = []
                    for el in all_samples:
                        if sample == el: continue
                        if "MC" in el:
                            comments.append(["<InputData", "Type", "MC",   '"'+el+'"'])
                        elif "DATA" in el:
                            comments.append(["<InputData", "Type", "DATA", '"'+el+'"'])
                    if "Puppi" in collection:
                        comments.append(["<Item Name", "JetCollection",           'Value', 'jetsAk4CHS'])
                        comments.append(["<Item Name", "TopJetCollection",        'Value', 'jetsAk8CHSSubstructure_SoftDropCHS'])
                        comments.append(["<Item Name", "additionalBranches",      'Value', 'hotvrGen hotvrPuppi jetsAk4Puppi'])
                    if "HOTVR" in collection:
                        comments.append(["<Item Name", "JetCollection",           'Value', 'jetsAk4CHS'])
                        comments.append(["<Item Name", "TopJetCollection",        'Value', 'jetsAk8CHSSubstructure_SoftDropCHS'])
                        # comments.append(["<Item Name", "GenTopJetCollection",     'Value', 'genjetsAk8SubstructureSoftDrop'])
                        comments.append(["<Item Name", "TopPuppiJetCollection",   'Value', 'jetsAk8PuppiSubstructure_SoftDropPuppi'])
                        comments.append(["<Item Name", "additionalBranchesPuppi", 'Value', 'jetsAk4Puppi'])
                    if "CHS" in collection:
                        comments.append(["<Item Name", "additionalBranchesPuppi", 'Value', 'jetsAk4Puppi'])
                        comments.append(["<Item Name", "additionalBranches",      'Value', 'hotvrGen hotvrPuppi jetsAk4Puppi'])
                        comments.append(["<Item Name", "TopPuppiJetCollection",   'Value', 'jetsAk8PuppiSubstructure_SoftDropPuppi'])
                    comment_lines(path, filename, comments, remove=True)
                    changes = []
                    # TODO change anmalara when creating xml. change also email
                    changes.append(["Mail=", "USER@mail.desy.de", "USER@mail.desy.de", os.environ["USER"]+"@mail.desy.de"])
                    changes.append(["<!ENTITY", "/nfs/dust/cms/user/USER", "USER", os.environ["USER"]])
                    changes.append(["<!ENTITY", "CMSSW_BASE", "CMSSW_BASE", os.environ["CMSSW_BASE"]])
                    change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
                    changes.append(["<ConfigSGE", "Workdir", "workdir_"+outdir, "workdir_"+outdir+"_"+sample])
                    changes.append(["<ConfigParse", 'FileSplit="20"', 'FileSplit="20"', 'FileSplit="'+newNumber(year,sample,ConfigFile,syst)+'"'])
                    changes.append(["<!ENTITY", "OUTDIR", outdir , outdir+"/"+folders])
                    if "Selection" in ConfigFile: # TODO How does this need to change with the invisiblechannel?
                        changes.append(["<!ENTITY", "SYSTEM", "Preselection/All/leptonchannel/nominal/" , "Preselection/"+folders])
                    if "SignalRegion" in ConfigFile:
                        changes.append(["<!ENTITY", "SYSTEM", "Selection/All/leptonchannel/nominal/" , "Selection/"+folders])
                    if "Puppi" in collection:
                        changes.append(["<!ENTITY", "isCHS",    '"true"', '"false"'])
                        changes.append(["<!ENTITY", "isHOTVR",  '"true"', '"false"'])
                        changes.append(["<Item Name", "METName", 'slimmedMETs', 'slimmedMETsPuppi'])
                        changes.append(["<Item Name", "additionalBranchesPuppi", '"additionalBranchesPuppi"', '"additionalBranches"     '])
                    if "CHS" in collection:
                        changes.append(["<!ENTITY", "isPuppi",  '"true"', '"false"'])
                        changes.append(["<!ENTITY", "isHOTVR",  '"true"', '"false"'])
                    if "HOTVR" in collection:
                        changes.append(["<!ENTITY", "isPuppi",  '"true"', '"false"'])
                        changes.append(["<!ENTITY", "isCHS",    '"true"', '"false"'])
                        changes.append(["<!ENTITY", "isCHS",    '"true"', '"false"'])
                        changes.append(["<Item Name", "METName", 'slimmedMETs', 'slimmedMETsPuppi'])
                    if "muon" in channel:
                        changes.append(["<!ENTITY", "muonchannel",   '"false"', '"true"'])
                    if "electron" in channel:
                        changes.append(["<!ENTITY", "electronchannel",   '"false"', '"true"'])
                    if "invisible" in channel:
                        changes.append(["<!ENTITY", "invisiblechannel",   '"false"', '"true"'])
                    changes.append(["<!ENTITY", "YEAR", 'defaultValue', year])
                    changes.append(["<Cycle", "TargetLumi", 'defaultValue', str(lumi)])
                    for var in YearVars:
                        changes.append(["<!ENTITY", var, 'defaultValue', YearVars[var][year]])
                    for syst_ in ["JER","JEC"]:
                        if syst_.lower() in syst.lower():
                            changes.append(["<!ENTITY", syst_.lower()+"smear_direction", '"nominal"', '"'+syst.lower().replace(syst_.lower()+"_","")+'"' ])
                    change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
