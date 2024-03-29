
from Utils import *
from fileManipulation import *
from xml.dom.minidom import parseString


YearVars = {}
YearVars['JEC_Version']             = {'UL16preVFP': 'Summer19UL16APV_V7',
                                       'UL16postVFP': 'Summer19UL16_V7',
                                       'UL17': 'Summer19UL17_V5',
                                       'UL18': 'Summer19UL18_V5',
                                       }
YearVars['JER_Version']             = {'UL16preVFP': 'Summer20UL16APV_JRV3',
                                       'UL16postVFP': 'Summer20UL16_JRV3',
                                       'UL17': 'Summer19UL17_JRV2',
                                       'UL18': 'Summer19UL18_JRV2',
                                       }
YearVars['lumi_file']               = {'UL16preVFP': os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/UL16preVFP/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16preVFP_normtag.root',
                                       'UL16postVFP': os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/UL16postVFP/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16postVFP_normtag.root',
                                       'UL17': os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/UL17/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON_normtag.root',
                                       'UL18': os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/UL18/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_normtag.root',
                                       }

YearVars['MCBtagEfficiencies']      = {'UL16preVFP': os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/BTag/SF_UL16preVFP.root',
                                       'UL16postVFP': os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/BTag/SF_UL16postVFP.root',
                                       'UL17': os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/BTag/SF_UL17.root',
                                       'UL18': os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/BTag/SF_UL18.root',
                                       }

YearVars['BTagCalibration_FixedWP'] = {'UL16preVFP': os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/BTag/subjet_deepCSV_106XUL16preVFP_v1.csv',
                                       'UL16postVFP': os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/BTag/subjet_deepCSV_106XUL16postVFP_v1.csv',
                                       'UL17': os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/BTag/subjet_DeepCSV_106X_UL17_SF.csv',
                                       'UL18': os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/BTag/subjet_deepCSV_106XUL18_v1.csv',
                                       }


def newNumber(year,sample,ConfigFile,syst,channel):
    newNumber = 20
    if 'DATA' in sample:
        if 'UL16' in year:
            newNumber = 50
            if any(x in sample for x in ['DATA_SingleMuon_RunF', 'DATA_SingleElectron_RunC', 'DATA_SingleElectron_RunD', 'DATA_MET_RunB', 'DATA_MET_RunC', 'DATA_MET_RunF']):
                newNumber = 250
            if any(x in sample for x in ['DATA_SingleElectron_RunE']):
                newNumber = 200
            if any(x in sample for x in ['DATA_SingleElectron_RunB']):
                newNumber = 180
            if any(x in sample for x in ['DATA_SinglePhoton_RunF']):
                newNumber = 150
            if any(x in sample for x in ['DATA_MET_RunD', 'DATA_MET_RunE', 'DATA_MET_RunG']):
                newNumber = 170
        if year=='UL17':
            newNumber = 100
            if any(x in sample for x in ['DATA_MET_RunD', 'DATA_SinglePhoton_RunE','DATA_SinglePhoton_RunF','DATA_MET_RunF','DATA_MET_RunE']):
                newNumber = 50
            if any(x in sample for x in ['DATA_MET_RunB', 'DATA_SinglePhoton_RunC','DATA_SingleElectron_RunB','DATA_SingleElectron_RunD']):
               newNumber = 150
        if year=='UL18':
            newNumber = 50
            if any(x in sample for x in ['DATA_SingleMuon_RunC_UL18', 'DATA_SingleMuon_RunB_UL18', 'DATA_MET_RunB_UL18']):
                newNumber = 100
            if any(x in sample for x in ['DATA_SingleElectron_RunD_UL18','DATA_SingleElectron_RunA_UL18','DATA_MET_RunD_UL18']):
                newNumber = 20
    if 'MC_DY' in sample:
        newNumber = 15 if 'UL18' in sample else 60
        if any(x in sample for x in ['MC_DY_HT400to600_UL17', 'MC_DY_HT600to800_UL17', 'MC_DY_inv_HT2500toInf', 'MC_DY_HT2500toInf_UL16postVFP', 'MC_DY_HT1200to2500_UL16postVFP']):
            newNumber = 30
        if any(x in sample for x in ['MC_DY_HT800to1200_UL17', 'MC_DY_HT1200to2500_UL17', 'MC_DY_HT2500toInf_UL17', 'MC_DY_inv_HT100to200_UL18']):
            newNumber = 70
        if any(x in sample for x in ['MC_DY_HT400to600_UL16postVFP', 'MC_DY_HT200to400_UL17', 'MC_DY_inv_HT600to800_UL17','MC_DY_inv_HT400to600_UL17', 'MC_DY_HT800to1200_UL16postVFP', 'MC_DY_HT100to200_UL16postVFP', 'MC_DY_HT100to200_UL18']):
            newNumber = 100
        if any(x in sample for x in ['MC_DY_inv_HT200to400_UL17', 'MC_DY_inv_HT800to1200_UL16postVFP']):
            newNumber = 140
        if any(x in sample for x in ['MC_DY_HT200to400_UL16postVFP', 'MC_DY_inv_HT1200to2500_UL16postVFP']):
            newNumber = 180
        if any(x in sample for x in ['MC_DY_HT100to200_UL17', 'MC_DY_inv_HT100to200_UL16postVFP']):
            newNumber = 240
    if 'MC_TT' in sample:
        newNumber = 15 if 'UL16' in year else 60
    if 'MC_W' in sample:
        newNumber = 200
        if any(x in sample for x in ['MC_WZTo2L2Q_UL16postVFP', 'MC_WZ_UL16postVFP']):
            newNumber = 40
        if any(x in sample for x in ['MC_WJets', 'MC_WW_UL16postVFP', 'MC_WJetsToLNu_HT1200to2500', 'MC_WJetsToLNu_HT400to600']):
            newNumber = 13 if 'inv' in channel else 30
        if any(x in sample for x in ['MC_WJetsToLNu_HT200to400', 'MC_WJetsToLNu_HT600to800', 'MC_WJetsToLNu_HT100to200', 'MC_WJetsToLNu_HT2500toInf_UL16postVFP', 'MC_WJetsToLNu_HT800to1200']):
            newNumber = 20
        if any(x in sample for x in ['MC_WJetsToLNu_HT100to200_UL16postVFP']):
            newNumber = 300
    if 'MC_ZZ' in sample:
        newNumber = 240
        if any(x in sample for x in ['MC_ZZTo4Q_UL16postVFP', 'MC_ZZ_UL16postVFP']):
            newNumber = 100 if 'inv' in channel else 150
        if any(x in sample for x in ['MC_ZZ_UL17']):
            newNumber = 500 if 'inv' in channel else 850
        if any(x in sample for x in ['MC_ZZTo2L2Nu_UL16postVFP']):
            newNumber = 100 if 'inv' in channel else 250
    if 'MC_ZJets' in sample:
        newNumber = 40
    if 'MC_QCD' in sample:
        newNumber = 40
        if any(x in sample for x in ['MC_QCD_HT500to700', 'MC_QCD_HT700to1000']):
            newNumber = 50
        if any(x in sample for x in ['MC_QCD_HT300to500']):
            newNumber = 130
        if any(x in sample for x in ['MC_QCD_HT200to300']):
            newNumber = 130
        if any(x in sample for x in ['MC_QCD_HT100to200_UL17']):
            newNumber = 185
        if any(x in sample for x in ['MC_QCD_HT100to200_UL18']):
            newNumber = 300
    if 'MC_ZprimeToZH' in sample:
        newNumber = 100
    if not 'Preselection' in ConfigFile and not 'SF' in ConfigFile and not 'LeptonIDStudies' in ConfigFile:
        newNumber = 1 if 'MC_DY' in sample else (10 if 'MC_TT' in sample or 'DATA' in sample else 1000)
        if 'MC_WJets' in sample: newNumber = 10
    if syst!='nominal':
        newNumber = int(0.9*newNumber)
    if 'Preselection' in ConfigFile:
        newNumber = int(0.9*newNumber)
    if 'LeptonIDStudies' in ConfigFile:
        if not 'MC_ZprimeToZH' in sample: newNumber = int(newNumber/3)
    if 'HEMIssueStudy' in ConfigFile:
        newNumber = int(100)
    # if 'UL17' in ConfigFile:
    #     isFast = False
    #     isFast = True
    #     isToReduce = any(x in sample for x in ['MC_ZJetsToQQ', 'MC_WJetsToQQ', 'MC_QCD'])
    #     newNumber = 10 if isToReduce else 50
    #     if isFast: newNumber = 5 if isToReduce else 10
    #     changes.append(NFileLine(newNumber))
    return str(max(1,int(newNumber/(defaulTimePerJob/TimePerJob))))
    # return str(max(1,int(newNumber/(defaulTimePerJob/1))))


@timeit
def CreateConfigFiles(year, samples, all_samples, collections, channels, systematics, controls, original_dir, SubmitDir, ConfigFile, Path_SFRAME, lumi):
    with open(original_dir+'config/'+ConfigFile, 'r') as search:
        for line in search:
            if '<ConfigSGE' in line:
                ConfigSGE = parseString(line).getElementsByTagName('ConfigSGE')[0]
                global TimePerJob; global defaulTimePerJob; defaulTimePerJob = 3.
                TimePerJob = int(ConfigSGE.attributes['TIME'].value) if ConfigSGE.hasAttribute('TIME') else defaulTimePerJob
    outdir = ConfigFile[0:ConfigFile.find('Config')]

    for collection in collections:
        for channel in channels:
            for syst in systematics:
                if ('Muon' in syst and not 'muon' in channel): continue
                folders = collection+'/'+channel+'channel/'+syst+'/'
                a = os.system('mkdir -p '+Path_SFRAME+outdir+'/'+folders)
                path = SubmitDir+folders
                if not os.path.exists(path):
                    os.makedirs(path)
                for sample in samples:
                    if DoControl(controls, collection+channel+syst+sample, channel, sample):
                        continue
                    filename = outdir+'_'+sample+'.xml'
                    a = os.system('cp '+original_dir+'config/'+ConfigFile+' '+path+filename)
                    a = os.system('cp '+original_dir+'JobConfig.dtd '+path)
                    comments = []
                    for el in all_samples:
                        if sample == el: continue
                        if 'MC' in el:
                            comments.append(['<InputData', 'Type', 'MC',   '"'+el+'"'])
                        elif 'DATA' in el:
                            comments.append(['<InputData', 'Type', 'DATA', '"'+el+'"'])
                    if 'Puppi' in collection:
                        comments.append(['<Item Name', 'JetCollection',           'Value', 'jetsAk4CHS'])
                        comments.append(['<Item Name', 'TopJetCollection',        'Value', 'jetsAk8CHSSubstructure_SoftDropCHS'])
                        comments.append(['<Item Name', 'additionalBranches',      'Value', 'hotvrGen hotvrPuppi jetsAk4Puppi'])
                    if 'HOTVR' in collection:
                        comments.append(['<Item Name', 'JetCollection',           'Value', 'jetsAk4CHS'])
                        comments.append(['<Item Name', 'TopJetCollection',        'Value', 'jetsAk8CHSSubstructure_SoftDropCHS'])
                        # comments.append(['<Item Name', 'GenTopJetCollection',     'Value', 'genjetsAk8SubstructureSoftDrop'])
                        comments.append(['<Item Name', 'TopPuppiJetCollection',   'Value', 'jetsAk8PuppiSubstructure_SoftDropPuppi'])
                        comments.append(['<Item Name', 'additionalBranchesPuppi', 'Value', 'jetsAk4Puppi'])
                    if 'CHS' in collection:
                        comments.append(['<Item Name', 'additionalBranchesPuppi', 'Value', 'jetsAk4Puppi'])
                        comments.append(['<Item Name', 'additionalBranches',      'Value', 'hotvrGen hotvrPuppi jetsAk4Puppi'])
                        comments.append(['<Item Name', 'TopPuppiJetCollection',   'Value', 'jetsAk8PuppiSubstructure_SoftDropPuppi'])
                    comment_lines(path, filename, comments, remove=True)
                    changes = []
                    # Change anmalara when creating xml. change also email
                    changes.append(['Mail=', 'USER@mail.desy.de', 'USER@mail.desy.de', os.environ['USER']+'@mail.desy.de'])
                    changes.append(['<!ENTITY', '/nfs/dust/cms/user/USER', 'USER', os.environ['USER']])
                    changes.append(['<!ENTITY', 'CMSSW_BASE', 'CMSSW_BASE', os.environ['CMSSW_BASE']])
                    change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
                    changes.append(['<ConfigSGE', 'Workdir', 'workdir_'+outdir, 'workdir_'+outdir+'_'+sample])
                    changes.append(['<ConfigParse', 'FileSplit="20"', 'FileSplit="20"', 'FileSplit="'+newNumber(year,sample,ConfigFile,syst,channel)+'"'])
                    changes.append(['<!ENTITY', 'OUTDIR', outdir, outdir+'/'+folders])
                    if 'Preselection' in ConfigFile or 'PDFReweight' in ConfigFile or 'HEMIssueStudy' in ConfigFile or 'LeptonIDStudies' in ConfigFile:
                        if 'invisible' in channel:
                            changes.append(['<!ENTITY', 'original_pdfname', 'defaultValue', 'NNPDF31_nnlo_as_0118_nf_4'])
                        else: changes.append(['<!ENTITY', 'original_pdfname', 'defaultValue', 'NNPDF31_lo_as_0130'])
                    if 'Selection' in ConfigFile:
                        changes.append(['<!ENTITY', 'SYSTEM', 'Preselection/All/leptonchannel/nominal/', 'Preselection/'+folders.replace('MuonScale_up','nominal').replace('MuonScale_down','nominal')])
                        changes.append(['<!ENTITY', 'SYSTEM', 'Preselection/All/invisiblechannel/nominal/', 'Preselection/'+folders.replace('MuonScale_up','nominal').replace('MuonScale_down','nominal')])
                    if 'SignalRegion' in ConfigFile:
                        changes.append(['<!ENTITY', 'SYSTEM', 'Selection/All/leptonchannel/nominal/', 'Selection/'+folders])
                        changes.append(['<!ENTITY', 'SYSTEM', 'Selection/All/invisiblechannel/nominal/', 'Selection/'+folders])
                    if 'Puppi' in collection:
                        changes.append(['<!ENTITY', 'isCHS',    '"true"', '"false"'])
                        changes.append(['<!ENTITY', 'isHOTVR',  '"true"', '"false"'])
                        changes.append(['<Item Name', 'METName', 'slimmedMETs', 'slimmedMETsPuppi'])
                        changes.append(['<Item Name', 'additionalBranchesPuppi', '"additionalBranchesPuppi"', '"additionalBranches"     '])
                    if 'CHS' in collection:
                        changes.append(['<!ENTITY', 'isPuppi',  '"true"', '"false"'])
                        changes.append(['<!ENTITY', 'isHOTVR',  '"true"', '"false"'])
                    if 'HOTVR' in collection:
                        changes.append(['<!ENTITY', 'isPuppi',  '"true"', '"false"'])
                        changes.append(['<!ENTITY', 'isCHS',    '"true"', '"false"'])
                        changes.append(['<!ENTITY', 'isCHS',    '"true"', '"false"'])
                        changes.append(['<Item Name', 'METName', 'slimmedMETs', 'slimmedMETsPuppi'])
                    if 'muon' in channel:
                        changes.append(['<!ENTITY', 'muonchannel',   '"false"', '"true"'])
                    if 'electron' in channel:
                        changes.append(['<!ENTITY', 'electronchannel',   '"false"', '"true"'])
                    if 'invisible' in channel:
                        changes.append(['<!ENTITY', 'invisiblechannel',   '"false"', '"true"'])
                    if 'charm' in channel:
                        changes.append(['<!ENTITY', 'charmchannel',   '"false"', '"true"'])
                    changes.append(['<!ENTITY', 'YEAR', 'defaultValue', year])
                    changes.append(['<Cycle', 'TargetLumi', 'defaultValue', str(lumi)])
                    for var in YearVars:
                        changes.append(['<!ENTITY', var, 'defaultValue', YearVars[var][year]])
                    for syst_ in ['JER','JEC']:
                        if syst_.lower() in syst.lower():
                            changes.append(['<!ENTITY', syst_.lower()+'smear_direction', '"nominal"', '"'+syst.lower().replace(syst_.lower()+'_','')+'"' ])
                    if 'MuonScale' in syst:
                        changes.append(['<!ENTITY', 'MuonScaleVariations', '"nominal"', '"'+syst.replace('MuonScale_','')+'"' ])
                    change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
