from Utils import *
sys.path.insert(0, os.getenv('CMSSW_BASE')+'/src/UHH2/common/UHH2-datasets')
from CrossSectionHelper import MCSampleValuesHelper

def CheckXsecInfo(ConfigFile='PreselectionConfig.xml'):
    helper = MCSampleValuesHelper()
    all_samples_config = VariablesBase().AllSubSamples_List
    all_samples_lumi = VariablesBase().AllSubSamples_List

    with open(ModuleRunnerBase().ConfigDir+ConfigFile, 'r') as config:
        lines = config.readlines()
    for line in lines:
        if not '<InputData' in line and not '<!ENTITY' in line: continue
        config, lumi = (None, None)
        if '<!ENTITY' in line and 'SYSTEM' in line:
            sample = line.split()[1]
            config = line.split()[3]
            sample, config = (sample.split('"')[0], config.split('"')[1])
        if '<InputData' in line:
            sample = line.split()[2]
            lumi = line.split()[3]
            sample, lumi = (sample.split('"')[1], float(lumi.split('"')[1]))
            if 'Run' in sample:
                all_samples_lumi.remove(sample)
                continue
        if config==lumi: continue
        if 'RunA' in sample and not '18' in sample:
            continue
        ref_sample = sample
        sample = sample.replace('MC_','').replace('DATA_','')
        sample, year = (sample[:sample.rfind('UL')-1], sample[sample.rfind('UL'):] )
        sample = sample.replace('DY_HT','DYJetsToLL_M-50_HT-').replace('DY_inv_HT','ZJetsToNuNu_HT-').replace('WJetsToLNu_HT','WJetsToLNu_HT-').replace('ZprimeToZH_M','ZprimeToZHToZlepHinc-').replace('ZprimeToZH_inv_M','ZprimeToZHToZinvHinc-')
        if '18' in year and 'SingleElectron' in sample:
            sample = sample.replace('SingleElectron','EGamma')
        if config!=None:
            config_ = helper.get_xml(sample,'13TeV',year)
            if not config_ in config:
                print sample, year, config, config_
            else:
                all_samples_config.remove(ref_sample)
        if lumi!=None:
            lumi_ = helper.get_lumi(sample,'13TeV',year, kFactor=True)
            diff_ = math.fabs(lumi-lumi_)
            if diff_>1e-03:
                print sample, year, lumi, lumi_, diff_
            else:
                all_samples_lumi.remove(ref_sample)
    print all_samples_config
    print all_samples_lumi

def main():
    CheckXsecInfo()

if __name__ == '__main__':
    main()
