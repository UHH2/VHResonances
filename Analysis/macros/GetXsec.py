#!/usr/bin/env python
from ModuleRunner import *
# to include CrossSectionHelper:
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-datasets'))
from CrossSectionHelper import MCSampleValuesHelper
import argparse

class GetXsec(VariablesBase):
    ''' Class to Get xsec info '''
    def __init__(self):
        VariablesBase.__init__(self)
        self.helper = MCSampleValuesHelper()
        # print 'Samples_Year_Dict', self.Samples_Year_Dict, '\n'
        # print 'SubSamples_Year_Dict', self.SubSamples_Year_Dict, '\n'
        # print 'Processes_Year_Dict', self.Processes_Year_Dict, '\n'
        # print 'AllSubSamples_List', self.AllSubSamples_List, '\n'
        # print 'AllProcesses_List', self.AllProcesses_List, '\n'

    def PrintValues(self, info_ ='lumi'):
        samples = []
        for mode in ['', '_inv']:
            for sample in self.SignalSamples:
                if ('_inv' in mode)!=('_inv' in sample): continue
                samples.append(sample)
        samples.extend(list(filter(lambda x: not 'inv' in x, self.Generic_SubSamples_Dict_['MC_DY'])))
        for s in ['MC_TTbar', 'MC_WW', 'MC_WZ', 'MC_ZZ', 'MC_WJets']:
            samples.extend(self.Generic_SubSamples_Dict_[s])
        samples.extend(list(filter(lambda x: 'inv' in x, self.Generic_SubSamples_Dict_['MC_DY'])))

        for year in self.years:
            for sample in samples:
                kFactor = not ('DY' in sample or 'WJet' in sample or 'ZJet' in sample)
                sample_ = sample
                sample_ = sample_.replace('MC_','')
                sample_ = sample_.replace('ZprimeToZH_M', 'ZprimeToZHToZlepHinc-')
                sample_ = sample_.replace('ZprimeToZH_inv_M', 'ZprimeToZHToZinvHinc-')
                sample_ = sample_.replace('DY_HT', 'DYJetsToLL_M-50_HT-')
                sample_ = sample_.replace('DY_inv_HT', 'ZJetsToNuNu_HT-')
                sample_ = sample_.replace('WJetsToLNu_HT', 'WJetsToLNu_HT-')
                sample += '_'+year
                spaces = ' '*(40-len(sample))
                if info_ =='lumi':
                    lumi = str(self.helper.get_lumi(sample_,'13TeV',year, kFactor=kFactor))
                    spaces_lumi = ' '*(15-len(lumi))
                    print('<InputData Type="MC"    Version="'+sample+'"'+spaces+'Lumi="'+lumi+'"'+spaces_lumi+'NEventsMax="&NEVT;" Cacheable="&CACHEABLE;">  &'+sample+';'+spaces+'<InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>')
                if info_ =='xml':
                    xml = str(self.helper.get_xml(sample_,'13TeV',year))
                    print('<!ENTITY '+sample+spaces+'SYSTEM  "CMSSW_BASE/src/UHH2/common/UHH2-datasets/'+xml+'">')

        for year in self.years:
            for mode in ['DATA_SingleElectron', 'DATA_SingleMuon', 'DATA_MET']:
                for sample in sorted(self.Generic_SubSamples_Dict_[mode]):
                    if not sample[-1] in self.RunPeriods_Dict[year]: continue
                    if year =='UL18' and 'SinglePhoton' in sample: continue
                    sample_ = sample
                    if year =='UL18':
                        sample_ = sample_.replace('SingleElectron', 'EGamma')
                    sample_ = sample_.replace('DATA_','')
                    sample += '_'+year
                    spaces = ' '*(40-len(sample))
                    if info_ =='lumi':
                        lumi = "1"
                        spaces_lumi = ' '*(15-len(lumi))
                        print('<InputData Type="DATA"  Version="'+sample+'"'+spaces+'Lumi="'+lumi+'"'+spaces_lumi+'NEventsMax="&NEVT;" Cacheable="&CACHEABLE;">  &'+sample+';'+spaces+'<InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>')
                    if info_ =='xml':
                        xml = str(self.helper.get_xml(sample_,'13TeV',year))
                        print('<!ENTITY '+sample+spaces+'SYSTEM  "CMSSW_BASE/src/UHH2/common/UHH2-datasets/'+xml+'">')





if __name__ == '__main__':

    XsecInfo = GetXsec()
    XsecInfo.PrintValues('xml')
    XsecInfo.PrintValues('lumi')
