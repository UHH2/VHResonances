from Utils import *
sys.path.insert(0, os.getenv('CMSSW_BASE')+'/src/UHH2/common/UHH2-datasets')
from CrossSectionHelper import MCSampleValuesHelper

helper = MCSampleValuesHelper()
energy = '13TeV'
ref_year = 'UL18'
HT = ['100to200','200to400','400to600','600to800','800to1200','1200to2500','2500toInf']
masses = ['600', '800', '1000', '1200', '1400', '1600', '1800', '2000', '2500', '3000', '3500', '4000', '4500', '5000', '5500', '6000', '7000', '8000']
DY = ["DYJetsToLL_M-50_HT-"+ht for ht in HT]
Znn = ["ZJetsToNuNu_HT-"+ht for ht in HT]
Wln = ["WJetsToLNu_HT-"+ht for ht in HT]
TT = ['TTToHadronic','TTToSemiLeptonic', 'TTTo2L2Nu']
VV = ['WW', 'WZ', 'ZZ']
Zprime = ["ZprimeToZHToZlepHinc-"+m for m in masses]
Zprime_inv = ["ZprimeToZHToZinvHinc-"+m for m in masses]

tables = [  (1, DY+ TT+ VV+ Zprime),
            (2, Znn+ Wln+ Zprime_inv)
]
for table, samples in tables:
    for sample in samples:
        name = str(sample)
        name = name.replace('JetsToLL_M-50', '')
        name = name.replace('ZJetsToNuNu', 'DY')
        name = name.replace('WJetsToLNu', 'WJets')
        name = name.replace('ZprimeToZHToZlepHinc-', 'ZprimeM-')
        name = name.replace('ZprimeToZHToZinvHinc-', 'ZprimeM-')
        name = name.replace('_','\_')
        gen = helper.get_xml(sample,energy,ref_year).split('/')[3].replace('_narrow_M','').replace('0To','0to').replace(sample,'').replace('CP5','').replace('PSweights','').replace('Summer20UL18','').replace('_v1','').replace('_v2','').replace('.xml','').replace('_','').replace('-','')
        gen = gen.replace(sample.replace('_HT-','HT'),'')
        if 'madgraph' in gen and not 'MLM' in gen: gen += 'MLM'
        if gen !='pythia8': gen = gen.replace('pythia8','')
        xsec= helper.get_xs(sample,energy,ref_year)
        xsec *= helper.get_br(sample, energy, ref_year)
        xsec *= helper.get_kfactor(sample, energy, ref_year)
        xsec = str(round(xsec,2))
        nev_16APV = "{:.3e}".format(helper.get_nevt(sample,energy,'UL16preVFP'))
        nev_16 = "{:.3e}".format(helper.get_nevt(sample,energy,'UL16postVFP'))
        nev_17 = "{:.3e}".format(helper.get_nevt(sample,energy,'UL17'))
        nev_18 = "{:.3e}".format(helper.get_nevt(sample,energy,'UL18'))
        print(name+' & '+gen+' & '+xsec+' & '+nev_16APV+' & '+nev_16+' & '+nev_17+' & '+nev_18+' \\\\ ')
    print('\n')

samples = Zprime+ Zprime_inv + DY+ TT[-1:] +TT[:-1]+ VV+ Wln + Znn

# for year in ['UL16preVFP','UL16postVFP','UL17','UL18']:
#     for sample in samples:
#         print(sample.replace('-HT-','HT')+'_'+year,  str(helper.get_lumi(sample,energy,year,kFactor=True)))

for year in ['UL16preVFP','UL16postVFP','UL17','UL18']:
    for sample in samples:
        info = helper.get_xml(sample,energy,year, info="Source")
        if info=='': continue
        info = info.split('/')
        info = list(filter(lambda x: 'MiniAOD' in x, info))[0]
        info = info.replace('RunIISummer20'+year.replace('preVFP','').replace('postVFP','')+'MiniAOD'+('APV' if year=='UL16preVFP' else '' )+'v2-','').replace('-v1','').replace('-v2','')
        if info=='106X_mcRun2_asymptotic_preVFP_v11': continue
        if info=='106X_mcRun2_asymptotic_v17': continue
        if info=='106X_mc2017_realistic_v9': continue
        if info=='106X_upgrade2018_realistic_v16_L1v1': continue
        print(info, helper.get_xml(sample,energy,year, info="Source").split('/'))