from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = 'Simulation'
TDR.extraText2 = 'Work in progress'

# ForThesis(TDR)

'''
Module to plot Tagger shapes
- Need the Selectionas input
'''

class PlotTaggerShapes(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.Samples = ['DY', 'WJets', 'ZprimeToZH']
        # self.Samples = ['ZprimeToZH_M2000']
        # self.Samples = ['DY', 'ZprimeToZH_M2000']
        # self.Discriminators = ['H4qvsQCD','H4qvsQCD_MD','HccvsQCD','HccvsQCD_MD','ZHccvsQCD','ZHccvsQCD_MD']
        self.Discriminators = ['DeepBoosted_H4qvsQCD','DeepBoosted_ZHccvsQCD_MD', 'ParticleNet_XccvsQCD_MD']

        self.color  = {'Background': ROOT.kOrange+1,
                       'Signal':     ROOT.kAzure+2,
                       }
        self.outdir = self.Path_ANALYSIS+'Analysis/OtherPlots/TaggerShape/'
        os.system('mkdir -p '+self.outdir)

    def LoadHistos(self):
        self.hists = {}
        for sample in ['Signal','Background']:
            for channel in self.Channels+['lepton']:
                for disc in self.Discriminators:
                    name_ = sample+channel+disc
                    self.hists[name_] = ROOT.TH1D(name_,name_,30, -0.01, 1.01)
                    self.hists[name_].SetDirectory(0)
        for sample in self.Samples:
            for channel in self.Channels:
                for year in self.years:
                    for collection in self.Collections:
                        filename = self.Path_STORAGE+year+'/Selection/'+collection+'/'+channel+'channel/nominal/work*/uhh2.AnalysisModuleRunner.MC.MC_'+sample+'*.root'
                        if channel == 'invisible':
                            filename = filename.replace('ZprimeToZH', 'ZprimeToZH_inv')
                        elif sample=='WJets': continue
                        for fname in glob.glob(filename):
                            print fname
                            f_ = ROOT.TFile(fname)
                            t_ = f_.Get('AnalysisTree')
                            for ev in t_:
                                if ev.ZprimeCandidate.size()!=1 :
                                    continue
                                for zp in ev.ZprimeCandidate:
                                    for disc in self.Discriminators:
                                        name = ('Signal' if 'ZprimeToZH' in sample else 'Background')+channel+disc
                                        val = zp.discriminator('btag_'+disc)
                                        self.hists[name].Fill(val,ev.weight_GLP)
                                        if not 'inv' in channel:
                                            self.hists[name.replace(channel,'lepton')].Fill(val,ev.weight_GLP)
                            f_.Close()

    def PlotHistos(self):
        self.LoadHistos()
        for disc in self.Discriminators:
            for channel in self.Channels+['lepton']:
                TDR.cms_lumi = self.lumi_map['RunII']['lumiPlot']+' fb^{-1}' if TDR.extraText!='Simulation' else 'MC RunII'
                canv = tdrCanvas(disc+channel, -0.05, 1.05, 0.0001, 20, disc.replace('_MD','').replace('DeepBoosted', 'deepak8'),'A.U.')
                canv.SetLogy(1)
                leg = tdrLeg(0.70,0.70,0.89,0.89, 0.030, 42, ROOT.kBlack);
                for sample in ['Signal', 'Background']:
                    color = self.color[sample]
                    h_ = self.hists[sample+channel+disc]
                    h_.Scale(1./h_.Integral())
                    h_.SetLineWidth(2)
                    tdrDraw(h_, 'hist', 0, color, 1, color, 1001, color)
                    h_.SetFillColorAlpha(color, 0.6 if 'Background' in sample else 0.45)
                    leg.AddEntry(h_, sample, 'f')
                canv.SaveAs(self.outdir+'TaggerShape_'+disc+'_'+channel+'.pdf')



if __name__ == '__main__':
    PlotBkg = PlotTaggerShapes()
    PlotBkg.PlotHistos()
