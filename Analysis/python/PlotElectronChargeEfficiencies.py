from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = 'Simulation'
TDR.extraText2 = 'Work in progress'

'''
Module to plot Signal Efficiencies
- Need the Selectionas input
'''

class PlotElectronChargeEfficiencies(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.MassPoints = np.array(self.MassPoints)
        self.SampleMask = (self.MassPoints>=1000)*(self.MassPoints<=6000)
        self.SampleMask = (self.MassPoints>=600)*(self.MassPoints<=8000)
        self.SamplesPlot = self.MassPoints[self.SampleMask]
        self.Samples = ['M'+str(m) for m in self.SamplesPlot]
        self.years = self.years+['RunII']
        self.years = ['UL18']
        self.Channels = ['electron']

        self.outdir = self.Path_ANALYSIS+'Analysis/OtherPlots/ElectronChargeEfficiencies/'
        os.system('mkdir -p '+self.outdir)
        self.LoadValues()

    def LoadValues(self):
        infos = OrderedDict([
            # ('weights',    {'color': ROOT.kBlack,   'legend': 'Inclusive'}),
            ('NLeptonSel', {'color': ROOT.kAzure+2, 'legend': 'Selection: 2 lep + kin. cuts'}),
            ('Ele_charge', {'color': ROOT.kRed+1,   'legend': 'Selection: opposite charge'}),
        ])
        for year in self.years:
            for collection in self.Collections:
                for channel in self.Channels:
                        vals = OrderedDict()
                        for folder in infos.keys():
                            vals[folder] = []
                        for sample in self.Samples:
                            if DoControl([''], year+channel, channel, sample): continue
                            fname = self.Signal+('_inv' if 'inv' in channel else '')+'_'+sample+'_'+year
                            fname = self.Path_SFRAME+year+'/LeptonIDStudies/'+collection+'/'+channel+'channel/nominal/workdir_LeptonIDStudies_'+fname+'/uhh2.AnalysisModuleRunner.MC.'+fname+'_0.root'
                            f_ = ROOT.TFile(fname)
                            for folder in infos.keys():
                                h_ = f_.Get('ZprimeCandidate_'+folder+'/sum_event_weights')
                                vals[folder].append(h_.GetBinContent(1))
                            f_.Close()
                        self.canv = tdrCanvas(year+collection+channel+'Cuts', self.SamplesPlot[0]-200, self.SamplesPlot[-1]+200, 0, 1.4, "M(Z') (GeV)", 'Signal selection efficiency')
                        self.leg = tdrLeg(0.40, 0.70, 0.89, 0.89, 0.040)
                        grs = []
                        for folder, info in infos.items():
                            eff = np.array(vals[folder])
                            eff /= np.array(vals[infos.keys()[0]])
                            gr = ROOT.TGraphErrors(len(self.SamplesPlot), array('d',self.SamplesPlot), array('d',eff))
                            gr.SetLineWidth(2)
                            color = info['color']
                            tdrDraw(gr, 'lp', ROOT.kFullDotLarge, mcolor=color, lcolor=color, fstyle=0)
                            self.leg.AddEntry(gr, info['legend'],'lp')
                            grs.append(gr)
                        self.canv.SaveAs(self.outdir+'SignalEfficiencies_'+year+'_'+collection+'_'+channel+'.pdf')

if __name__ == '__main__':
    PlotBkg = PlotElectronChargeEfficiencies()
    # PlotBkg.PlotEfficiencyCuts()