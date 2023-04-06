from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = 'Simulation'
TDR.extraText2 = 'Work in progress'

'''
Module to plot Signal Efficiencies
- Need the Selectionas input
'''

class PlotSignalEfficiencies(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.MassPoints = np.array(self.MassPoints)
        self.Samples = ['M'+str(m) for m in self.MassPoints]
        self.SampleMask = (self.MassPoints>=1000)*(self.MassPoints<=6000)
        self.SamplesPlot = self.MassPoints[self.SampleMask]
        self.years = self.years+['RunII']
        self.years = ['RunII']
        # self.Channels = ['invisible']

        self.outdir = self.Path_ANALYSIS+'Analysis/OtherPlots/SignalEfficiencies/'
        os.system('mkdir -p '+self.outdir)

        self.colors = {
            'Inc'    : (ROOT.kBlack,    'H #rightarrow inclusive'),
            'bb'     : (ROOT.kAzure+2,  'H #rightarrow b#bar{b}'),
            'WW'     : (ROOT.kGreen+2,  'H #rightarrow VV* (4q merged)'),
            'cc'     : (ROOT.kGreen+3,  'H #rightarrow c#bar{c}'),
            'gg'     : (ROOT.kOrange+1, 'H #rightarrow gg'),
            'else'   : (ROOT.kOrange-2, 'H #rightarrow VV* (other)'),
            'tautau' : (ROOT.kRed+1,    'H #rightarrow #tau^{+}#tau^{#font[122]{\55}}'),
        }

        self.decays = []
        self.Hdecays = ['bb', 'cc', 'WW', 'tautau', 'else', 'gg', 'Inc']
        self.Zdecays = ['ee', 'mumu', 'else']
        for H in self.Hdecays:
            self.decays.append(('Hto'+H).replace('toelse', 'else'))
        for Z in ['ee', 'mumu', 'else']:
            self.decays.append('Z'+Z)
        for H in self.Hdecays:
            if 'Inc' in H: continue
            for Z in ['ee', 'mumu', 'else']:
                self.decays.append(('Z'+Z+'Hto'+H).replace('toelse', 'else'))
                self.decays.append('Z'+Z)

        self.CutsList = ['Triggers', 'Lepton selection', 'Angular cuts', 'b tag veto', 'ParticleNet', 'DeepAk8']
        self.Cuts = {
            'weights': {
                'module'   : 'Preselection',
                'color'    : ROOT.kBlack,
                'style'    : ROOT.kSolid,
                'folder'   : 'weights',
                },
            'Triggers': {
                'module'   : 'Preselection',
                'color'    : ROOT.kRed+1,
                'style'    : ROOT.kSolid,
                'folder'   : 'Trigger',
                },
            'Lepton selection': {
                'module'   : 'Preselection',
                'color'    : ROOT.kOrange+1,
                'style'    : ROOT.kSolid,
                'folder'   : 'NLeptonSel',
                },
            'Angular cuts': {
                'module'   : 'Selection',
                'color'    : ROOT.kOrange-2,
                'style'    : ROOT.kSolid,
                'folder'   : 'MuonScale',
                },
            'b tag veto': {
                'module'   : 'Selection',
                'color'    : ROOT.kAzure+2,
                'style'    : ROOT.kSolid,
                'folder'   : 'ZprimeSelection',
                },
            'PTMassCut': {
                'module'   : 'Selection',
                'color'    : ROOT.kBlack,
                'style'    : ROOT.kSolid,
                'folder'   : 'PTMassCut',
                },
            'ScaleFactors': {
                'module'   : 'Selection',
                'color'    : ROOT.kBlack,
                'style'    : ROOT.kSolid,
                'folder'   : 'ScaleFactors',
                },
            'ParticleNet': {
                'module'   : 'SignalRegion',
                'color'    : ROOT.kGreen+3,
                'style'    : ROOT.kSolid,
                'folder'   : 'PN_ZHccvsQCD_MD_SR',
                },
            'DeepAk8': {
                'module'   : 'SignalRegion',
                'color'    : ROOT.kGreen+2,
                'style'    : ROOT.kDashed,
                'folder'   : 'DeepAk8_ZHccvsQCD_MD_SR',
                },
        }

        self.values = {}
        self.LoadValues()

    def LoadValues(self):
        for year in self.years:
            for collection in self.Collections:
                for channel in self.Channels:
                    for cut in self.Cuts:
                        module = self.Cuts[cut]['module']
                        folder = self.Cuts[cut]['folder']
                        for decay in self.decays:
                            unique_name = year+collection+channel+folder+decay
                            if not unique_name in self.values:
                                self.values[unique_name] = np.array([0.]*len(self.Samples))
                            unique_name = year+collection+'chargedlepton'+folder+decay
                            if not unique_name in self.values:
                                self.values[unique_name] = np.array([0.]*len(self.Samples))
                        for index, sample in enumerate(self.Samples):
                            if DoControl([''], year+channel, channel, sample): continue
                            fname = self.Path_STORAGE+year+'/'+module+'/'+collection+'/'+channel+'channel/nominal/uhh2.AnalysisModuleRunner.MC.'+self.Signal+('_inv' if 'inv' in channel else '')+'_'+sample+'_'+year+'_noTree.root'
                            f_ = ROOT.TFile(fname)
                            for decay in self.decays:
                                unique_name = year+collection+channel+folder+decay
                                h_ = f_.Get('ZprimeCandidate_'+folder+'/sum_event_weights'+('_'+decay if not 'Inc' in decay else ''))
                                val = h_.GetBinContent(1)
                                self.values[unique_name][index] += val
                                if not 'inv' in channel:
                                    if cut=='weights' and 'muon' in channel: continue
                                    self.values[unique_name.replace(channel,'chargedlepton')][index] += val
                            f_.Close()

    def CreateCanvas(self, CanvName, isLow=False, isBR=False, isSR=False):
        ymin, ymax = (0.005,40)
        if isLow:
            ymin = 0.0001 if isSR else 0.001
        if isBR:
            ymin,ymax = (0.0,1.2)
        self.canv = tdrCanvas(CanvName, self.SamplesPlot[0]-200, self.SamplesPlot[-1]+200, ymin, ymax, "M(Z') (GeV)", 'Signal selection efficiency')
        self.canv.SetLogy(not isBR)
        lowLeg = 'HiggsDecays' in CanvName and not isBR
        self.leg = tdrLeg(0.40,0.65 if lowLeg else 0.70, 0.93 if lowLeg else 0.95,0.90, 0.040, 42, ROOT.kBlack)
        self.leg.SetNColumns(2)

    def PlotEfficiencyCuts(self):
        plotName = 'Cuts_'
        for year in self.years:
            for collection in self.Collections:
                for channel in self.Channels+['chargedlepton']:
                    TDR.cms_lumi = self.lumi_map[year]['lumiPlot']+' fb^{-1}' if TDR.extraText!='Simulation' else 'MC '+year
                    self.CreateCanvas(year+collection+channel+'Cuts')
                    grs = []
                    for cut in self.CutsList:
                        folder  = self.Cuts[cut]['folder']
                        color   = self.Cuts[cut]['color']
                        style   = self.Cuts[cut]['style']
                        eff = self.values[year+collection+channel+folder+'HtoInc'][self.SampleMask]
                        eff /= self.values[year+collection+channel+'weights'+'HtoInc'][self.SampleMask]
                        if 'SignalRegion' == self.Cuts[cut]['module']:
                            eff *= self.values[year+collection+channel+'PTMassCut'+'HtoInc'][self.SampleMask]
                            eff /= self.values[year+collection+channel+'ScaleFactors'+'HtoInc'][self.SampleMask]
                        gr = ROOT.TGraphErrors(len(self.SamplesPlot), array('d',self.SamplesPlot), array('d',eff))
                        gr.SetLineWidth(2)
                        tdrDraw(gr, 'lp', ROOT.kFullDotLarge, color, style, color, 0, color);
                        self.leg.AddEntry(gr, cut,'lp')
                        grs.append(gr)
                    self.canv.SaveAs(self.outdir+'SignalEfficiencies_'+plotName+year+'_'+collection+'_'+channel+'.pdf')

    def PlotEfficiencyHiggsDecays(self, isBR=False, isSR=False):
        plotName = 'HiggsDecays_'+('BR_' if isBR else '')+('SR_' if isSR else '')
        cut  = 'PTMassCut' if not isSR else 'ParticleNet'
        folder  = self.Cuts[cut]['folder']
        for year in self.years:
            for collection in self.Collections:
                for channel in self.Channels:
                    self.CreateCanvas(year+collection+channel+plotName, isLow= True, isBR = isBR, isSR = isSR)
                    grs = []
                    for decay in self.Hdecays:
                        if 'Inc' in decay and isBR: continue
                        color   = self.colors[decay][0]
                        legName = self.colors[decay][1]
                        decay = ('Hto' if not 'else' in decay else 'H')+decay
                        decaynorm = decay if isBR else 'HtoInc'
                        if isBR and not 'Inc' in decay:
                            decaynorm = ('Zee' if 'ele' in channel else ('Zmumu' if 'muon' in channel else '') ) +decaynorm
                        eff = self.values[year+collection+channel+folder+decay][self.SampleMask]
                        eff /= self.values[year+collection+channel+'weights'+decaynorm][self.SampleMask]
                        gr = ROOT.TGraphErrors(len(self.SamplesPlot), array('d',self.SamplesPlot), array('d',eff))
                        gr.SetLineWidth(2)
                        tdrDraw(gr, 'lp', ROOT.kFullDotLarge, color, ROOT.kSolid, color, 0, color);
                        self.leg.AddEntry(gr, legName,'lp')
                        grs.append(gr)
                    self.canv.SaveAs(self.outdir+'SignalEfficiencies_'+plotName+year+'_'+collection+'_'+channel+'.pdf')



if __name__ == '__main__':
    PlotBkg = PlotSignalEfficiencies()
    PlotBkg.PlotEfficiencyCuts()
    PlotBkg.PlotEfficiencyHiggsDecays(isBR = False,isSR = False)
    PlotBkg.PlotEfficiencyHiggsDecays(isBR = False,isSR = True)
    PlotBkg.PlotEfficiencyHiggsDecays(isBR = True, isSR = False)
    PlotBkg.PlotEfficiencyHiggsDecays(isBR = True, isSR = True)
