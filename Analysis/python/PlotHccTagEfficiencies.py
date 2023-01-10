from Utils import *
from array import array
import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = 'Simulation'
TDR.extraText2 = 'Work in progress'

# ForThesis(TDR)

ROOT.gInterpreter.ProcessLine('#include "'+os.environ['CMSSW_BASE']+'/src/UHH2/common/include/JetHists.h"')


'''
Module to study Hcc Efficiencies
- Need the Preselection output as input
- Needed for Selection module
- SFs stored as default in ScaleFactors folder. Be carefull since they would override the previous version.
- Runs automatically for all years and channels.
- Default collection is Puppi.
'''

colors = {'UL16preVFP':  ROOT.kAzure-2,
          'UL16postVFP': ROOT.kGreen+2,
          'UL17':        ROOT.kRed+1,
          'UL18':        ROOT.kOrange+1,
          'lepton':      ROOT.kFullCircle,
          'muon':        ROOT.kFullTriangleDown,
          'electron':    ROOT.kFullTriangleUp,
          'invisible':   ROOT.kFullSquare,

}

class PlotHccTagEfficiencies(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.years = ['UL16preVFP', 'UL16postVFP','UL17','UL18']
        TDR.cms_lumi     = self.lumi_map['RunII']['lumiPlot']+' fb^{-1}'
        self.outdir      = self.Path_ANALYSIS+'Analysis/OtherPlots/HccTag/'
        self.nameXaxis   = 'p_{T} [GeV]'
        self.nameYaxis   = 'Efficiency'
        self.Xaxis_min   = 0
        self.Xaxis_max   = 7
        self.Yaxis_min   = 0.001
        self.Yaxis_max   = 1.2
        self.lepton      = 'lepton'
        self.histoname   = 'HccMCEff'
        self.defaultFlav = 'FlavB'
        self.defaultCut  = 'PTMassCut'
        self.Cuts        = ['PTMassCut']
        self.Flavours    = ['Hcc', 'HccMD', 'PN_MD']
        self.Modes       = ['tot', 'Twp', 'Mwp', 'Lwp']
        self.histos      = {}
        self.recreate    = False
        # self.recreate    = True

        os.system('mkdir -p '+self.outdir)
        os.system('mkdir -p '+self.outdir.replace('OtherPlots','ScaleFactors'))
        self.MergeSignalSamples()

    def MergeSignalSamples(self):
        for year in self.years:
            for channel in self.Channels:
                filename = self.Path_STORAGE+year+'/Selection/Puppi/'+channel+'channel/nominal/'+self.PrefixrootFile+'MC.'+self.Signal+'_'+year+'_noTree.root'
                if 'invisible' in channel: filename = filename.replace(self.Signal+'_',self.Signal+'_inv_')
                if self.recreate: os.system('rm '+filename)
                if (not os.path.exists(filename)): os.system('hadd -f '+filename+' '+filename.replace(year+'_noTree.root','M*'+year+'_noTree.root'))

    def LoadHistos(self):
        for year in self.years:
            # Load MC_DY histos
            for channel in self.Channels:
                filename = self.Path_STORAGE+year+'/Selection/Puppi/'+channel+'channel/nominal/'+self.PrefixrootFile+'MC.'+self.Signal+'_'+year+'_noTree.root'
                # filename = self.Path_STORAGE+year+'/Selection/Puppi/'+channel+'channel/nominal/'+self.PrefixrootFile+'MC.'+self.Signal+'_M3000_'+year+'_noTree.root'
                if 'invisible' in channel: filename = filename.replace(self.Signal+'_',self.Signal+'_inv_')
                file_ = ROOT.TFile(filename)
                for cut,flavor,mode in list(itertools.product(self.Cuts, self.Flavours, self.Modes)):
                    hname = year+channel+cut+flavor+mode
                    self.histos[hname] = file_.Get('ZprimeCandidate_'+cut+'/'+flavor+'_ptvsmatch_'+mode).Clone(hname)
                    self.histos[hname].SetDirectory(0)
                    # Merge all the channels
                    hname_lep = hname.replace(channel,self.lepton)
                    if hname_lep in self.histos:
                        self.histos[hname_lep].Add(self.histos[hname])
                    else:
                        self.histos[hname_lep] = self.histos[hname].Clone(hname_lep)
                        self.histos[hname_lep].SetDirectory(0)
                file_.Close()
            # Calculate MC_DY efficiencies
            for channel in self.Channels + [self.lepton]:
                for cut,flavor,mode in list(itertools.product(self.Cuts, self.Flavours, self.Modes)):
                    hname = year+channel+cut+flavor+mode
                    hname_MTwp = hname.replace(mode,'MTwp')
                    if not hname_MTwp+'sum' in self.histos:
                        self.histos[hname_MTwp+'sum'] = ROOT.TH1D(hname_MTwp+'sum',hname_MTwp+'sum', 7,array('d',range(0,8)))
                    self.histos[hname+'sum'] = ROOT.TH1D(hname+'sum',hname+'sum', 7,array('d',range(0,8)))
                    for x in range(1,self.histos[hname].GetNbinsX()+1):
                        val = 0
                        err = 0
                        for y in range(1,self.histos[hname].GetNbinsY()+1):
                            val += self.histos[hname].GetBinContent(x,y)
                            err += self.histos[hname].GetBinError(x,y)*self.histos[hname].GetBinError(x,y)
                        self.histos[hname+'sum'].SetBinContent(x,val)
                        self.histos[hname+'sum'].SetBinError(x,np.sqrt(err))
                        if mode=='Twp':
                            self.histos[hname_MTwp+'sum'].SetBinContent(x,val)
                            self.histos[hname_MTwp+'sum'].SetBinError(x,np.sqrt(err))
                        if mode=='Mwp':
                            self.histos[hname_MTwp+'sum'].SetBinContent(x,val+self.histos[hname_MTwp+'sum'].GetBinContent(x))
                            self.histos[hname_MTwp+'sum'].SetBinError(x,np.sqrt(err+self.histos[hname_MTwp+'sum'].GetBinError(x)*self.histos[hname_MTwp+'sum'].GetBinError(x)))


                    self.histos[hname+'sum'+'ratio'] = self.histos[hname+'sum'].Clone(hname+'sum'+'ratio')
                    self.histos[hname+'sum'+'ratio'].Divide(self.histos[hname+'sum'+'ratio'],self.histos[hname.replace(mode,'tot'+'sum')],1,1,'B')
                    self.histos[hname+'sum'+'ratio'].SetDirectory(0)

                    self.histos[hname+'ratio'] = self.histos[hname].Clone(hname+'ratio')
                    self.histos[hname+'ratio'].Divide(self.histos[hname+'ratio'],self.histos[hname.replace(mode,'tot')],1,1,'B')
                    self.histos[hname+'ratio'].SetDirectory(0)

                    if mode=='Mwp':
                        self.histos[hname_MTwp+'sum'+'ratio'] = self.histos[hname_MTwp+'sum'].Clone(hname_MTwp+'sum'+'ratio')
                        self.histos[hname_MTwp+'sum'+'ratio'].Divide(self.histos[hname_MTwp+'sum'+'ratio'],self.histos[hname_MTwp.replace('MTwp','tot'+'sum')],1,1,'B')
                        self.histos[hname_MTwp+'sum'+'ratio'].SetDirectory(0)

    def ResetCanvas(self, name='canv'):
        self.canv = tdrCanvas(name, self.Xaxis_min, self.Xaxis_max, self.Yaxis_min, self.Yaxis_max, self.nameXaxis,self.nameYaxis)
        # self.leg = tdrLeg(0.15,0.70,0.95,0.89, 0.025, 42, ROOT.kBlack)
        # self.leg.SetNColumns( len(self.Channels) + 1)
        self.leg1 = tdrLeg(0.35,0.70,0.70,0.89, 0.03, 42, ROOT.kBlack)
        self.leg1.SetNColumns(2)
        self.leg2 = tdrLeg(0.70,0.70,0.95,0.89, 0.03, 42, ROOT.kBlack)
        self.leg2.SetNColumns(2)
        self.objs = {}
        for l,c in colors.items():
            obj = ROOT.TGraph()
            self.objs[l] = obj
            if l in self.years:
                obj.SetMarkerColor(c)
                obj.SetLineColor(c)
                self.leg1.AddEntry(obj, l, 'lp')
            else:
                obj.SetMarkerColor(ROOT.kBlack)
                obj.SetMarkerStyle(c)
                self.leg2.AddEntry(obj, l, 'lp')
        GettdrCanvasHist(self.canv).GetXaxis().SetBinLabel(int(GettdrCanvasHist(self.canv).GetNbinsX()*(1-0.5)/7),'200')
        GettdrCanvasHist(self.canv).GetXaxis().SetBinLabel(int(GettdrCanvasHist(self.canv).GetNbinsX()*(2-0.5)/7),'250')
        GettdrCanvasHist(self.canv).GetXaxis().SetBinLabel(int(GettdrCanvasHist(self.canv).GetNbinsX()*(3-0.5)/7),'300')
        GettdrCanvasHist(self.canv).GetXaxis().SetBinLabel(int(GettdrCanvasHist(self.canv).GetNbinsX()*(4-0.5)/7),'350')
        GettdrCanvasHist(self.canv).GetXaxis().SetBinLabel(int(GettdrCanvasHist(self.canv).GetNbinsX()*(5-0.5)/7),'400')
        GettdrCanvasHist(self.canv).GetXaxis().SetBinLabel(int(GettdrCanvasHist(self.canv).GetNbinsX()*(6-0.5)/7),'450')
        GettdrCanvasHist(self.canv).GetXaxis().SetBinLabel(int(GettdrCanvasHist(self.canv).GetNbinsX()*(7-0.5)/7),'500')

    def PlotHistos(self):
        # For Debug
        # for [hname,hist] in self.histos.items():
        #     self.canv = tdrCanvas('canv_'+hname, self.Xaxis_min, self.Xaxis_max, -2.4, 2.4, self.nameXaxis,'eta')
        #     self.canv.SetRightMargin(0.15)
        #     self.canv.SetLogz(True)
        #     hist.SetMinimum(self.Yaxis_min)
        #     hist.SetMaximum(self.Yaxis_max)
        #     hist.Draw('colz')
        #     self.canv.SaveAs(self.outdir+hname+'.pdf')
        for match in ['sum','HbbMatch','HWWMatch','HccMatch']:
            for flavor,mode in list(itertools.product(self.Flavours, self.Modes+['MTwp'])):
                if mode=='tot': continue
                if match!='sum' and mode=='MTwp': continue
                self.ResetCanvas(match+flavor+mode)
                for year in self.years:
                    for channel in self.Channels + [self.lepton]:
                        hname_ = year+channel+self.defaultCut+flavor+mode+'ratio'
                        if match == 'sum':
                            hname_ = hname_.replace('ratio', 'sumratio')
                        hname = hname_+match
                        self.histos[hname] = ROOT.TH1D(hname,hname, 7,array('d',range(0,8)))
                        for x in range(1,self.histos[hname].GetNbinsX()+1):
                            if match == 'sum':
                                self.histos[hname].SetBinContent(x,self.histos[hname_].GetBinContent(x))
                                self.histos[hname].SetBinError(x,self.histos[hname_].GetBinError(x))
                            else:
                                self.histos[hname].SetBinContent(x,self.histos[hname_].GetBinContent(x,ROOT.StringToMatching(match)+1))
                                self.histos[hname].SetBinError(x,self.histos[hname_].GetBinError(x,ROOT.StringToMatching(match)+1))
                        for col in colors:
                            if col in hname:
                                if 'UL' in col: color = colors[col]
                                else : point = colors[col]
                        tdrDraw(self.histos[hname], 'P', point, color, 1, color, 0, color)
                        # self.leg.AddEntry(self.histos[hname], hname.replace(self.defaultCut,'').replace('sum','').replace('ratio','').replace(flavor+mode,(flavor+mode).replace('Hcc','')).replace('Match',''), 'lp')
                self.canv.SaveAs(self.outdir+'Years_lepton_'+flavor+'_'+mode+'_'+match+'_Signal.pdf')

    def SaveRootFiles(self):
        for year in self.years:
            file_ = ROOT.TFile(self.outdir.replace('OtherPlots','ScaleFactors')+'SF_'+year+'.root', 'RECREATE')
            for flavor in self.Flavours:
                self.histos[year+self.lepton+self.defaultCut+flavor].Write(self.histoname+flavor+'Eff')
            file_.Close()


def main():
    PlotSyst = PlotHccTagEfficiencies()
    PlotSyst.LoadHistos()
    PlotSyst.PlotHistos()
    # PlotSyst.SaveRootFiles()


if __name__ == '__main__':
    main()
