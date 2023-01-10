from Utils import *
from array import array
import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = 'Simulation'
TDR.extraText2 = 'Work in progress'

# ForThesis(TDR)

ROOT.gInterpreter.ProcessLine('#include "'+os.environ['CMSSW_BASE']+'/src/UHH2/common/include/JetHists.h"')


'''
Module to study BTag Efficiencies
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

class PlotBTagEfficiencies(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.years = ['UL16preVFP', 'UL16postVFP','UL17','UL18']
        TDR.cms_lumi     = self.lumi_map['RunII']['lumiPlot']+' fb^{-1}'
        self.outdir      = self.Path_ANALYSIS+'Analysis/OtherPlots/BTag/'
        self.nameXaxis   = 'p_{T} [GeV]'
        self.nameYaxis   = 'Efficiency'
        self.Xaxis_min   = 0
        self.Xaxis_max   = 1000-1
        self.Yaxis_min   = 0.001
        self.Yaxis_max   = 1.2
        self.lepton      = 'lepton'
        self.histoname   = 'BTagMCEff'
        self.defaultFlav = 'FlavB'
        self.defaultCut  = 'JetDiLeptonPhiAngular'
        self.Cuts        = ['NLeptonSel', 'NBoostedJet','JetDiLeptonPhiAngular']
        self.Flavours    = ['FlavB', 'FlavC', 'FlavUDSG']
        self.Modes       = ['Passing', 'Total']
        self.histos      = {}
        self.recreate    = False
        # self.recreate    = True

        os.system('mkdir -p '+self.outdir)
        os.system('mkdir -p '+self.outdir.replace('OtherPlots','ScaleFactors'))
        self.MergeSignalSamples()

    def MergeSignalSamples(self):
        for year in self.years:
            for channel in self.Channels:
                filename = self.Path_STORAGE+year+'/Preselection/Puppi/'+channel+'channel/nominal/'+self.PrefixrootFile+'MC.'+self.Signal+'_'+year+'_noTree.root'
                if 'invisible' in channel: filename = filename.replace(self.Signal+'_',self.Signal+'_inv_')
                if self.recreate: os.system('rm '+filename)
                if (not os.path.exists(filename)): os.system('hadd -f '+filename+' '+filename.replace(year+'_noTree.root','M*'+year+'_noTree.root'))

    def LoadHistos(self):
        for year in self.years:
            # Load MC_DY histos
            for channel in self.Channels:
                filename = self.Path_STORAGE+year+'/Preselection/Puppi/'+channel+'channel/nominal/'+self.PrefixrootFile+'MC.'+self.Signal+'_'+year+'_noTree.root'
                if 'invisible' in channel: filename = filename.replace(self.Signal+'_',self.Signal+'_inv_')
                file_ = ROOT.TFile(filename)
                for cut,flavor,mode in list(itertools.product(self.Cuts, self.Flavours, self.Modes)):
                    hname = year+channel+cut+flavor+mode
                    self.histos[hname] = file_.Get(self.histoname.replace('MC','')+'_'+cut+'/'+self.histoname+flavor+mode).Clone(hname)
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
                for cut,flavor in list(itertools.product(self.Cuts, self.Flavours)):
                        hname = year+channel+cut+flavor
                        self.histos[hname] = self.histos[hname+'Passing'].Clone(hname)
                        self.histos[hname].Divide(self.histos[hname+'Total'])
                        self.histos[hname].SetDirectory(0)

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
                # self.leg.AddEntry(self.histos[hname+'pt'], hname.replace(self.defaultCut,'').replace('FlavUDSG','_light'), 'lp')

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
        for flavor in self.Flavours:
            self.ResetCanvas(flavor)
            for year in self.years:
                for channel in self.Channels + [self.lepton]:
                    hname = year+channel+self.defaultCut+flavor
                    self.histos[hname+'pt'] = ROOT.TH1D(hname+'pt',hname+'pt', ROOT.kBTagMCEffBinsPt.size()-1,array('d',list(ROOT.kBTagMCEffBinsPt)))
                    for x in range(1,self.histos[hname+'pt'].GetNbinsX()+1):
                        self.histos[hname+'pt'].SetBinContent(x,self.histos[hname].GetBinContent(x,1))
                        self.histos[hname+'pt'].SetBinError(x,self.histos[hname].GetBinError(x,1))
                    for col in colors:
                        if col in hname:
                            if 'UL' in col: color = colors[col]
                            else : point = colors[col]
                    tdrDraw(self.histos[hname+'pt'], 'P', point, color, 1, color, 0, color)
                    # self.leg.AddEntry(self.histos[hname+'pt'], hname.replace(self.defaultCut,'').replace('FlavUDSG','_light'), 'lp')
            self.canv.SaveAs(self.outdir+'Years_lepton_'+flavor+'_Signal.pdf')

    def SaveRootFiles(self):
        for year in self.years:
            file_ = ROOT.TFile(self.outdir.replace('OtherPlots','ScaleFactors')+'SF_'+year+'.root', 'RECREATE')
            for flavor in self.Flavours:
                self.histos[year+self.lepton+self.defaultCut+flavor].Write(self.histoname+flavor+'Eff')
            file_.Close()


def main():
    PlotSyst = PlotBTagEfficiencies()
    PlotSyst.LoadHistos()
    PlotSyst.PlotHistos()
    PlotSyst.SaveRootFiles()


if __name__ == '__main__':
    main()
