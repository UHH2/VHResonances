from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = 'Simulation'
TDR.extraText2 = 'Work in progress'

'''
Module to plot Tagger shapes
- Need the Selectionas input
'''


class FitSignalShapes(VariablesBase):
    def __init__(self, xmin=-1, xmax=-1):
        VariablesBase.__init__(self)
        self.xmin = xmin
        self.xmax = xmax
        self.years = ['RunII']
        # self.Channels = ['muon']
        # self.Channels = ['electron']
        # self.Channels = ['invisible']
        # self.Folder = 'DeepAk8_ZHccvsQCD_MD_SR'
        self.Folder = 'PN_ZHccvsQCD_MD_SR'
        self.Syst_list = self.Systematics[1:]+self.Systematics_Scale
        self.Syst = [self.Systematics[0]]+[x+v for x in self.Syst_list for v in ['Up', 'Down'] if (not 'murmuf' in x and not 'PDF' in x)]
        self.Syst_list = []
        self.Syst = ['nominal']
        # self.SignalSamples = [self.Signal+mode+'_M'+str(mass) for mass in [1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000] for mode in ['','_inv']]
        self.SignalSamples = [self.Signal+mode+'_M'+str(mass) for mass in [3000] for mode in ['','_inv']]

        self.color  = {'EGE': ROOT.kAzure+2,
                       'CB': ROOT.kOrange+1,
                       }

        self.outdir = self.Path_ANALYSIS+'Analysis/OtherPlots/SignalShapes/'
        os.system('mkdir -p '+self.outdir)
        self.LoadHistos()

    def LoadHistos(self):
        self.hists = {}
        for year in self.years:
            for channel in self.Channels:
                for collection in self.Collections:
                    for sample in self.SignalSamples:
                        for sys in self.Syst:
                            unique_name = year+'_'+channel+'_'+collection+'_'+sample+'_'+sys
                            if DoControl([''], unique_name, channel, sample): continue
                            folder = sys.replace('Up','_up').replace('Down','_down') if sys.replace('Up','').replace('_Down','') in self.Systematics else 'nominal'
                            fname = self.Path_STORAGE+'/'+year+'/SignalRegion/'+collection+'/'+channel+'channel/'+folder+'/'+self.PrefixrootFile+'MC.'+sample+'_'+year+'_noTree.root'
                            f_ = ROOT.TFile(fname)
                            self.hists[unique_name] = f_.Get('ZprimeCandidate_'+self.Folder+'/Zprime_mass'+('' if not 'inv' in channel else '_transversal')+'_rebin100')
                            self.hists[unique_name].SetDirectory(0)
                            self.hists[unique_name].Scale(self.BRs[channel]*round(ROOT.xsec_ref['default_value'],3))
                            if 'inv' in channel:
                                for i in range(self.hists[unique_name].GetNbinsX()):
                                    self.hists[unique_name].SetBinError(i,self.hists[unique_name].GetBinError(i)*1.15)
                            f_.Close()

    def DoFits(self):
        self.funcs = {}
        self.chi2 = {}
        self.errors = {}
        for year in self.years:
            for channel in self.Channels:
                isInv = 'inv' in channel
                isEle = 'ele' in channel
                for collection in self.Collections:
                    for sample in self.SignalSamples:
                        if DoControl([''], year+'_'+channel+'_'+collection+'_'+sample, channel, sample): continue
                        mass = float(sample.replace('MC_ZprimeToZH_M','').replace('MC_ZprimeToZH_inv_M',''))
                        p_xmin = mass-20*(mass*0.02+20.)
                        p_xmax =  mass+20*(mass*0.02+20.)
                        f_xmin, f_xmax = p_xmin, p_xmax
                        if isInv :
                            f_xmin = 0.7742*mass-113;
                            f_xmax = 1.0857*mass+46;
                            if (mass==1000): f_xmin, f_xmax = (1000., 1300.)
                            if (mass==1200): f_xmin, f_xmax = (1000., 1700.)
                            if (mass==1400): f_xmin, f_xmax = (1000., 1900.)
                            if (mass==1600): f_xmin, f_xmax = (1000., 2100.)
                            if (mass==1800): f_xmin, f_xmax = (1000., 2500.)
                            if (mass==2000): f_xmin, f_xmax = (1200., 2600.)
                            if (mass==2500): f_xmin, f_xmax = (1000., 3400.)
                            if (mass==3000): f_xmin, f_xmax = (1500., 3500.) #OK
                            if (mass==3500): f_xmin, f_xmax = (1100., 3900.)
                            if (mass==4000): f_xmin, f_xmax = (1800., 4800.)
                            if (mass==4500): f_xmin, f_xmax = (1300., 4900.)
                            if (mass==5000): f_xmin, f_xmax = (1400., 5900.)
                            if (mass==5500): f_xmin, f_xmax = (2500., 6000.)
                            if (mass==6000): f_xmin, f_xmax = (2800., 7000.)
                            if (mass==7000): f_xmin, f_xmax = (3000., 7600.)
                            if (mass==8000): f_xmin, f_xmax = (4000., 8800.)
                        else:
                            f_xmin = 0.7742*mass-113;
                            f_xmax = 1.0857*mass+46;

                            if (mass==1000): f_xmin, f_xmax = (700, 1400)
                            if (mass==1200): f_xmin, f_xmax = (700, 1600)
                            if (mass==1400): f_xmin, f_xmax = (1000, 1600 if isEle else 1800)
                            if (mass==1600): f_xmin, f_xmax = (1000, 1900)
                            if (mass==1800): f_xmin, f_xmax = (1100, 2000)
                            if (mass==2000): f_xmin, f_xmax = (1200, 2200)
                            if (mass==2500): f_xmin, f_xmax = (1800, 2700)
                            if (mass==3000): f_xmin, f_xmax = (1700, 3300)#OK
                            if (mass==3500): f_xmin, f_xmax = (2400, 3900)
                            if (mass==4000): f_xmin, f_xmax = (3000, 4400 if isEle else 4500)
                            if (mass==4500): f_xmin, f_xmax = (3400, 4900)
                            if (mass==5000): f_xmin, f_xmax = (4300 if isEle else 3300, 5300 if isEle else 5400)
                            if (mass==5500): f_xmin, f_xmax = (4400, 6000)
                            if (mass==6000): f_xmin, f_xmax = (4600, 6600)
                            if (mass==7000): f_xmin, f_xmax = (5700, 7600)
                            if (mass==8000): f_xmin, f_xmax = (6800, 8800)
                        if self.xmin!=-1 or self.xmax!=-1:
                            f_xmin, f_xmax = (self.xmin, self.xmax)
                        for sys in self.Syst:
                            unique_name = year+'_'+channel+'_'+collection+'_'+sample+'_'+sys
                            if DoControl([''], unique_name, channel, sample): continue
                            hist = self.hists[unique_name]
                            self.funcs[unique_name+'CB']  = ROOT.TF1(unique_name+'CB',  CrystalBall_Fit, f_xmin, f_xmax, 5)
                            self.funcs[unique_name+'EGE'] = ROOT.TF1(unique_name+'EGE', ExpGaussExp_Fit, f_xmin, f_xmax, 5)
                            self.funcs[unique_name+'CB'].SetParameters(mass,100.,0.7,5.,hist.Integral())
                            self.funcs[unique_name+'CB'].SetParLimits(0, mass*0.7, mass*1.05)
                            self.funcs[unique_name+'CB'].SetParLimits(1, 10., 500)
                            self.funcs[unique_name+'CB'].SetParLimits(2, 0, 10)
                            self.funcs[unique_name+'CB'].SetParLimits(3, 0, 30)
                            self.funcs[unique_name+'CB'].SetParLimits(4, 0.001, 100)

                            self.funcs[unique_name+'EGE'].SetParameters(mass,100.,2.,2.,20.)
                            self.funcs[unique_name+'EGE'].SetParLimits(0, mass*0.7, mass*1.05)
                            self.funcs[unique_name+'EGE'].SetParLimits(1, 10., 500)
                            self.funcs[unique_name+'EGE'].SetParLimits(2, 0.001, 10)
                            self.funcs[unique_name+'EGE'].SetParLimits(3, 0, 10)
                            self.funcs[unique_name+'EGE'].SetParLimits(4, 0.001, 100)

                            for func_name in ['CB','EGE']:
                                fname = unique_name+func_name
                                func = self.funcs[fname]
                                func.SetLineColor(self.color[func_name])
                                fitRes = hist.Fit(fname, 'RQMS')
                                self.errors['band'+fname], self.errors['band_pull'+fname], self.errors['pull'+fname] = ComputeHistWithCL(unique_name+func_name, func, fitRes, hist, cl=0.95)
                                self.chi2[fname] = (func.GetChisquare(), func.GetNDF(), func.GetProb())
                        TDR.lumi_13TeV  = str(round(float(self.lumi_map[year]['lumi_fb']),1))+' fb^{-1}' if TDR.extraText!='Simulation' else 'MC '+year
                        ymax = 0.3 if not isInv else 1
                        # ymax = 10 if not isInv else 50
                        canv = tdrDiCanvas(year+'_'+channel+'_'+collection+'_'+sample, p_xmin, p_xmax, 0.001, ymax, 0.5, 1.5, "M(Z') [GeV]" if not isInv else "M_{T}(Z') [GeV]", 'Events', 'Hist/Fit')
                        rt.gStyle.SetOptFit(0)
                        canv.cd(1)
                        # .SetLogy(1)
                        # leg = tdrLeg(0.40,0.60,0.92,0.89, 0.035, 42, ROOT.kBlack)
                        leg = {}
                        leg['CB']  = tdrLeg(0.60,0.89-0.045*8,0.92,0.89, 0.04, 42, ROOT.kBlack)
                        leg['EGE'] = tdrLeg(0.60,0.50-0.045*8,0.92,0.50, 0.04, 42, ROOT.kBlack)
                        for sys in ['nominal']:
                            unique_name = year+'_'+channel+'_'+collection+'_'+sample+'_'+sys
                            if DoControl([''], unique_name, channel, sample): continue
                            hist = self.hists[unique_name]
                            for func_name in ['CB','EGE']:
                                isCB = 'CB' == func_name
                                fname = unique_name+func_name
                                color = self.color[func_name]
                                func = self.funcs[fname]
                                func.SetLineWidth(2)
                                func.SetNpx(1000)
                                canv.cd(1)
                                band = self.errors['band'+fname]
                                chi2_red = self.chi2[fname][0]/self.chi2[fname][1]
                                prob_new = self.chi2[fname][2]
                                legName = 'CrystalBall' if isCB else 'ExpGaussExp'
                                leg[func_name].AddEntry(func, legName, 'l')
                                leg[func_name].AddEntry(ROOT.TObject(), '#chi^{2}/n.d.f = '+str(int(chi2_red*100)/100.), '')
                                leg[func_name].AddEntry(ROOT.TObject(), 'p-value = '+str(int(prob_new*100))+' %', '')
                                legName = 'Norm = '+ str(round(func.GetParameter(4),2))
                                leg[func_name].AddEntry(ROOT.TObject(), legName, '')
                                for par in range(0,4):
                                    if par == 0: parName = '#mu'
                                    if par == 1: parName = '#sigma'
                                    if par == 2: parName = '#alpha' if isCB else 'k_{L}'
                                    if par == 3: parName = 'n' if isCB else 'k_{H}'
                                    legName = parName+' = '+ str(round(func.GetParameter(par),2))+' #pm '+str(round(func.GetParError(par),2))
                                    leg[func_name].AddEntry(ROOT.TObject(), legName, '')
                                tdrDraw(band, 'e3', 9, ROOT.kWhite, 1, color+1, 1001, color+1)
                                band.SetMarkerSize(0)
                                band.SetFillColorAlpha(color+1,0.35)
                                func.Draw('same')
                            for func_name in ['EGE','CB']:
                                isCB = 'CB' == func_name
                                fname = unique_name+func_name
                                color = self.color[func_name]
                                pull = self.errors['pull'+fname]
                                band_pull = self.errors['band_pull'+fname]
                                canv.cd(2)
                                tdrDraw(pull, '', ROOT.kFullCircle, color, 1, color, 0, color)
                                tdrDraw(band_pull, 'e3', 9, ROOT.kWhite, 1, color+1, 1001, color+1)
                                band_pull.SetMarkerSize(0)
                                band_pull.SetFillColorAlpha(color+1,0.35)
                                line_ =  rt.TLine(p_xmin, 1, p_xmax, 1)
                                line_.SetLineColor(ROOT.kBlack)
                                line_.SetLineStyle(ROOT.kSolid)
                                line_.Draw('same')
                            canv.cd(1)
                            tdrDraw(hist, '', ROOT.kFullCircle, ROOT.kBlack, 1, ROOT.kBlack, 0, ROOT.kBlack)
                        if self.xmin!=-1 or self.xmax!=-1:
                            canv.SaveAs(self.outdir+'SignalShape_'+year+'_'+channel+'_'+collection+'_'+sample+'_'+str(self.xmin)+'_'+str(self.xmax)+'_'+self.Folder+'.pdf')
                        else:
                            canv.SaveAs(self.outdir+'SignalShape_'+year+'_'+channel+'_'+collection+'_'+sample+'_'+self.Folder+'.pdf')

if __name__ == '__main__':
    PlotBkg = FitSignalShapes()
    PlotBkg.DoFits()
    # for xmin in range(1300,1700+1,100):
    #     for xmax in range(3300,4000+1,100):
    #         if xmin >= xmax: continue
    #         PlotBkg = FitSignalShapes(xmin,xmax)
    #         PlotBkg.DoFits()
