from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = 'Simulation'
TDR.extraText2 = 'Work in progress'

f_xmin = 1000
f_xmax = 6000
fNorms = {}

def CrystalBall_Norm(x, par):
    return CrystalBall_Fit(x, par)#*fNorms[str(par[0])+str(par[1])+str(par[2])+str(par[3])]

def ExpGaussExp_Norm(x, par):
    return ExpGaussExp_Fit(x, par)#*fNorms[str(par[0])+str(par[1])+str(par[2])+str(par[3])]


'''
Module to plot Tagger shapes
- Need the Selectionas input
'''

class PlotSignalShapes(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.Samples = ['M'+str(m) for m in self.MassPointsReduced if m!=1000 and m!=1200]
        self.Samples_reduced = ['M1600', 'M2000', 'M2500', 'M3000', 'M3500', 'M4000', 'M4500', 'M5000']
        self.Syst_list = self.Systematics[1:]+self.Systematics_Scale
        self.Syst_list = list(filter(lambda x: not 'PDF' in x and not 'murmuf' in x ,self.Syst_list))
        self.Syst = ['']+[x+v for x in self.Syst_list for v in ['Up', 'Down']]
        # self.Channels = ['muon', 'electron']
        self.Folder = 'PN_ZHccvsQCD_MD'

        self.color  = {'M1600': ROOT.kAzure+2,
                       'M2000': ROOT.kAzure+7,
                       'M2500': ROOT.kGreen+3,
                       'M3000': ROOT.kGreen+2,
                       'M3500': ROOT.kOrange+1,
                       'M4000': ROOT.kOrange-2,
                       'M4500': ROOT.kRed+2,
                       'M5000': ROOT.kRed+1,
                       }
        self.style  = {'':           ROOT.kSolid,
                       'up':         ROOT.kDashed,
                       'down':       ROOT.kDotted,
                       'Nominal':    ROOT.kSolid,
                       'Syst. up':   ROOT.kDashed,
                       'Syst. down': ROOT.kDotted,
                       }
        self.outdir = self.Path_ANALYSIS+'Analysis/OtherPlots/SignalShapes/'
        os.system('mkdir -p '+self.outdir)

        self.LoadFuncions()

    def LoadFuncions(self):
        self.funcs = {}
        for year in ['RunII']:
            for channel in self.Channels:
                for collection in self.Collections:
                    for sample in self.Samples:
                        for sys in self.Syst:
                            if DoControl([''], year+channel+sys, channel, sample): continue
                            name = year+channel+collection+sample+sys
                            self.funcs[name] = ROOT.TF1(name, ExpGaussExp_Norm if 'inv' in channel else CrystalBall_Norm, f_xmin, f_xmax, 5)
                    fname = self.Path_ANALYSIS+'/Analysis/Limits/nominal/'+year+'/'+collection+'/'+channel+'channel/'+self.Folder+'/datacards/SignalProperties_'+year+'_'+self.Folder+'.txt'
                    with open(fname, 'r') as f_:
                        for line in f_.readlines():
                            # if not any(x in line for x in self.Syst_list): continue
                            if 'Up' in line and not any(x in line for x in self.Syst_list): continue
                            if 'Down' in line and not any(x in line for x in self.Syst_list): continue
                            # if 'Down' in line and not 'btag' in line: continue
                            if 'corrected' in line:
                                sys = line.split()[0].split('0')[-1]
                                sample = line.split()[0].replace(sys,'')
                                if not any(x in sample for x in self.Samples): continue
                                if not any(x in sys for x in self.Syst): continue
                                if DoControl([''], year+channel+sys, channel, sample): continue
                                self.funcs[year+channel+collection+sample+sys].SetParameter(4,float(line.split()[-1]))
                                self.funcs[year+channel+collection+sample+'stat4'] = 2
                            if not 'param' in line: continue
                            par = line.split()[0].replace('_'+channel+'_'+year,'')
                            var = float(line.split()[2])
                            sys = par.split('0')[-1]
                            sample = par.split('_')[2].replace(sys,'')
                            if not any(x in sample for x in self.Samples): continue
                            if not any(x in sys for x in self.Syst): continue
                            if DoControl([''], year+channel+sys, channel, sample): continue
                            self.funcs[year+channel+collection+sample+sys].SetParameter(int(par.split('_')[1][1]),var)
                            if sys == '':
                                self.funcs[year+channel+collection+sample+'stat'+par.split('_')[1][1]] = float(line.split()[3])
                    for sample in self.Samples:
                        for sys in self.Syst:
                            if DoControl([''], year+channel+sys, channel, sample): continue
                            name = year+channel+collection+sample+sys
                            # norm = self.funcs[year+channel+collection+sample+sys].GetParameter(4)
                            # self.funcs[year+channel+collection+sample+sys].SetParameter(4,1)
                            # self.funcs[year+channel+collection+sample+sys].SetParameter(4,norm/1000)
                            # func = self.funcs[name]
                            # fname = str(func.GetParameter(0))+str(func.GetParameter(1))+str(func.GetParameter(2))+str(func.GetParameter(3))
                            # f_ = ROOT.TF1('', ExpGaussExp_Fit if 'inv' in channel else CrystalBall_Fit, f_xmin,f_xmax, 5)
                            # f_.SetParameters(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2),func.GetParameter(3), 1)
                            # fNorms[fname] = 1./f_.Integral(0,10000)

    def PlotFuncions(self):
        for channel in self.Channels:
            TDR.cms_lumi = self.lumi_map['RunII']['lumiPlot']+' fb^{-1}'
            isInv = 'inv' in channel
            canv = tdrCanvas(channel, f_xmin, f_xmax, 0.0001, 30 if not isInv else 100, "M(Z') [GeV]" if not isInv else "M_{T}(Z') [GeV]" ,'A.U.', isExtraSpace=True)
            hframe = ROOT.gROOT.FindObject('hframe')
            hframe.GetXaxis().SetNdivisions(5)
            leg2 = tdrLeg(0.42,0.66,0.65,0.89, 0.035, 42, ROOT.kBlack)
            leg = tdrLeg(0.60,0.60,0.92,0.89, 0.035, 42, ROOT.kBlack)
            leg.SetNColumns(2)
            fs = []
            for sys in ['Nominal','Syst. up', 'Syst. down']:
                f_ = ROOT.TF1()
                f_.SetLineColor(ROOT.kBlack)
                f_.SetLineStyle(self.style[sys])
                fs.append(f_)
                leg2.AddEntry(f_, sys,'l')
            for year in ['RunII']:
                for collection in self.Collections:
                    for sample in self.Samples_reduced:
                        for sys in self.Syst:
                            if DoControl([''], year+channel+sys, channel, sample): continue
                            func = self.funcs[year+channel+collection+sample+sys]
                            func.SetLineColor(self.color[sample])
                            func.SetLineStyle(self.style['down' if 'Down' in sys else ('up' if 'Up' in sys else '')])
                            func.SetNpx(1000)
                            func.Draw('same')
                            if sys=='': leg.AddEntry(func, "Z' "+str(round(float(sample.replace('M',''))/1000.,1))+' TeV', 'l')
            canv.SaveAs(self.outdir+'SignalShape_'+channel+'.pdf')

    def PlotParameters(self):
        for par in range(0,6):
            if par == 0: ymin, ymax, ParName, funcName = (0,6000, '#mu [GeV]',    'pol1')
            if par == 1: ymin, ymax, ParName, funcName = (0,350,  '#sigma [GeV]', 'pol2')
            if par == 2: ymin, ymax, ParName, funcName = (0,2,    '#alpha',       'pol2')
            if par == 3: ymin, ymax, ParName, funcName = (0,3,    'n',            'pol2')
            if par == 4: ymin, ymax, ParName, funcName = (0,30,   'Norm',         'pol2')
            if par == 5: ymin, ymax, ParName, funcName = (0,2,    'lnN',          'pol0')
            TDR.cms_lumi = self.lumi_map['RunII']['lumiPlot']+' fb^{-1}'
            yName = 'parameter '+str(par+1)
            if par == 4: yName = 'Normalisation'
            if par == 5: yName = 'lnN'
            canv = tdrCanvas(str(par), f_xmin+200, f_xmax-800, ymin, ymax, "M(Z') [GeV]" ,yName, isExtraSpace=True)
            leg = tdrLeg(0.42,0.66,0.65,0.89, 0.035, 42, ROOT.kBlack)
            LegendNames = {}
            graphs = []
            for channel in self.Channels:
                isInv = 'inv' in channel
                isMuo = 'muon' in channel
                color = ROOT.kAzure+2 if isInv else (ROOT.kRed+1 if isMuo else ROOT.kGreen+2)
                markerStyle = ROOT.kFullCircle if isInv else (ROOT.kFullTriangleUp if isMuo else ROOT.kFullTriangleDown)
                if par==2:
                    ParName = 'k_{L}' if isInv else '#alpha'
                if par==3:
                    ParName = 'k_{H}' if isInv else 'n'
                for year in ['RunII']:
                    for collection in self.Collections:
                        x_bins = []
                        y_bins = []
                        x_up = []
                        y_stat = []
                        y_up = []
                        y_down = []
                        for sample in self.Samples:
                            x_bins.append(float(sample.replace('M','')))
                            nominal = self.funcs[year+channel+collection+sample].GetParameter(par) if par!=5 else 1
                            up = 0
                            down = 0
                            for sys in self.Syst:
                                if DoControl([''], year+channel+sys, channel, sample): continue
                                if par!=5:
                                    toadd = nominal-self.funcs[year+channel+collection+sample+sys].GetParameter(par)
                                    if sys =='':
                                        stat = self.funcs[year+channel+collection+sample+'stat'+str(par)]
                                    else:
                                        if toadd>0:
                                            up += toadd**2
                                        else:
                                            down += toadd**2
                                else:
                                    stat = 0
                                    up   += (1-lnN(self.funcs[year+channel+collection+sample].GetParameter(4),self.funcs[year+channel+collection+sample+sys].GetParameter(4)))**2
                                    down += (1-lnN(self.funcs[year+channel+collection+sample].GetParameter(4),self.funcs[year+channel+collection+sample+sys].GetParameter(4)))**2

                            up = np.sqrt(up+ stat**2)
                            down = np.sqrt(down+ stat**2)
                            x_up.append(0.)
                            if isInv and par ==4:
                                nominal /= 3.
                                stat    /= 3.
                                up      /= 3.
                                down    /= 3.
                            if isInv and par ==2:
                                nominal *= 3.
                                stat    *= 3.
                                up      *= 3.
                                down    *= 3.
                            y_bins.append(nominal)
                            y_stat.append(stat)
                            y_up.append(up)
                            y_down.append(down)
                        graph_stat = ROOT.TGraphAsymmErrors(len(x_bins), array('d',x_bins),array('d',y_bins),array('d',x_up),array('d',x_up),array('d',y_stat),array('d',y_stat))
                        graph_syst = ROOT.TGraphAsymmErrors(len(x_bins), array('d',x_bins),array('d',y_bins),array('d',x_up),array('d',x_up),array('d',y_down),array('d',y_up))
                        func_stat = rt.TF1('par'+str(par)+year+channel+collection+sample+'stat',funcName,x_bins[0]-100, x_bins[-1]+100)
                        func_syst = rt.TF1('par'+str(par)+year+channel+collection+sample+'syst',funcName,x_bins[0]-100, x_bins[-1]+100)
                        graphs.append(func_stat)
                        graphs.append(func_syst)
                        style_stat = ROOT.kSolid
                        style_syst = ROOT.kDashed
                        func_stat.SetLineColor(color)
                        func_syst.SetLineColor(color)
                        func_stat.SetLineStyle(style_stat)
                        func_syst.SetLineStyle(style_syst)
                        rt.gStyle.SetOptFit(0)
                        fitRes_stat = graph_stat.Fit(func_stat,'RMQS')
                        fitRes_syst = graph_syst.Fit(func_syst,'RMQS0')
                        band_stat, _, _ = ComputeHistWithCL('par'+str(par)+year+channel+collection+sample+'_stat', func_stat, fitRes_stat, graph_stat.GetHistogram(), cl=0.95)
                        band_syst, _, _ = ComputeHistWithCL('par'+str(par)+year+channel+collection+sample+'_syst', func_syst, fitRes_syst, graph_syst.GetHistogram(), cl=0.95)
                        tdrDraw(band_stat, 'e4', 9, ROOT.kWhite, 1, color+1, 1001, color+1)
                        band_stat.SetMarkerSize(0)
                        band_syst.SetMarkerSize(0)
                        band_stat.SetFillColorAlpha(color+1,0.35)
                        band_syst.SetFillColorAlpha(color+1,0.35)
                        tdrDraw(graph_syst, 'pe', markerStyle, color, style_syst, color, 1001, color)
                        tdrDraw(graph_stat, 'pe', markerStyle, color, style_stat, color, 1001, color)
                        if channel not in LegendNames:
                            leg.AddEntry(graph_syst, channel+' -- '+ParName+(' *0.33' if isInv and par ==4 else (' *3' if isInv and par ==2 else '')), 'lp')
                        graphs.append(graph_stat)
                        graphs.append(graph_syst)
                        graphs.append(band_stat)
                        graphs.append(band_syst)
            canv.SaveAs(self.outdir+'SignalParameter_'+str(par)+'_withFit.pdf')



if __name__ == '__main__':
    Plotter = PlotSignalShapes()
    Plotter.PlotFuncions()
    Plotter.PlotParameters()
