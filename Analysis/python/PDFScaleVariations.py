from Utils import *
from array import array
import json

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

ForThesis(TDR)

class PDFScaleVariations(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        TDR.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/PDFScaleVariations/"
        os.system("mkdir -p "+self.outdir)
        self.histfolder = "ZprimeCandidate_DeepAk8_ZHccvsQCD_MD_SR"
        self.hname = "Zprime_mass_rebin100"
        self.years = self.years+["RunII"]
        self.years = ["RunII"]
        with open(self.Path_ANALYSIS+"Analysis/OtherPlots/SignalNormalization/normalization.json", 'r') as f_:
            self.NormInfo = json.load(f_)


    def PlotWeights(self):
        self.histos = {}
        weight_lenght = 110
        canv  = tdrCanvas("weights", -1, weight_lenght+1, 0.9, 1.1, "index", "Variation") # The range is quite variable depending on the year
        norm = 0
        shift = 9
        hists = []
        year = "2016"
        ch   = "muon"
        mass = 1400
        fname = self.Path_STORAGE+year+"/Selection/Puppi/"+ch+"channel/nominal/workdir_Selection_MC_ZprimeToZH_M"+str(mass)+"_"+year+"/"+self.PrefixrootFile+"MC."+self.Signal+"_M"+str(mass)+"_"+year+"_0.root"
        if not os.path.isfile(fname):
            print "NOT FOUND,", fname
            return
        f_ = rt.TFile(fname)
        t_ = f_.Get("AnalysisTree")
        hist = rt.TH1D("weights", "weights", weight_lenght+1, 0, weight_lenght)
        hist_wrtNominal = rt.TH1D("weights wrt nominal", "weights", weight_lenght, 0, weight_lenght)
        for ev in t_:
            for ind in range(weight_lenght):
                fill = ev.genInfo.systweights().at(ind)
                ref = ev.genInfo.systweights().at(0) if ind<shift else ev.genInfo.systweights().at((int(ind-shift)/100)*100 +shift)
                hist.Fill(ind+1, fill)
                hist_wrtNominal.Fill(ind+1, fill/ref)
                hist.SetBinError(ind+1,0)
                hist_wrtNominal.SetBinError(ind+1,0)
            break
        hist.SetDirectory(0)
        hist_wrtNominal.SetDirectory(0)
        hists.append(hist)
        hists.append(hist_wrtNominal)
        tdrDraw(hist, "p", rt.kFullCircle, rt.kBlack, rt.kSolid, rt.kBlack, 0)
        tdrDraw(hist_wrtNominal, "p", rt.kFullCircle, rt.kBlue+1, rt.kSolid, rt.kBlue+1, 0)
        canv.SaveAs(self.outdir+"Weights.pdf")




    def PlotVariations(self):
        self.Norm = {}
        for year in self.years:
            self.Norm[year] = {}
            for ch in self.Channels:
                hname = self.hname.replace("mass", "mass_transversal") if "inv" in ch else self.hname
                self.Norm[year][ch] = {}
                for mass in self.MassPointsReduced:
                # for mass in [5000]:
                    self.Norm[year][ch][mass] = {"Up": {}, "Down": {}, "Nominal": 0}
                    sample = self.Signal+("_inv" if "inv" in ch else "")+"_M"+str(mass)+"_"+year
                    fname = self.Path_STORAGE+year+"/SignalRegion/Puppi/"+ch+"channel/nominal/"+self.PrefixrootFile+"MC."+sample+"_"+"noTree.root"
                    if not os.path.isfile(fname):
                        print "NOT FOUND,", fname
                        continue
                    f_ = rt.TFile(fname)
                    canv  = tdrCanvas(year+ch+str(mass), mass*0.5, mass*1.5, 1, 1e4, "M(Zprime)", "Events")
                    canv.SetLogy(1)
                    leg = tdrLeg(0.60,0.70,0.95,0.89, 0.025, 42, ROOT.kBlack)

                    nominal = f_.Get(self.histfolder+"/"+hname)
                    self.Norm[year][ch][mass]["Nominal"] = nominal.Integral()
                    nominal.SetLineWidth(2)
                    tdrDraw(nominal, "hist", rt.kFullCircle, rt.kGreen+2, rt.kSolid, rt.kGreen+2, 0)
                    leg.AddEntry(nominal, "nominal", "lp")

                    ScaleVarTotUp = nominal.Clone("ScaleVariationUp")
                    ScaleVarTotDown = nominal.Clone("ScaleVariationDown")
                    hists = []
                    hists.append(nominal)
                    for var in ["upup", "upnone", "noneup", "nonedown", "downnone", "downdown"]:
                        h = f_.Get(self.histfolder.replace("ZprimeCandidate", "ZprimeCandidate_murmuf_"+var)+"/"+hname)
                        # h.Scale(nominal.Integral()/h.Integral()) TODO??????
                        hists.append(h)
                    for bin in range(0,nominal.GetNbinsX()+1):
                        y_vals = []
                        for h in hists:
                            y_vals.append(h.GetBinContent(bin))
                        ScaleVarTotUp.SetBinContent(bin, np.max(y_vals))
                        ScaleVarTotDown.SetBinContent(bin, np.min(y_vals))
                    self.Norm[year][ch][mass]["Up"]["QCD Scale"] = ScaleVarTotUp.Integral()/self.Norm[year][ch][mass]["Nominal"]
                    self.Norm[year][ch][mass]["Down"]["QCD Scale"] = ScaleVarTotDown.Integral()/self.Norm[year][ch][mass]["Nominal"]
                    ScaleVarTotUp.SetLineWidth(2)
                    ScaleVarTotDown.SetLineWidth(2)
                    tdrDraw(ScaleVarTotUp, "hist", rt.kFullCircle, rt.kRed+1, rt.kSolid, rt.kRed+1, 0)
                    tdrDraw(ScaleVarTotDown, "hist", rt.kFullCircle, rt.kOrange+1, rt.kSolid, rt.kOrange+1, 0)
                    leg.AddEntry(ScaleVarTotUp, "Scale Up", "lp")
                    leg.AddEntry(ScaleVarTotDown, "Scale Down", "lp")


                    all_hists = []
                    for PDF in ["PDF", "NNPDF", "NNPDF31_lo_as_0130", "PDF4LHC15_nnlo_100"]:
                    # for PDF in ["PDF"]:
                        PDFVarTotUp = nominal.Clone("PDFVariationUp")
                        PDFVarTotDown = nominal.Clone("PDFVariationDown")
                        hists = []
                        for var in range(0,100):
                            h = f_.Get(self.histfolder.replace("ZprimeCandidate", "ZprimeCandidate_"+PDF+"_"+str(var))+"/"+hname)
                            if not h: continue
                            # ratio = nominal.Integral()/h.Integral()
                            ratio = (nominal.Integral()/h.Integral())/(self.NormInfo[year]["Puppi"][ch][sample]["weights"]/self.NormInfo[year]["Puppi"][ch][sample][PDF+"_"+str(var)+"_weights"])
                            # if abs(ratio-1)>0.01: print ratio
                            # h.Scale(ratio)
                            hists.append(h)
                        for bin in range(0,nominal.GetNbinsX()+1):
                            y_vals = []
                            for h in hists:
                                y_vals.append(h.GetBinContent(bin))
                            y_vals = np.array(y_vals)
                            var = np.std(y_vals)
                            # var = np.std(y_vals[abs(y_vals/nominal.GetBinContent(bin)-1) < 1])
                            # if np.isnan(var) : var = 0
                            PDFVarTotUp.SetBinContent(bin,   nominal.GetBinContent(bin) + var)
                            PDFVarTotDown.SetBinContent(bin, nominal.GetBinContent(bin) - var)
                        self.Norm[year][ch][mass]["Up"][PDF] = PDFVarTotUp.Integral()/self.Norm[year][ch][mass]["Nominal"]
                        self.Norm[year][ch][mass]["Down"][PDF] = PDFVarTotDown.Integral()/self.Norm[year][ch][mass]["Nominal"]
                        PDFVarTotUp.SetLineWidth(2)
                        PDFVarTotDown.SetLineWidth(2)
                        tdrDraw(PDFVarTotUp, "hist", rt.kFullCircle, rt.kBlue+1, rt.kSolid, rt.kBlue+1, 0)
                        tdrDraw(PDFVarTotDown, "hist", rt.kFullCircle, rt.kAzure+7, rt.kSolid, rt.kAzure+7, 0)
                        leg.AddEntry(PDFVarTotUp, PDF+" Up", "lp")
                        leg.AddEntry(PDFVarTotDown, PDF+" Down", "lp")
                        all_hists.append(PDFVarTotUp)
                        all_hists.append(PDFVarTotDown)

                    canv.SaveAs(self.outdir+year+ch+"_M"+str(mass)+".pdf")
                    f_.Close()
                    graph = rt.TGraph()


        canv_Var  = tdrCanvas("Variatios", self.MassPointsReduced[0]-50, self.MassPointsReduced[-1]+100, 0, 2, "M(Zprime)", "Variations")
        leg_Var = tdrLeg(0.60,0.70,0.95,0.89, 0.025, 42, ROOT.kBlack)
        graphs = []
        colors = {"2016":       ROOT.kGreen+2,
                  "2017":       ROOT.kRed+1,
                  "2018":       ROOT.kOrange+1,
                  "RunII":      ROOT.kBlue+1,
                  "lepton":     ROOT.kFullCircle,
                  "muon":       ROOT.kFullTriangleDown,
                  "electron":   ROOT.kFullTriangleUp,
                  "invisible":  ROOT.kFullSquare,
                  "PDF" :                ROOT.kGreen+2, # NNPDF31_nnlo_hessian_pdfas
                  "NNPDF" :              ROOT.kBlue+1, # NNPDF31_nnlo_as_0118_nf_4_mc_hessian, reweighted for invisible
                  "NNPDF31_lo_as_0130" : ROOT.kRed+1, # Largest scale smaller shape
                  "PDF4LHC15_nnlo_100" : ROOT.kOrange+1,
        }
        # for year in self.years:
        for year in ["RunII"]:
            for ch in self.Channels:
                # for PDF in ["PDF", "NNPDF", "NNPDF31_lo_as_0130", "PDF4LHC15_nnlo_100"]:
                # for PDF in ["PDF", "NNPDF", ]:
                # for PDF in ["PDF", "NNPDF31_lo_as_0130"]:
                for PDF in ["PDF", "PDF4LHC15_nnlo_100"]:
                    y_ = []
                    y_errUp = []
                    y_errDown = []
                    for mass in self.MassPointsReduced:
                        y_errUp.append(self.Norm[year][ch][mass]["Up"][PDF])
                        y_errDown.append(self.Norm[year][ch][mass]["Down"][PDF])
                    gr_Up = rt.TGraph(len(self.MassPointsReduced), array('d',self.MassPointsReduced), array('d',y_errUp))
                    gr_Down = rt.TGraph(len(self.MassPointsReduced), array('d',self.MassPointsReduced), array('d',y_errDown))
                    graphs.append(gr_Up)
                    graphs.append(gr_Down)
                    if "inv" in ch:
                        gr_Up.SetMarkerSize(1.5)
                        gr_Down.SetMarkerSize(1.5)
                    else:
                        gr_Up.SetMarkerSize(2)
                        gr_Down.SetMarkerSize(2)
                    gr_Up.SetLineWidth(2)
                    gr_Down.SetLineWidth(2)
                    tdrDraw(gr_Up, "PL", colors[ch], colors[PDF], ROOT.kSolid, colors[PDF], 1000, colors[PDF] )
                    tdrDraw(gr_Down, "PL", colors[ch], colors[PDF], ROOT.kDashed, colors[PDF], 1000, colors[PDF] )
                    if "inv" in ch:
                        name = PDF.replace("_as_0130","").replace("_100","")
                        if "PDF"==PDF[-3:]:
                            name = name.replace("PDF","PDF (NNPDF31_nnlo)")
                        leg_Var.AddEntry(gr_Up, name, "lp")

        canv_Var.SaveAs(self.outdir+"Variations.pdf")



def main():
    Var = PDFScaleVariations()
    # Var.PlotWeights()
    Var.PlotVariations()

if __name__ == '__main__':
    main()
