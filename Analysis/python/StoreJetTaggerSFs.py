from Utils import *
from array import array


import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText = "Work in progress"

# ForThesis(TDR)

'''
Module to plot and store Hcc SFs
'''

def StoreJetTaggerSF():

    colors = {"UL16preVFP":   ROOT.kAzure+2,
              "UL16postVFP":  ROOT.kGreen+2,
              "UL17":         ROOT.kRed+1,
              "UL18":         ROOT.kOrange+1,
              "FlavC":        (ROOT.kFullCircle, 1.8, rt.kSolid),
              "FlavB":        (ROOT.kFullSquare, 1.6, rt.kSolid),
              "FlavL":        (ROOT.kFullTriangleUp, 2.0, rt.kSolid),
    }


    SF_Dict = {
        "UL16preVFP": {
            "FlavC": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.161,  1.161,  0.877,  1.215],
                "SF_err_up":   [+0.218*2, +0.218, +0.120, +0.232],
                "SF_err_down": [-0.210*2,-0.210,-0.090,-0.230],
                },
            "FlavB": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.035,  1.035,  1.066,  1.013],
                "SF_err_up":   [+0.101*2, +0.101, +0.095, +0.117],
                "SF_err_down": [-0.094*2, -0.094, -0.083, -0.111],
                },
            "FlavL": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.000,  1.000,  1.000,  1.000],
                "SF_err_up":   [+0.200, +0.200, +0.200, +0.200],
                "SF_err_down": [-0.200, -0.200, -0.200, -0.200],
                },
            },
        "UL16postVFP": {
            "FlavC": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.164,  1.164,  1.078,  1.095],
                "SF_err_up":   [+0.149*2, +0.149, +0.196, +0.228],
                "SF_err_down": [-0.168*2, -0.168, -0.191, -0.269],
                },
            "FlavB": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.033,  1.033,  1.098,  1.061],
                "SF_err_up":   [+0.062*2, +0.062, +0.081, +0.084],
                "SF_err_down": [-0.054*2, -0.054, -0.073, -0.060],
                },
            "FlavL": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.000,  1.000,  1.000,  1.000],
                "SF_err_up":   [+0.200, +0.200, +0.200, +0.200],
                "SF_err_down": [-0.200, -0.200, -0.200, -0.200],
                },
            },
        "UL17": {
            "FlavC": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.316,  1.316,  1.296,  1.115],
                "SF_err_up":   [+0.308*2, +0.308, +0.314, +0.193],
                "SF_err_down": [-0.303*2, -0.303, -0.314, -0.188],
                },
            "FlavB": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [0.957,  0.957,  1.018,  0.976],
                "SF_err_up":   [+0.062*2, +0.062, +0.056, +0.039],
                "SF_err_down": [-0.067*2, -0.067, -0.055, -0.043],
                },
            "FlavL": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.000,  1.000,  1.000,  1.000],
                "SF_err_up":   [+0.200, +0.200, +0.200, +0.200],
                "SF_err_down": [-0.200, -0.200, -0.200, -0.200],
                },
            },
        "UL18": {
            "FlavC": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [0.981, 0.981, 0.946, 0.951],
                "SF_err_up":   [+0.203*2, +0.203, +0.107, +0.143],
                "SF_err_down": [-0.136*2, -0.136, -0.108, -0.122],
                },
            "FlavB": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [0.904,  0.904,  1.005,  0.988],
                "SF_err_up":   [+0.068*2, +0.068, +0.036, +0.030],
                "SF_err_down": [-0.080*2, -0.080, -0.036, -0.038],
                },
            "FlavL": {
                "pt":          [(200, 450),(450, 500), (500, 600), (600, 3000)],
                "SF":          [1.000,  1.000,  1.000,  1.000],
                "SF_err_up":   [+0.200, +0.200, +0.200, +0.200],
                "SF_err_down": [-0.200, -0.200, -0.200, -0.200],
                },
            },
        }

    dir = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/ScaleFactors/Taggers/"
    os.system("mkdir -p "+dir)
    flavs = ["FlavC", "FlavB", "FlavL"]
    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]
    graphs = []
    TDR.cms_lumi = VariablesBase().lumi_map['RunII']['lumiPlot']+' fb^{-1}'
    canvs = {}
    legs = {}
    lines = {}
    for val in [1,1.2,0.8]: 
        lines[val] = rt.TLine(180, val, 1020, val)
    for flavor in flavs:
        canvs[flavor] = tdrCanvas("JetTaggerSF"+flavor, 180, 1020, 0.6, 2.2, "p_{T} [GeV]", "Scale Factor")
        legs[flavor] = tdrLeg(0.70,0.60,0.85,0.85, 0.045, 42, ROOT.kBlack)
        #legs[flavor].SetNColumns(3)
        # tdrHeader(legs[flavor], "Cat. flavour-"+flavor.replace("Flav","").lower().replace("l","light"), textAlign = 22)
        for line in lines.values():
            tdrDrawLine(line, rt.kBlack, rt.kDotted)

    for year in years:
        file = ROOT.TFile(dir+"Tagger_SF_"+year+".root", "RECREATE")

        for flavor in flavs:
            pts = SF_Dict[year][flavor]["pt"]
            SFs = SF_Dict[year][flavor]["SF"]
            SFsErr_Up = SF_Dict[year][flavor]["SF_err_up"]
            SFsErr_Down = SF_Dict[year][flavor]["SF_err_down"]
            pt_bins = sorted(set([pt for bin in pts for pt in bin]))

            name_hist = "SF_Nominal_"+flavor+"_"+year
            name_histUp = "SF_Up_"+flavor+"_"+year
            name_histDown = "SF_Down_"+flavor+"_"+year

            hist     = ROOT.TH1F(name_hist,name_hist,len(pt_bins)-1, array('d',pt_bins))
            histUp   = ROOT.TH1F(name_histUp,name_histUp,len(pt_bins)-1, array('d',pt_bins))
            histDown = ROOT.TH1F(name_histDown,name_histDown,len(pt_bins)-1, array('d',pt_bins))

            x_bins = []
            x_up   = []

            for x in range(hist.GetNbinsX()):
                bin = x+1
                pt = hist.GetBinCenter(bin)
                x_bins.append(hist.GetBinCenter(bin))
                x_up.append(hist.GetBinWidth(bin)/2)
                if pt<pts[x][0] or pt>pts[x][1]: raise ValueError("Unexpected behaviour")
                hist.SetBinContent(bin,SFs[x])
                histUp.SetBinContent(bin,SFs[x]+SFsErr_Up[x])
                histDown.SetBinContent(bin,SFs[x]+SFsErr_Down[x])

            x_bins[0] = 350
            x_up[0] = 150
            x_bins[-1] = 800
            x_up[-1] = 200
            # for x in range(1,hist.GetNbinsX()+1):
            #     print hist.GetBinCenter(x), hist.GetBinContent(x), histUp.GetBinContent(x), histDown.GetBinContent(x)


            hist.Write(name_hist)
            histUp.Write(name_histUp)
            histDown.Write(name_histDown)

            canvs[flavor].cd()
            graph = rt.TGraphAsymmErrors(len(x_bins), array('d',x_bins),array('d',SFs),array('d',x_up),array('d',x_up),array('d',np.abs(SFsErr_Down)),array('d',np.array(SFsErr_Up)))
            graphs.append(graph)
            graph.SetLineWidth(3)
            graph.SetMarkerSize(colors[flavor][1])
            tdrDraw(graph, "P", marker= colors[flavor][0], mcolor=colors[year], lstyle=colors[flavor][2], lcolor= colors[year], fstyle=0, fcolor=colors[year])
            legs[flavor].AddEntry(graph, year, "lp")

        file.Close()

    for flavor in flavs:
        canvs[flavor].SaveAs(dir+"JetTaggerSFs"+flavor+".pdf")


def main():
    StoreJetTaggerSF()

if __name__ == '__main__':
    main()
