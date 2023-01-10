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

    colors = {"2016":  ROOT.kGreen+2,
              "2017":  ROOT.kRed+1,
              "2018":  ROOT.kOrange+1,
              "FlavC": (ROOT.kFullCircle, 1.8, rt.kSolid),
              "FlavB": (ROOT.kFullTriangleUp, 2.0, rt.kDashed),
              "FlavL": (ROOT.kFullSquare, 1.6, rt.kDotted),
    }


    SF_Dict = {
        "2016": {
            "FlavC": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.069,+1.122,+1.121,+1.074,+0.911,+0.998],
                "SF_err_up":   [+0.092,+0.086,+0.085,+0.067,+0.046,+0.049],
                "SF_err_down": [-0.083,-0.077,-0.075,-0.061,-0.042,-0.046],
                },
            "FlavB": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.083,+1.070,+1.141,+1.056,+1.186,+1.082],
                "SF_err_up":   [+0.087,+0.085,+0.074,+0.067,+0.061,+0.070],
                "SF_err_down": [-0.086,-0.086,-0.074,-0.068,-0.062,-0.071],
                },
            "FlavL": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.136,+1.029,+1.057,+1.165,+1.265,+1.265],
                "SF_err_up":   [+0.100,+0.078,+0.068,+0.074,+0.054,+0.060],
                "SF_err_down": [-0.096,-0.075,-0.067,-0.072,-0.053,-0.059],
                },
            },
        "2017": {
            "FlavC": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.116,+1.091,+1.110,+0.964,+1.046,+1.050],
                "SF_err_up":   [+0.086,+0.076,+0.084,+0.077,+0.064,+0.049],
                "SF_err_down": [-0.078,-0.069,-0.075,-0.070,-0.059,-0.046],
                },
            "FlavB": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.241,+1.308,+1.326,+1.079,+1.181,+1.009],
                "SF_err_up":   [+0.112,+0.100,+0.093,+0.102,+0.073,+0.060],
                "SF_err_down": [-0.111,-0.099,-0.093,-0.105,-0.071,-0.061],
                },
            "FlavL": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.346,+1.417,+1.363,+1.744,+1.476,+1.693],
                "SF_err_up":   [+0.101,+0.084,+0.100,+0.117,+0.080,+0.065],
                "SF_err_down": [-0.098,-0.083,-0.097,-0.108,-0.077,-0.063],
                },
            },
        "2018": {
            "FlavC": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.019,+0.889,+1.060,+1.096,+1.087,+1.072],
                "SF_err_up":   [+0.087,+0.068,+0.091,+0.103,+0.080,+0.048],
                "SF_err_down": [-0.074,-0.058,-0.081,-0.091,-0.072,-0.045],
                },
            "FlavB": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.282,+1.464,+1.071,+1.116,+1.062,+1.172],
                "SF_err_up":   [+0.105,+0.098,+0.105,+0.118,+0.078,+0.057],
                "SF_err_down": [-0.102,-0.102,-0.105,-0.118,-0.078,-0.058],
                },
            "FlavL": {
                "pt":          [(200, 250),(250, 300),(300, 350),(350, 400),(400, 500),(500, 2000)],
                "SF":          [+1.152,+1.027,+1.245,+0.963,+1.109,+1.102],
                "SF_err_up":   [+0.099,+0.085,+0.111,+0.106,+0.074,+0.057],
                "SF_err_down": [-0.100,-0.091,-0.106,-0.102,-0.071,-0.055],
                },
            },
        }

    dir = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/ScaleFactors/Taggers/"
    os.system("mkdir -p "+dir)

    graphs = []
    TDR.cms_lumi = self.lumi_map['RunII']['lumiPlot']+' fb^{-1}'
    canvs = {}
    legs = {}
    for flavor in ["FlavC", "FlavB", "FlavL"]:
        canvs[flavor] = tdrCanvas("JetTaggerSF"+flavor, 180, 620, 0.6, 2.6, "p_{T} [GeV]", "Scale Factor")
        legs[flavor] = tdrLeg(0.70,0.60,0.85,0.85, 0.045, 42, ROOT.kBlack)
        #legs[flavor].SetNColumns(3)
        tdrHeader(legs[flavor], "Cat. flavour-"+flavor.replace("Flav","").lower().replace("l","light"), textAlign = 22)


    for year in ["2016", "2017", "2018"]:
        file = ROOT.TFile(dir+"Tagger_SF_"+year+".root", "RECREATE")

        for flavor in ["FlavC", "FlavB", "FlavL"]:
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

            x_bins.pop()
            x_up.pop()
            x_bins.append(550)
            x_up.append(50)
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
            tdrDraw(graph, "P", colors[flavor][0], colors[year], colors[flavor][2], colors[year], 0, colors[year])
            legs[flavor].AddEntry(graph, year, "lp")

        file.Close()

    for flavor in ["FlavC", "FlavB", "FlavL"]:
        canvs[flavor].SaveAs(dir+"JetTaggerSFs"+flavor+".pdf")


def main():
    StoreJetTaggerSF()

if __name__ == '__main__':
    main()
