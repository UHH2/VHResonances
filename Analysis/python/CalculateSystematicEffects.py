from Utils import *

from math import sqrt

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

ForThesis(TDR)

'''
Module to plot LeptonID Efficiency

- Need the LeptonIDStudiesModule output as input

'''


class CalculateSystematicEffects(VariablesBase):
    def __init__(self, year = "RunII", channel = "muonchannel", histFolder="btag_DeepBoosted_H4qvsQCD", studies = "nominal", collection="Puppi"):
        VariablesBase.__init__(self)
        TDR.cms_lumi = self.lumi_map[year]['lumiPlot']+' fb^{-1}'
        self.year = year
        self.channel = channel.replace("channel","")
        self.isInv   = "inv" in self.channel
        self.histFolder = histFolder
        self.studies = studies
        self.collection = collection
        self.inputdir = self.Path_ANALYSIS+"Analysis/Limits/"+self.studies+"/"+self.year+"/"+self.collection+"/"+self.channel+"channel/"+self.histFolder+"/datacards/"
        self.InputFilename  = "SignalProperties_"+self.year+"_"+self.histFolder+".txt"
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/Systematics/"
        os.system("mkdir -p "+self.outdir)

    def ReadFile(self):
        self.parametesMap = {}
        # print self.inputdir+self.InputFilename
        with open(self.inputdir+self.InputFilename, "r") as f:
            for line in f.readlines():
                if "param" in line:
                    self.parametesMap[line.split()[0]] = (float(line.split()[2]),float(line.split()[3]))
                    # print line.split()[0], self.parametesMap[line.split()[0]]
                if "events" in line:
                    name = "norm"+line.split()[0]+"_"+self.channel+"_"+self.year
                    if not name in self.parametesMap: self.parametesMap[name] = (0,0)
                    if "corrected" in line: self.parametesMap[name] = (self.parametesMap[name][0], float(line.split()[-1]))
                    else: self.parametesMap[name] = (float(line.split()[-1]), self.parametesMap[name][1])
                    # print name, self.parametesMap[name]
        # prettydic(self.parametesMap)

    def DoPlots(self):
        gr_all = {}
        dummy = array('d',np.array([0]*len(self.MassPointsReduced)))
        for par in ["norm", "sg_p0_", "sg_p1_", "sg_p2_", "sg_p3_"]:
            y_min = 2*1e-03
            y_max = (50 if self.isInv else 15)*float(self.lumi_map[year]["lumi_fb"])/float(self.lumi_map["RunII"]["lumi_fb"])
            if par=="sg_p0_": y_max = 9000
            if par=="sg_p1_": y_max = (300 if self.isInv else 250)
            if par=="sg_p2_": y_max = 1.5
            if par=="sg_p3_": y_max = (200 if self.isInv else 10)

            canv = tdrCanvas(par, 1200, 5200, y_min, y_max, "M(Zprime)", par.replace("_",""))
            values = []
            stat = []
            syst = []
            for mass in self.MassPointsReduced:
                nominal = self.parametesMap[par+"M"+str(mass)+"_"+self.channel+"_"+self.year][0]
                nominal_var = self.parametesMap[par+"M"+str(mass)+"_"+self.channel+"_"+self.year][1]
                if par=="norm": nominal_var -= nominal
                max_ = -1
                for sys in self.Systematics+ self.Systematics_Scale:
                    if sys=="nominal": continue
                    # if "murmuf" in sys: continue
                    for var in ["Up","Down"]:
                        if DoControl([""], collection+channel+sys+self.Signal+("_inv" if self.isInv else "")+"_M"+str(mass), channel, self.Signal+("_inv" if self.isInv else "")+"_M"+str(mass)): continue
                        val  = abs(nominal-self.parametesMap[par+"M"+str(mass)+sys+var+"_"+self.channel+"_"+self.year][0])
                        max_ = max(max_,sqrt(val**2+nominal_var**2))
                values.append(nominal)
                stat.append(nominal_var)
                syst.append(max_)
            values = array('d',np.array(values))
            stat = array('d',np.array(stat))
            syst = array('d',np.array(syst))
            gr_syst = rt.TGraphErrors(len(self.MassPointsReduced), array('d',self.MassPointsReduced), values, dummy, syst)
            gr_syst1 = rt.TGraphErrors(len(self.MassPointsReduced), array('d',self.MassPointsReduced), values, dummy, syst)
            gr = rt.TGraphErrors(len(self.MassPointsReduced), array('d',self.MassPointsReduced), values, dummy, stat)
            gr_all[par+"syst"] = gr_syst
            gr_all[par] = gr
            func_name = "pol2"
            if "norm" in par and "ele" in channel: func_name = "pol4"
            if "p0" in par: func_name = "pol1"
            if "p1" in par: func_name = "pol2"
            if "p2" in par: func_name = "pol0"
            if "p3" in par: func_name = "pol1"
            func = rt.TF1("func"+par,func_name,self.MassPointsReduced[0], self.MassPointsReduced[-1])
            gr_all[par+"func"] = func
            func.SetLineColor(rt.kBlack)
            rt.gStyle.SetOptFit(0)
            gr_syst1.Fit(func,"RQ")
            rt.TVirtualFitter.GetFitter().GetConfidenceIntervals(gr_syst1)
            tdrDraw(gr_syst1, "l3", rt.kDot, rt.kOrange-2, rt.kSolid, rt.kOrange-2, 1001,rt.kOrange-2)
            func.Draw("same")
            tdrDraw(gr_syst, "lp", rt.kFullCircle, rt.kRed, rt.kSolid, rt.kRed, 1001, rt.kRed)
            tdrDraw(gr, "lp", rt.kFullCircle, rt.kBlack, rt.kSolid, rt.kBlack, 1001, rt.kBlack)
            canv.SaveAs(self.outdir+"variations_"+par+year+channel+histFolder+".pdf")

    def DoCalculations(self):
        SystVariations = {}
        for par in ["norm", "sg_p0_", "sg_p1_", "sg_p2_", "sg_p3_"]:
            for mass in self.MassPointsReduced:
                nominal = self.parametesMap[par+"M"+str(mass)+"_"+self.channel+"_"+self.year][0]
                nominal_var = self.parametesMap[par+"M"+str(mass)+"_"+self.channel+"_"+self.year][1] if par!="norm" else nominal
                variations = {}
                for syst in self.Systematics+self.Systematics_Scale:
                    if syst == "nominal": continue
                    for var in ["Up", "Down"]:
                        systvar = self.parametesMap[par+"M"+str(mass)+syst+var+"_"+self.channel+"_"+self.year]
                        variations[syst+var] = systvar[0]
                        SystVariations[str(mass)+par+syst+var] = (round(abs(nominal-systvar[0])/(nominal_var+systvar[1]),3),round(abs(1.-systvar[0]/nominal)*100,3))
                        # print str(mass)+par+syst+var, SystVariations[str(mass)+par+syst+var], nominal, systvar[0], nominal_var, systvar[1]
                        # SystVariations.setdefault(syst,[]).append(round(abs(nominal-variations[syst+var])/nominal_var,3))
                # PrintFormattedLine([mass,par,nominal,nominal_var,round(np.mean(variations.values()),2),round(np.std(variations.values()),2),round(np.max(abs(nominal-np.array(variations.values()))/nominal_var),3)])
        for x in SystVariations:
            if SystVariations[x][0]>1.9 and SystVariations[x][1]>5:
                print x, SystVariations[x]
        for mass in self.MassPointsReduced:
            for par in ["norm", "sg_p0_", "sg_p1_", "sg_p2_", "sg_p3_"]:
                list_ = [SystVariations[x] for x in SystVariations if str(mass) in x and par in x]
                print mass, par, np.max(list_), np.mean(list_)
        for syst in self.Systematics+self.Systematics_Scale:
            if syst == "nominal": continue
            for par in ["norm", "sg_p0_", "sg_p1_", "sg_p2_", "sg_p3_"]:
                list_ = [SystVariations[x] for x in SystVariations if syst in x and par in x]
                print syst, par, np.max(list_), np.mean(list_)





if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years if len(args.years)!=0 else ["RunII"]
    histFolders = args.histFolders if len(args.histFolders)!=0 else ["DeepAk8_ZHccvsQCD_MD2"]
    Channels    = args.Channels if len(args.Channels)!=0 else ["invisiblechannel"]
    Collections = args.Collections if len(args.Collections)!=0 else ["Puppi"]

    studies = "nominal"

    for year in years:
        for histFolder in histFolders:
            for channel in Channels:
                for collection in Collections:
                    CSE = CalculateSystematicEffects(year=year, studies=studies, histFolder=histFolder, channel=channel, collection=collection)
                    CSE.ReadFile()
                    CSE.DoPlots()
                    # CSE.DoCalculations()
