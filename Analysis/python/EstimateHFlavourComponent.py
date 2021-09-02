from Utils import *

'''
Module to extract Tagger info. The output will be used later on to make quick checks
Need to choose if you want the whole dataset or not via the isFast option.
'''


class EstimateHFlavorComponent(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/"
        os.system("mkdir -p "+self.outdir)

    def GetJetFlavor(self,event, jet, radius=0.8):
        n_c = 0
        n_b = 0
        n_l = 0
        for gp in event.GenParticles:
            if deltaR(jet,gp)<=radius:
                if (abs(gp.pdgId()) == rt.b ): n_b+=1
                if (abs(gp.pdgId()) == rt.c ): n_c+=1
                if (abs(gp.pdgId()) <= rt.s ): n_l+=1
                if (abs(gp.pdgId()) == rt.g ): n_l+=1
        if (n_b >= 1): return rt.bb
        elif (n_c >= 1): return rt.cc
        else: return rt.light

    def GetQuarkComponent(self,event, jet, radius=0.8):
        n_c = 0
        n_b = 0
        n_l = 0
        for gp in event.GenParticles:
            if deltaR(jet,gp)<=radius:
                if (abs(gp.pdgId()) == rt.b ): n_b+=1
                if (abs(gp.pdgId()) == rt.c ): n_c+=1
                if (abs(gp.pdgId()) <= rt.s ): n_l+=1
                if (abs(gp.pdgId()) == rt.g ): n_l+=1
        return (n_c,n_b,n_l)

    def GetJetFlavor2(self,event, jet, radius=0.8):
        list_ = []
        for gp in event.GenParticles:
            if deltaR(jet,gp)<=radius:
                list_.append(abs(gp.pdgId()))
        return list_


    def GetFlavourComposition(self):
        vars = {}
        sels = {}
        print "year\tchannel\tmass\t",
        for decay in ["Hcc","Hbb","HWW","Hgg","Htautau", "Helse","nomatch"]:
            for flav in ["", "_bb","_cc","_LL"]:
                print decay+flav+"\t",
        print ""
        vars["tot"] = []
        for year, channel, mass in list(itertools.product(["2016","2017","2018"], self.Channels, self.MassPointsReduced)):
            sample = self.Signal+("_inv" if "inv" in channel else "")+"_M"+str(mass)+"_"+year
            if DoControl([""], year+channel+sample, channel, sample): continue
            # print year,channel,sample
            vars.setdefault(year,{}).setdefault(channel,{}).setdefault(sample,[])
            sels.setdefault(year,{}).setdefault(channel,{}).setdefault(sample,[])
            BR = {"tot":0, "WW4q":0, "ZZ4q":0, "WW":0, "ZZ":0, "ll":0, "tau":0, "light":0, "else":0}
            BR2 = {}
            for i in ["0","1","2","3","4"]:
                for flav in ["", "_bb","_cc","_LL"]:
                    BR2["Wc"+i+flav] = 0
                    BR2["Zc"+i+flav] = 0
            for filename in glob.glob(self.Path_STORAGE+year+"/Selection/Puppi/"+channel+"channel/nominal/workdir_Selection_"+sample+"/*.root"):
                f_ = ROOT.TFile(filename)
                t_ = f_.Get("AnalysisTree")
                i = 0
                notmatched = 0
                for ev in t_:
                    if ev.ZprimeCandidate.size()!=1 :
                        print "Unexpected number of ZprimeCandidate."
                        continue
                    for zp in ev.ZprimeCandidate:
                        jet = zp.H()
                        num = jet.btag_MassDecorrelatedDeepBoosted_probHcc()+jet.btag_MassDecorrelatedDeepBoosted_probZcc()
                        den = num+jet.btag_DeepBoosted_raw_score_qcd()
                        score = num/den if den!=0 else 0
                        decay = rt.ZprimeDecayToString(int(ev.HDecay))
                        if abs(zp.H().hadronFlavour()) == 5: flavour = rt.bb
                        elif abs(zp.H().hadronFlavour()) == 4: flavour = rt.cc
                        else: flavour = rt.light
                        store = decay
                        if flavour==rt.bb:
                            store += "_bb"
                        if flavour==rt.cc:
                            store += "_cc"
                        if flavour==rt.light:
                            store += "_LL"
                        if decay =="Helse":
                            BR["tot"] += 1
                            for gp in ev.GenParticles:
                                if abs(gp.pdgId())!=25: continue
                                d1 = gp.daughter(ev.GenParticles, 1)
                                d2 = gp.daughter(ev.GenParticles, 2)
                                dec = abs(d1.pdgId())
                                if dec==24 or dec ==23:
                                    d11 = abs(d1.daughter(ev.GenParticles, 1).pdgId())
                                    d12 = abs(d1.daughter(ev.GenParticles, 2).pdgId())
                                    d21 = abs(d2.daughter(ev.GenParticles, 1).pdgId())
                                    d22 = abs(d2.daughter(ev.GenParticles, 2).pdgId())
                                    nc = 0
                                    nl = 0
                                    for dd in [d11,d12,d21,d22]:
                                        if dd ==4: nc+=1
                                        if dd>=11 and dd<=15 ==4: nl+=1
                                    if nl ==2 and dec==23:
                                        isll = True
                                    elif nl ==1 and dec==24:
                                        isll = True
                                    else:
                                        isll = False
                                if dec==24:
                                    BR2["Wc"+str(nc)+store.replace(decay,"")] +=1
                                    if isll:
                                        BR["ll"] += 1
                                    elif d11 <= 4 and d12 <= 4 and d21 <= 4 and d22 <= 4:
                                        BR["WW4q"] += 1
                                    else: BR["WW"] += 1
                                elif dec==23:
                                    BR2["Zc"+str(nc)+store.replace(decay,"")] +=1
                                    if isll:
                                        BR["ll"] += 1
                                    elif d11 <= 4 and d12 <= 4 and d21 <= 4 and d22 <= 4:
                                        BR["ZZ4q"] += 1
                                    else:
                                        BR["ZZ"] += 1
                                elif dec==15:
                                    BR["tau"] += 1
                                elif dec<4 or dec==21:
                                    BR["light"] += 1
                                else:
                                    BR["else"] += 1
                        vars[year][channel][sample].append(store)
                        if score<0.8: continue
                        sels[year][channel][sample].append(store)
            print "TOT=", BR["tot"],
            for x in BR:
                print "match=",x, ":", round(BR[x]*100./BR["tot"],2),
            for x in BR2:
                print "match=",x, ":", round(BR2[x]*100./BR["tot"],2),
            print ""
            # print mass, i, notmatched
            var = vars[year][channel][sample]
            vars["tot"].extend(var)
            sel = sels[year][channel][sample]
            tot_var = len(var)
            tot_sel = len(sel)
            if tot_var==0: continue
            if tot_sel==0: continue
            print year, "\t", channel, "\t", mass, "\t",
            for decay in ["Hcc","Hbb","HWW","Hgg","Htautau", "Helse","nomatch"]:
                for flav in ["", "_bb","_cc","_LL"]:
                    skim_var = len(list(filter(lambda x: decay in x and flav in x , var)))
                    skim_sel = len(list(filter(lambda x: decay in x and flav in x , sel)))
                    print round(skim_var*100./tot_var,2), "\t", round(skim_sel*100./tot_sel,2), "\t", round(skim_sel*100./skim_var if skim_var!=0 else 0,2), " -- ",
            print ""
        for x in BR:
            print x, BR[x]*100/BR["tot"],
        print ""
        tot = len( vars["tot"])
        for decay in ["Hcc","Hbb","HWW","Hgg","Htautau", "Helse","nomatch"]:
            for flav in ["", "_bb","_cc","_LL"]:
                tot_df = len(list(filter(lambda x: decay in x and flav in x , vars["tot"])))
                print round(tot_df*100./tot,2), "\t",
        print ""

    def GetSFEffect(self):
        vars = {}
        for year, channel, mass in list(itertools.product(["2016", "2017"], self.Channels, self.MassPointsReduced)):
            sample = self.Signal+("_inv" if "inv" in channel else "")+"_M"+str(mass)+"_"+year
            if DoControl([""], year+channel+sample, channel, sample): continue
            vars.setdefault(year,{}).setdefault(channel,{}).setdefault(sample,[])
            f_ = ROOT.TFile(self.Path_STORAGE+year+"/SignalRegion/Puppi/"+channel+"channel/nominal/"+self.PrefixrootFile+"MC."+sample+"_noTree.root")
            hname = "Zprime_mass_transversal_rebin100" if "inv" in channel else "Zprime_mass_rebin100"
            nom = f_.Get("ZprimeCandidate_DeepAk8_ZHccvsQCD_MD_SR/"+hname)
            var_up = f_.Get("ZprimeCandidate_taggerSF_up_DeepAk8_ZHccvsQCD_MD_SR/"+hname)
            var_down = f_.Get("ZprimeCandidate_taggerSF_down_DeepAk8_ZHccvsQCD_MD_SR/"+hname)
            print year,channel,sample,
            print "NORM", abs(1-var_up.Integral()/nom.Integral())*100., abs(1-var_down.Integral()/nom.Integral())*100,
            print "MEAN", abs(1-var_up.GetMean()/nom.GetMean())*100., abs(1-var_down.GetMean()/nom.GetMean())*100,
            print "STD", abs(1-var_up.GetStdDev()/nom.GetStdDev())*100., abs(1-var_down.GetStdDev()/nom.GetStdDev())*100





def main():
    EHFC = EstimateHFlavorComponent()
    EHFC.GetFlavourComposition()
    # EHFC.GetSFEffect()

if __name__ == '__main__':
    main()
