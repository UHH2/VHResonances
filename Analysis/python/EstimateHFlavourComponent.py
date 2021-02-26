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


    def GetJetFlavor2(self,event, jet, radius=0.8):
        list_ = []
        for gp in event.GenParticles:
            if deltaR(jet,gp)<=radius:
                list_.append(abs(gp.pdgId()))
        return list_


    def GetFlavourComposition(self):
        vars = {}
        elses = {}
        for decay in ["Hcc","Hbb","HWW","Hgg","Helse","nomatch"]:
            for flav in ["", "_bb","_cc","_LL"]:
                print decay+flav+"\t",
        print ""
        for year, channel, mass in list(itertools.product(["2018"], self.Channels, self.MassPointsReduced)):
            sample = self.Signal+("_inv" if "inv" in channel else "")+"_M"+str(mass)+"_"+year
            if DoControl([""], year+channel+sample, channel, sample): continue
            # print year,channel,sample
            vars.setdefault(year,{}).setdefault(channel,{}).setdefault(sample,[])
            elses.setdefault(year,{}).setdefault(channel,{}).setdefault(sample,[])
            for filename in glob.glob(self.Path_STORAGE+year+"/Selection/Puppi/"+channel+"channel/nominal/workdir_Selection_"+sample+"/*.root"):
                f_ = ROOT.TFile(filename)
                t_ = f_.Get("AnalysisTree")
                i = 0
                for ev in t_:
                    if ev.ZprimeCandidate.size()!=1 :
                        print "Unexpected number of ZprimeCandidate."
                        continue
                    for zp in ev.ZprimeCandidate:
                        jet = zp.H()
                        # num = jet.btag_MassDecorrelatedDeepBoosted_probHcc()+jet.btag_MassDecorrelatedDeepBoosted_probZcc()
                        # num = jet.btag_DeepBoosted_probHcc()
                        # score = (num)/(num+jet.btag_DeepBoosted_raw_score_qcd())
                        # if score<0.8: continue
                        # vars[year][channel][sample].append(self.GetJetFlavor(ev,zp.H()))
                        decay = rt.ZprimeDecayToString(int(ev.HDecay))
                        flavour = self.GetJetFlavor(ev,zp.H())
                        # if decay=="Helse" and flavour==rt.cc:
                        #     elses[year][channel][sample].append(self.GetJetFlavor2(ev,zp.H()))
                        store = decay
                        if flavour==rt.bb:
                            store += "_bb"
                        if flavour==rt.cc:
                            store += "_cc"
                        if flavour==rt.light:
                            store += "_LL"
                        vars[year][channel][sample].append(store)
            var = vars[year][channel][sample]
            els = list(itertools.chain(*elses[year][channel][sample]))
            tot = len(var)
            if tot==0: continue
            print sorted(set(els))
            # tot_c = len(list(filter(lambda x: x ==rt.cc , var)))
            # tot_b = len(list(filter(lambda x: x ==rt.bb , var)))
            # print round(tot_c*100./tot,2), round(tot_b*100./tot,2), round((tot-tot_c-tot_b)*100./tot,2)
            for decay in ["Hcc","Hbb","HWW","Hgg","Helse","nomatch"]:
                for flav in ["", "_bb","_cc","_LL"]:
                    tot_df = len(list(filter(lambda x: decay in x and flav in x , var)))
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
