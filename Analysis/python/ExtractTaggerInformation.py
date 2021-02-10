from Utils import *

'''
Module to extract Tagger info. The output will be used later on to make quick checks
Need to choose if you want the whole dataset or not via the isFast option.
'''


class Extractor(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.isFast = True
        self.isFast = False
        self.fName = "TaggerVariables_"+self.Signal
        self.fName = "TaggerVariables_"+self.MainBkg
        self.fName = "TaggerVariables_MC_TTbar"
        if not self.isFast: self.fName += "_all"
        # self.Samples = filter(lambda x: "MC_TTbar" in x or self.MainBkg in x or self.Signal in x, self.Processes_Year_Dict["2016"]) # 2016 as default. They are all the same
        if "MC_TTbar"   in self.fName: self.Samples = filter(lambda x: "MC_TTbar"   in x, self.Processes_Year_Dict["2016"]) # 2016 as default. They are all the same
        if self.MainBkg in self.fName: self.Samples = filter(lambda x: self.MainBkg in x, self.Processes_Year_Dict["2016"]) # 2016 as default. They are all the same
        if self.Signal  in self.fName: self.Samples = filter(lambda x: self.Signal  in x, self.Processes_Year_Dict["2016"]) # 2016 as default. They are all the same
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/"
        os.system("mkdir -p "+self.outdir)

    @timeit # circa 40 minutes
    def StoreVars(self):
        vars = {}
        for year, channel, sample in list(itertools.product(self.years, self.Channels, self.Samples)):
            if self.isFast and "electron" in channel: continue
            sample = sample.replace("2016",year)
            if DoControl([""], year+channel+sample, channel, sample): continue
            print year,channel,sample
            for filename in glob.glob(self.Path_STORAGE+year+"/Selection/Puppi/"+channel+"channel/nominal/workdir_Selection_"+sample+"/*.root"):
                f_ = ROOT.TFile(filename)
                t_ = f_.Get("AnalysisTree")
                i = 0
                for ev in t_:
                    if ev.ZprimeCandidate.size()!=1 :
                        print "Unexpected number of ZprimeCandidate."
                        continue
                    if i>1000 and self.isFast: break
                    if i%10000 == 0 and i!=0: print i,"out of", t_.GetEntriesFast(),  sample, channel, year
                    i +=1
                    for zp in ev.ZprimeCandidate:
                        vars.setdefault("year",[]).append(year)
                        vars.setdefault("channel",[]).append(channel)
                        vars.setdefault("sample",[]).append(sample)
                        vars.setdefault("weight_GLP",[]).append(ev.weight_GLP)
                        vars.setdefault("Zprime_mass",[]).append(zp.Zprime_mass())
                        vars.setdefault("Match",[]).append(rt.MatchingToString(rt.FloatToMatching(zp.discriminator("Match"))))
                        vars.setdefault("MatchingStatus",[]).append(rt.MatchingStatusToString(rt.FloatToMatching(zp.discriminator("MatchingStatus"))))
                        vars.setdefault("HDecay",[]).append(rt.ZprimeDecayToString(int(ev.HDecay)))
                        vars.setdefault("ZDecay",[]).append(rt.ZprimeDecayToString(int(int(ev.ZDecay))))
                        vars.setdefault("ZprimeDecay",[]).append(rt.ZprimeDecayToString(int(int(ev.ZprimeDecay))))
                        Z = zp.Z()
                        vars.setdefault("Z_pt",        []).append(Z.pt())
                        vars.setdefault("Z_DeltaR_ll", []).append(deltaR(zp.leptons()[0], zp.leptons()[1]))
                        jet = zp.H()
                        QCD = jet.btag_DeepBoosted_raw_score_qcd()
                        Hcc = jet.btag_DeepBoosted_probHcc()
                        Zcc = jet.btag_DeepBoosted_probZcc()
                        Hcc_MD = jet.btag_MassDecorrelatedDeepBoosted_probHcc()
                        Zcc_MD = jet.btag_MassDecorrelatedDeepBoosted_probZcc()
                        HccvsQCD = (Hcc)/(Hcc+QCD) if (Hcc+QCD)!=0 else 9999
                        HccvsQCD_MD = (Hcc_MD)/(Hcc_MD+QCD) if (Hcc_MD+QCD)!=0 else 9999
                        ZccvsQCD = (Zcc )/(Zcc+QCD) if (Zcc+QCD)!=0 else 9999
                        ZccvsQCD_MD = (Zcc_MD )/(Zcc_MD+QCD) if (Zcc_MD+QCD)!=0 else 9999
                        ZHccvsQCD = (Zcc+Hcc)/(Zcc+Hcc+QCD) if (Zcc+Hcc+QCD)!=0 else 9999
                        ZHccvsQCD_MD = (Zcc_MD+Hcc_MD)/(Zcc_MD+Hcc_MD+QCD) if (Zcc_MD+Hcc_MD+QCD)!=0 else 9999
                        vars.setdefault("jet_pt",           []).append(jet.pt())
                        vars.setdefault("jet_M",            []).append(jet.v4().M())
                        vars.setdefault("jet_SD",           []).append(jet.softdropmass())
                        vars.setdefault("jet_QCD",          []).append(QCD)
                        vars.setdefault("jet_QCD_MD",       []).append(jet.btag_MassDecorrelatedDeepBoosted_raw_score_qcd())
                        vars.setdefault("jet_Hcc",          []).append(Hcc)
                        vars.setdefault("jet_Hcc_MD",       []).append(Hcc_MD)
                        vars.setdefault("jet_Zcc",          []).append(Zcc)
                        vars.setdefault("jet_Zcc_MD",       []).append(Zcc_MD)
                        vars.setdefault("jet_Wqq",          []).append(jet.btag_DeepBoosted_probWqq())
                        vars.setdefault("jet_Wqq_MD",       []).append(jet.btag_MassDecorrelatedDeepBoosted_probWqq())
                        vars.setdefault("jet_Wcq",          []).append(jet.btag_DeepBoosted_probWcq())
                        vars.setdefault("jet_Wcq_MD",       []).append(jet.btag_MassDecorrelatedDeepBoosted_proWcq())
                        vars.setdefault("jet_Hbb",          []).append(jet.btag_DeepBoosted_probHbb())
                        vars.setdefault("jet_Hbb_MD",       []).append(jet.btag_MassDecorrelatedDeepBoosted_probHbb())
                        vars.setdefault("jet_Zbb",          []).append(jet.btag_DeepBoosted_probZbb())
                        vars.setdefault("jet_Zbb_MD",       []).append(jet.btag_MassDecorrelatedDeepBoosted_probZbb())
                        vars.setdefault("jet_bbvsL_MD",     []).append(jet.btag_MassDecorrelatedDeepBoosted_bbvsLight())
                        vars.setdefault("jet_ccvsL_MD",     []).append(jet.btag_MassDecorrelatedDeepBoosted_ccvsLight())
                        vars.setdefault("jet_Hqqqq",        []).append(jet.btag_DeepBoosted_probHqqqq())
                        vars.setdefault("jet_Hqqqq_MD",     []).append(jet.btag_MassDecorrelatedDeepBoosted_probHqqqq())
                        vars.setdefault("jet_ZbbvsQCD",     []).append(jet.btag_DeepBoosted_ZbbvsQCD())
                        vars.setdefault("jet_ZHbbvsQCD",    []).append(jet.btag_DeepBoosted_ZHbbvsQCD())
                        vars.setdefault("jet_HbbvsQCD",     []).append(jet.btag_DeepBoosted_HbbvsQCD())
                        vars.setdefault("jet_H4qvsQCD",     []).append(jet.btag_DeepBoosted_H4qvsQCD())
                        vars.setdefault("jet_HccvsQCD",     []).append(HccvsQCD)
                        vars.setdefault("jet_HccvsQCD_MD",  []).append(HccvsQCD_MD)
                        vars.setdefault("jet_ZccvsQCD",     []).append(ZccvsQCD)
                        vars.setdefault("jet_ZccvsQCD_MD",  []).append(ZccvsQCD_MD)
                        vars.setdefault("jet_ZHccvsQCD",    []).append(ZHccvsQCD)
                        vars.setdefault("jet_ZHccvsQCD_MD", []).append(ZHccvsQCD_MD)

                        vars.setdefault("subjet_size",   []).append(jet.subjets().size())
                        for j in range(0,2):
                            bb, b, lepb, uds, g, c = (-1,-1,-1,-1,-1,-1)
                            if j < jet.subjets().size():
                                sj = jet.subjets()[j]
                                bb = sj.btag_DeepFlavour_bb()
                                b = sj.btag_DeepFlavour_b()
                                lepb = sj.btag_DeepFlavour_lepb()
                                uds = sj.btag_DeepFlavour_uds()
                                g = sj.btag_DeepFlavour_g()
                                c = sj.btag_DeepFlavour_c()
                            vars.setdefault("subjet_"+str(j)+"_bb",   []).append(bb)
                            vars.setdefault("subjet_"+str(j)+"_b",    []).append(b)
                            vars.setdefault("subjet_"+str(j)+"_lepb", []).append(lepb)
                            vars.setdefault("subjet_"+str(j)+"_uds",  []).append(uds)
                            vars.setdefault("subjet_"+str(j)+"_g",    []).append(g)
                            vars.setdefault("subjet_"+str(j)+"_c",    []).append(c)

        df = pd.DataFrame(data=vars)
        joblib.dump(df, self.outdir+self.fName+".pkl")

    def LoadVars(self):
        self.df = pd.read_pickle(self.outdir+self.fName+".pkl")
        print self.df

def main():
    extractor = Extractor()
    extractor.StoreVars()
    extractor.LoadVars()

if __name__ == '__main__':
    main()
