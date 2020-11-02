from Utils import *

'''
Module to extract Tagger info. The output will be used later on to make quick checks
Need to choose if you want the whole dataset or not via the isFast option.
'''
# TODO invisible channel not fully implemented yet


class Extractor(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.isFast = True
        self.isFast = False
        self.fName = "TaggerVariables"
        if not self.isFast: self.fName += "_all"
        self.Samples = filter(lambda x: self.MainBkg in x or self.Signal in x, self.Processes_Year_Dict["2016"]) # 2016 as default. They are all the same
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/"
        os.system("mkdir -p "+self.outdir)

    @timeit # circa 40 minutes
    def StoreVars(self):
        vars = {}
        for year, channel, sample in list(itertools.product(self.years, self.Channels, self.Samples)):
            if self.isFast and "electron" in channel: continue
            sample = sample.replace("2016",year)
            if DoControl([""], year+channel+sample, channel, sample): continue
            if "invisible" in channel: continue #TODO
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
                        jet = zp.H()
                        vars.setdefault("jet_pt",        []).append(jet.pt())
                        vars.setdefault("jet_QCDb",      []).append(jet.btag_DeepBoosted_probQCDb())
                        vars.setdefault("jet_QCDbb",     []).append(jet.btag_DeepBoosted_probQCDbb())
                        vars.setdefault("jet_QCDc",      []).append(jet.btag_DeepBoosted_probQCDc())
                        vars.setdefault("jet_QCDcc",     []).append(jet.btag_DeepBoosted_probQCDcc())
                        vars.setdefault("jet_QCDqq",     []).append(jet.btag_DeepBoosted_probQCDothers())
                        vars.setdefault("jet_Tbqq",      []).append(jet.btag_DeepBoosted_probTbqq())
                        vars.setdefault("jet_Tbcq",      []).append(jet.btag_DeepBoosted_probTbcq())
                        vars.setdefault("jet_Tbq",       []).append(jet.btag_DeepBoosted_probTbq())
                        vars.setdefault("jet_Tbc",       []).append(jet.btag_DeepBoosted_probTbc())
                        vars.setdefault("jet_Wqq",       []).append(jet.btag_DeepBoosted_probWqq())
                        vars.setdefault("jet_Wcq",       []).append(jet.btag_DeepBoosted_probWcq())
                        vars.setdefault("jet_Zcc",       []).append(jet.btag_DeepBoosted_probZcc())
                        vars.setdefault("jet_Zqq",       []).append(jet.btag_DeepBoosted_probZqq())
                        vars.setdefault("jet_Zbb",       []).append(jet.btag_DeepBoosted_probZbb())
                        vars.setdefault("jet_Hbb",       []).append(jet.btag_DeepBoosted_probHbb())
                        vars.setdefault("jet_Hcc",       []).append(jet.btag_DeepBoosted_probHcc())
                        vars.setdefault("jet_Hqqqq",     []).append(jet.btag_DeepBoosted_probHqqqq())
                        vars.setdefault("jet_TvsQCD",    []).append(jet.btag_DeepBoosted_TvsQCD())
                        vars.setdefault("jet_WvsQCD",    []).append(jet.btag_DeepBoosted_WvsQCD())
                        vars.setdefault("jet_ZvsQCD",    []).append(jet.btag_DeepBoosted_ZvsQCD())
                        vars.setdefault("jet_ZbbvsQCD",  []).append(jet.btag_DeepBoosted_ZbbvsQCD())
                        vars.setdefault("jet_ZHbbvsQCD", []).append(jet.btag_DeepBoosted_ZHbbvsQCD())
                        vars.setdefault("jet_HbbvsQCD",  []).append(jet.btag_DeepBoosted_HbbvsQCD())
                        vars.setdefault("jet_H4qvsQCD",  []).append(jet.btag_DeepBoosted_H4qvsQCD())
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
        df.to_pickle(self.outdir+self.fName+".pkl")

    def LoadVars(self):
        self.df = pd.read_pickle(self.outdir+self.fName+".pkl")
        print self.df

def main():
    extractor = Extractor()
    extractor.StoreVars()
    extractor.LoadVars()

if __name__ == '__main__':
    main()
