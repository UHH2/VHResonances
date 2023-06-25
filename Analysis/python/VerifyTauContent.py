from Utils import *

# import tdrstyle_all as TDR
# TDR.writeExtraText = True
# TDR.extraText  = 'Simulation'
# TDR.extraText2 = 'Work in progress'

'''
Quick module to verify the content for signal.
It can also be used to verify the efficiencies of the cxx module.

- Mind which module you run it on.
- For each decay, the NEvents, SigEff and SigEff/BR are calculated.
- For the 'weights' step, the Eff/BR is expected to be 1.
- The output is print on screen. This is only a quick check.

'''


class VerifyTauContent(VariablesBase):
    def __init__(self, year = 'RunII', channel = 'muonchannel', collection='Puppi'):
        VariablesBase.__init__(self)
        self.year = year
        self.channel = channel.replace('channel','')
        self.collection = collection
        self.cuts = OrderedDict([
            ('Preselection', ['weights', 'NLeptonSel']),
            ('Selection',    ['Preselection', 'ScaleFactors']),
        ])
    
    def CalculateEfficiency(self, f_num, f_den, cut, scale=1000.):
        scale = float(scale)
        ref = f_den.Get('ZprimeCandidate_weights/sum_event_weights').GetBinContent(1)
        sum = 0
        for d in ['1','2']:
            h_ = f_num.Get('gen_'+cut+'/flavor_daughter'+d+'Z')
            for tau in [15, -15]:
                sum += h_.GetBinContent(h_.FindBin(tau))
        sum /=2 #double counting of daughters
        eff = str(round(scale*sum/ref,1))
        if len(str(scale))>6:
            eff += ' *1e-'+str(len(str(scale))-6)
        return eff

    def ReadFile(self):
        print self.year, self.channel, self.collection
        ref_module = 'Preselection'
        bkg_samples = ['DY']
        if 'invisible' in self.channel:
            bkg_samples += ['WJets']
        PrintFormattedLine(self.MassPoints+bkg_samples,5)
        for module, cuts in self.cuts.items():
            for cut in cuts:
                sums = []
                for mp in self.MassPoints:
                    massPoint = str(mp)
                    fname_base = self.Path_STORAGE+self.year+'/MODULE/'+self.collection+'/'+self.channel+'channel/nominal/'+self.PrefixrootFile+'MC.'+self.Signal+('_inv' if 'invisible' in channel else '')+'_M'+massPoint+'_'+self.year+'_noTree.root'
                    f_den = ROOT.TFile(fname_base.replace('MODULE',ref_module))
                    f_num = ROOT.TFile(fname_base.replace('MODULE',module))
                    sums.append(self.CalculateEfficiency(f_num, f_den, cut))
                    f_num.Close()
                    f_den.Close()
                for sample in bkg_samples:
                    fname_base = self.Path_STORAGE+self.year+'/MODULE/'+self.collection+'/'+self.channel+'channel/nominal/'+self.PrefixrootFile+'MC.MC_'+sample+'_'+self.year+'_noTree.root'
                    fname_den_base = fname_base.replace('MODULE',ref_module)
                    fname_num_base = fname_base.replace('MODULE',module)
                    if 'Preselection' in ref_module and not 'RunII' in self.year:
                        fname_den_base = fname_den_base.replace('.root', '_merge.root')
                    if 'Preselection' in module and not 'RunII' in self.year:
                        fname_num_base = fname_num_base.replace('.root', '_merge.root')
                    f_den = ROOT.TFile(fname_den_base)
                    f_num = ROOT.TFile(fname_num_base)
                    sums.append(self.CalculateEfficiency(f_num, f_den, cut, scale=10000000))
                    f_num.Close()
                    f_den.Close()
                PrintFormattedLine(sums,5)


if __name__ == '__main__':

    years = ['UL16preVFP','UL16postVFP','UL17','UL18', 'RunII']
    Channels = ['muonchannel','electronchannel', 'invisiblechannel']
    Collections = ['Puppi']

    years = ['RunII']
    # Channels = ['muonchannel']

    for year in years:
        for channel in Channels:
            for collection in Collections:
                VTC = VerifyTauContent(year=year, channel=channel, collection=collection)
                VTC.ReadFile()
