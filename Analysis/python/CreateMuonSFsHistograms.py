from Utils import *
import math
from collections import OrderedDict


'''
More details here:
https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016
https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017
https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018
'''

class CreateMuonSFsHistograms():
    def __init__(self):
        self.dir = os.environ['CMSSW_BASE']+'/src/UHH2/VHResonances/Analysis/ScaleFactors/Muons/'
        self.fname = 'Muon_mode_SF_year.root'
        self.hname = 'mode_SF'

    def ReconstructionSF(self, mode='Reconstruction'):
        SFs = {
            'UL16': {
                'eta': OrderedDict([
                    ((0.0,1.6), {'pt': OrderedDict([ ((50,100),(0.9914,0.0008)), ((100,150),(0.9936,0.0009)), ((150,200),(0.9930,0.0010)), ((200,300),(0.993,0.002)), ((300,400),(0.9900,0.0040)), ((400,600),(0.990,0.003)), ((600,1500),(0.989,0.004)), ((1500,3500),(0.80,0.30))]),}),
                    ((1.6,2.4), {'pt': OrderedDict([ ((50,100),(1.0000,0.0000)), ((100,150),(0.9930,0.0010)), ((150,200),(0.9910,0.0010)), ((200,300),(0.985,0.001)), ((300,400),(0.9810,0.0020)), ((400,600),(0.979,0.004)), ((600,1500),(0.978,0.005)), ((1500,3500),(0.90,0.20))]),}),
                    ])
                },
            'UL17': {
                'eta': OrderedDict([
                    ((0.0,1.6), {'pt': OrderedDict([ ((50,100),(0.9938,0.0006)), ((100,150),(0.9950,0.0007)), ((150,200),(0.9960,0.0010)), ((200,300),(0.996,0.001)), ((300,400),(0.9940,0.0010)), ((400,600),(1.003,0.006)), ((600,1500),(0.987,0.003)), ((1500,3500),(0.90,0.10))]),}),
                    ((1.6,2.4), {'pt': OrderedDict([ ((50,100),(1.0000,0.0000)), ((100,150),(0.9930,0.0010)), ((150,200),(0.9890,0.0010)), ((200,300),(0.986,0.001)), ((300,400),(0.9890,0.0010)), ((400,600),(0.983,0.003)), ((600,1500),(0.986,0.006)), ((1500,3500),(1.01,0.01))]),}),
                    ])
                },
            'UL18': {
                'eta': OrderedDict([
                    ((0.0,1.6), {'pt': OrderedDict([ ((50,100),(0.9943,0.0007)), ((100,150),(0.9948,0.0007)), ((150,200),(0.9950,0.0009)), ((200,300),(0.994,0.001)), ((300,400),(0.9914,0.0009)), ((400,600),(0.993,0.002)), ((600,1500),(0.991,0.004)), ((1500,3500),(1.00,0.10))]),}),
                    ((1.6,2.4), {'pt': OrderedDict([ ((50,100),(1.0000,0.0000)), ((100,150),(0.9930,0.0010)), ((150,200),(0.9900,0.0010)), ((200,300),(0.988,0.001)), ((300,400),(0.9810,0.0020)), ((400,600),(0.983,0.003)), ((600,1500),(0.978,0.006)), ((1500,3500),(0.98,0.03))]),}),
                    ])
                },
            }
        for year in ['UL16', 'UL17', 'UL18']:
            fname = self.fname.replace('mode',mode).replace('year',year)
            hname = self.hname.replace('mode',mode)
            hist = ROOT.TH2D(hname,hname,3,0,2.4,69,50,3500) # arbitrary choice of pt bins of 50GeV to simplify loop
            for eta in range(1, hist.GetNbinsX()+1):
                x  = hist.GetXaxis().GetBinCenter(eta)
                dx = hist.GetXaxis().GetBinWidth(eta)/2
                df = ref_eta_bin = 0
                for eta_bin, pts in SFs[year]['eta'].items():
                    if round(x-dx,2)>=round(eta_bin[0],2) and round(x+dx,2)<=round(eta_bin[1],2):
                        df = pts
                        ref_eta_bin = eta_bin
                        break
                # print x-dx, x+dx, ref_eta_bin
                for pt in range(1, hist.GetNbinsY()+1):
                    y  = hist.GetYaxis().GetBinCenter(pt)
                    dy = hist.GetYaxis().GetBinWidth(pt)/2
                    sf = ref_pt_bin = 0
                    for pt_bin, sfs in df['pt'].items():
                        if round(y-dy,2)>=round(pt_bin[0],2) and round(y+dy,2)<=round(pt_bin[1],2):
                            sf = sfs
                            ref_pt_bin = pt_bin
                            break
                    # print " ", y-dy, y+dy, ref_pt_bin, sf
                    hist.SetBinContent(eta,pt,sf[0])
                    hist.SetBinError(eta,pt,sf[1])
            file_ = ROOT.TFile(self.dir+fname, 'RECREATE')
            hist.Write()
            file_.Close()


def main():
    SFCreator = CreateMuonSFsHistograms()
    SFCreator.ReconstructionSF()



if __name__ == '__main__':
    main()
