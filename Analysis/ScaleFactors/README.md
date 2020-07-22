# Note on SF.

**BTag SFs**:
- BTag ScaleFactors added from:
  - 2016: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy subjet_DeepCSV_2016LegacySF_V1.csv
  - 2017: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X subjet_DeepCSV_94XSF_V4_B_F.csv
  - 2018: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X subjet_DeepCSV_102XSF_V1.csv

- BTag Efficiencies calculated within the framework.


**Muon SFs**:
- IDs:
  - 2016: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2016LegacyRereco
  - 2016: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiB2G-MUO for TkrHighPt
  - 2017: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
  - 2018: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2018

- Trigger:
  - 2016: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonWorkInProgressAndPagResults (IsoMu24_OR_IsoTkMu24, Mu50_OR_TkMu50, BtoF and G)
  - 2017: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
  - 2018: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2018

- Reco:
  - A conservative systematic is set for PT > 300 GeV. It is recommended to add systematic uncertainties of 0.5% (PT < 300 GeV) or 1% (PT > 300 GeV).
  - Same for 2016, 2017 and 2018
  - More details here https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceSelectionAndCalibrationsRun2#Special_systematic_uncertainties

- Tracking:
  - An additional P-dependent systematic uncertainty to account for reconstruction with HighPt ID is used
  - Same for 2016, 2017 and 2018
  - More details here https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceSelectionAndCalibrationsRun2#Special_systematic_uncertainties

**Electron SFs**:
- IDs:
  - 2016: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations (Fall17v2)
  - 2017: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations (Fall17v2)
  - 2018: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations (Fall17v2)

- Trigger:
  - 2016: Relative to HLT_Ele27_WPTight_Gsf OR HLT_Photon175 OR HLT_Ele115_CaloIdVT_GsfTrkIdT
  - 2017: Relative to HLT_Ele35_WPTight_Gsf OR HLT_Photon200 OR HLT_Ele115_CaloIdVT_GsfTrkIdT
  - 2018: Relative to HLT_Ele32_WPTight_Gsf OR HLT_Photon200 OR HLT_Ele115_CaloIdVT_GsfTrkIdT
  - Take from Alberto's Analysis. Corresponding to what stated in https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgHLTScaleFactorMeasurements.
  - Relative to HEEPV70ID.

- Reco:
  - 2016: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations
  - 2017: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations
  - 2018: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations

- HLT:
  - 1% unc applyed for 2017 as described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaRunIIRecommendations
  - No effect on 2016 and 2018





**Lumi via brilcalc**:
https://twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart

  Year | Run   | Lumi         | Command | Notes
  ---- | ----- | ------------ | ------- | ------
  2016 | B     | 5.711130445  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_272007-275376_13TeV_PromptReco_Collisions16_JSON_eraB.txt |
  2016 | C     | 2.572903489  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_275657-276283_13TeV_PromptReco_Collisions16_JSON_eraC.txt |
  2016 | D     | 4.242291557  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_276315-276811_13TeV_PromptReco_Collisions16_JSON_eraD.txt |
  2016 | E     | 4.025228137  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_276831-277420_13TeV_PromptReco_Collisions16_JSON_eraE.txt |
  2016 | F     | 3.104509132  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_277772-278808_13TeV_PromptReco_Collisions16_JSON_eraF.txt |
  2016 | G     | 7.575824256  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_278820-280385_13TeV_PromptReco_Collisions16_JSON_eraG.txt |
  2016 | H     | 8.650628380  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_280919-284044_13TeV_PromptReco_Collisions16_JSON_eraH.txt |
  2016 | BCDEF | 19.65606276  | Sum |
  2016 | GH    | 16.226452636 | Sum |
  2016 | Tot   | 35.882515396 | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt |
  2018 | A     | 13.977900815 | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Era/Prompt/Cert_315252-316995_13TeV_PromptReco_Collisions18_JSON_eraA.txt |
  2018 | A_1   | 8.950818835  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Era/Prompt/Cert_315252-316995_13TeV_PromptReco_Collisions18_JSON_eraA.txt | before muon HLT update, run<316361
  2018 | A_2   | 5.027081980  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Era/Prompt/Cert_315252-316995_13TeV_PromptReco_Collisions18_JSON_eraA.txt | after muon HLT update, run>=316361
  2016 | Tot   | 59.832475339 | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt |
