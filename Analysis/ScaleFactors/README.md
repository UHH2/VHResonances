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


**Electron SFs**:
- IDs:
  - 2016: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations (Fall17v2)
  - 2017: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations (Fall17v2)
  - 2018: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations (Fall17v2)




  Year | Run   | Lumi         | Command
  ---- | ----- | ------------ | -------
  2016 | B     | 5.711130445  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_272007-275376_13TeV_PromptReco_Collisions16_JSON_eraB.txt
  2016 | C     | 2.572903489  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_275657-276283_13TeV_PromptReco_Collisions16_JSON_eraC.txt
  2016 | D     | 4.242291557  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_276315-276811_13TeV_PromptReco_Collisions16_JSON_eraD.txt
  2016 | E     | 4.025228137  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_276831-277420_13TeV_PromptReco_Collisions16_JSON_eraE.txt
  2016 | F     | 3.104509132  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_277772-278808_13TeV_PromptReco_Collisions16_JSON_eraF.txt
  2016 | G     | 7.575824256  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_278820-280385_13TeV_PromptReco_Collisions16_JSON_eraG.txt
  2016 | H     | 8.650628380  | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Era/Prompt/Cert_280919-284044_13TeV_PromptReco_Collisions16_JSON_eraH.txt
  2016 | BCDEF | 19.65606276  | Sum
  2016 | GH    | 16.226452636 | Sum
  2016 | Tot   | 35.882515396 | brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt
