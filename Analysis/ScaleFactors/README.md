# Note on SF.

**BTag SFs**:
- BTag ScaleFactors added from:
  - 2016: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP
  - 2016: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP
  - 2017: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
  - 2018: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18

- BTag Efficiencies calculated within the framework.


**Muon SFs**:
- IDs:
  - Trigger and Isolation:
  - Taken from https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/
- Reco:
  - 2016 https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016
  - 2017 https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017
  - 2018 https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018

**Electron SFs**:
- IDs:
 - 2016preVFP: /eos/cms/store/group/phys_egamma/UL_SF_EGammaFiles/UL16preVFP/Electrons/egammaEffi.txt_Ele_Loose_preVFP_EGM2D.root
  - 2016postVFP: /eos/cms/store/group/phys_egamma/UL_SF_EGammaFiles/UL16postVFP/Electrons/egammaEffi.txt_Ele_Loose_postVFP_EGM2D.root
  - 2017: /eos/cms/store/group/phys_egamma/UL_SF_EGammaFiles/UL17/Electrons/egammaEffi.txt_EGM2D_Loose_UL17.root
  - 2018: /eos/cms/store/group/phys_egamma/UL_SF_EGammaFiles/UL18/Electrons/egammaEffi.txt_Ele_Loose_EGM2D.root

- Trigger:
  - 2016: Relative to HLT_Ele27_WPTight_Gsf OR HLT_Photon175 OR HLT_Ele115_CaloIdVT_GsfTrkIdT
  - 2017: Relative to HLT_Ele35_WPTight_Gsf OR HLT_Photon200 OR HLT_Ele115_CaloIdVT_GsfTrkIdT
  - 2018: Relative to HLT_Ele32_WPTight_Gsf OR HLT_Photon200 OR HLT_Ele115_CaloIdVT_GsfTrkIdT
  - Take from http://indico.cern.ch/event/1146225/contributions/4835158/attachments/2429997/4160813/HLTSFsWprime%20.pdf
  - Relative to Loose ID

- Reco:
  - 2016preVFP: /eos/cms/store/group/phys_egamma/UL_SF_EGammaFiles/UL16preVFP/Electrons/RecoSFs/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root
  - 2016postVFP: /eos/cms/store/group/phys_egamma/UL_SF_EGammaFiles/UL16postVFP/Electrons/RecoSFs/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root
  - 2017: /eos/cms/store/group/phys_egamma/UL_SF_EGammaFiles/UL17/Electrons/RecoSFs/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root
  - 2018: /eos/cms/store/group/phys_egamma/UL_SF_EGammaFiles/UL18/Electrons/RecoSFs/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root




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
