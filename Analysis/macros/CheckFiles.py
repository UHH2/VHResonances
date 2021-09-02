# import glob, os, sys, ROOT
from ModuleRunner import *
# ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT.kError) + ";")
# ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT. kFatal) + ";")

class CheckFileNumbers(ModuleRunner):
    def __init__(self, year, Collections, Channels, Systematics, Module):
        ModuleRunner.__init__(self,year)
        self.SetModule(Module, Collections, Channels, Systematics)


    def Count(self, doSFRAME=False):
        for collection, channel, syst,sample in self.SmartLoop(self.Samples):
            if DoControl([""],self.year+collection+channel+syst, channel, sample): continue
            folder = self.Module+"/"+collection+"/"+channel+"channel/"+syst+"/workdir_"+self.Module+"_"+sample
            config_path  = self.SubmitDir+folder
            storage_path = self.Path_STORAGE+year+"/"+folder
            sframeDir = self.Path_SFRAME+folder
            tot_xml = 0
            tot_root = 0
            for xml in glob.glob(config_path+"/*xml"):
                xml_name = xml.replace(config_path,"")
                if "Result" in xml_name: continue
                if self.Module in xml_name: continue
                tot_xml += 1
                number = xml_name.replace(sample,"").replace(".xml","").replace("/","").replace("_","")
                root_file = (sframeDir if doSFRAME else storage_path)+"/"+self.PrefixrootFile+("MC." if "MC" in xml_name else "DATA.")+sample+"_"+str(int(number)-1)+".root"
                root_file = root_file.replace("__","_")
                if not os.path.isfile(root_file):
                    print "FILE doesn't exist", root_file
                    print "mkdir -p", sframeDir, "; sframe_main",xml
                else:
                    ntuple = ROOT.TFile(str(root_file))
                    if ntuple.IsZombie() or ntuple.ReadKeys()==0:
                        # print root_file
                        print "sframe_main",xml
                    else: tot_root += 1
                    ntuple.Close()
            print folder, " "*(100-len(folder)) , "DONE. Counted: XML=", tot_xml, " "*(5-len(str(tot_xml))), "ROOT:", tot_root , " "*(5-len(str(tot_root))), ("MISSING: "+str(tot_xml-tot_root) if (tot_xml!=tot_root) else "")



if __name__ == '__main__':

    Years       = ["2016", "2017","2018"]
    Years       = ["2018"]
    Collections = ["Puppi"]
    Channels    = ["muon", "electron", "invisible"]
    # Channels    = ["invisible"]
    Systematics = ["nominal", "JEC_up", "JEC_down", "JER_up", "JER_down"]
    Systematics = ["nominal", "JEC_up", "JEC_down", "JER_up", "JER_down", "MuonScale_up", "MuonScale_down"]
    Systematics = ["nominal"]

    Module = "HEMIssueStudy"
    # Module = "Preselection"
    # Module = "Selection"
    # Module = "SignalRegion"

    for year in Years:
        CFN = CheckFileNumbers(year = year, Systematics = Systematics, Collections = Collections, Channels = Channels, Module = Module)
        CFN.Count()
        # CFN.Count(doSFRAME=True)
