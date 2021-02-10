from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

ForThesis(TDR)

'''
Module to plot Systematics

- Need the Selection or SignalRegion output as input
- Specify in the steer file which years, channels and histFolders you want to run over.

'''
# TODO invisible channel not fully implemented yet
# TODO fix yaxis


class CheckFileNumbers(VariablesBase):
    def __init__(self, year="", channel="", collection="collection", systematic="", module=""):
        VariablesBase.__init__(self)
        self.folder = year+"/"+module+"/"+collection+"/"+channel+"/"+systematic+"/"
        self.module = module

    def Count(self):
        config_path  = self.Path_ANALYSIS+"config/SubmittedJobs/"+self.folder
        storage_path = self.Path_STORAGE+self.folder
        tot_xml = 0
        tot_root = 0
        for workdir in glob.glob(config_path+"work*"):
            file_dir = workdir.replace(config_path,storage_path)
            workdir_name = workdir.replace(config_path,"").replace("workdir_"+self.module+"_","")
            for xml in glob.glob(workdir+"/*xml"):
                xml_name = xml.replace(workdir,"")
                if "Result" in xml_name: continue
                if self.module in xml_name: continue
                tot_xml += 1
                number = xml_name.replace(workdir_name+"_","").replace(".xml","").replace("/","")
                number = int(number)
                root_file = file_dir+"/"+self.PrefixrootFile+("MC." if "MC" in xml_name else "DATA.")+workdir_name+"_"+str(number-1)+".root"
                # print number, xml, root_file
                if not os.path.isfile(root_file):
                    # print "FILE doesn't exist"
                    print "mkdir -p", file_dir.replace(self.Path_ANALYSIS,self.Path_SFRAME), "; sframe_main",xml
                else: tot_root += 1
        print self.folder, " "*(70-len(self.folder)) , "DONE. Counted: XML=", tot_xml, " "*(5-len(str(tot_xml))), "ROOT:", tot_root


# MC_DY_HT400to600_2016_1.xml
# uhh2.AnalysisModuleRunner.MC.MC_DY_HT400to600_2016.root
if __name__ == '__main__':
    years       = ["2016", "2017", "2018"]
    Channels    = ["muonchannel", "electronchannel", "invisiblechannel"]
    Collections = ["Puppi"]
    Systematics = ["nominal", "JEC_up", "JEC_down", "JER_up", "JER_down", "MuonScale_up", "MuonScale_down"]
    Modules     = ["Preselection", "Selection", "SignalRegion"]

    years       = ["2016", "2018"]
    Channels    = ["invisiblechannel"]
    # Systematics = ["nominal"]
    Modules     = ["Preselection"]

    for year in years:
        for channel in Channels:
            for collection in Collections:
                for syst in Systematics:
                    for module in Modules:
                        CFN = CheckFileNumbers(year=year, channel=channel, collection=collection, systematic=syst, module=module)
                        CFN.Count()
