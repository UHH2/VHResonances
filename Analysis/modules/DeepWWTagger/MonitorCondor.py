import glob
import numpy as np

Collection = "HOTVR"
# Collection = "Puppi"
path = "/nfs/dust/cms/user/amalara//WorkingArea/File/NeuralNetwork/input_varariables/VHResonances/2017/"+Collection+"/leptonchannel/"

mode = "Vars"
mode = "Inputs"

Samples = {}
Samples["Vars"] = ["MC_ZZ", "MC_WZ", "MC_WZ_Zmatch", "MC_QCD_Pt15to30", "MC_QCD_Pt30to50", "MC_QCD_Pt50to80", "MC_QCD_Pt80to120", "MC_QCD_Pt120to170",
                   "MC_QCD_Pt2400to3200", "MC_QCD_Pt3200toInf", "MC_WJetsToQQ_HT600to800", "MC_WJetsToQQ_HT400to600", "MC_ZJetsToQQ_HT800toInf",
                   "MC_HZ_HiggsToWWZToLL", "MC_TTToHadronic", "MC_TTToSemiLeptonic", "MC_QCD_Pt170to300", "MC_QCD_Pt300to470", "MC_QCD_Pt470to600",
                   "MC_QCD_Pt600to800", "MC_QCD_Pt800to1000", "MC_QCD_Pt1000to1400", "MC_QCD_Pt1400to1800","MC_QCD_Pt1800to2400"]
Samples["Inputs"] = ["MC_HWW", "MC_Top", "MC_QCD", "MC_W", "MC_Z"]

memory = {}
time = {}

for sample in Samples[mode]:
    print sample
    memory[sample] = []
    time[sample] = []
    # for ext in ["*.o*","*.l*","*.e*"]:
    for ext in ["*.o*","*.l*"]:
        files = glob.glob(path+"MC_*/"+mode+"/log/*"+sample+ext)
        for f in files:
            lines = open(f).readlines()
            if ext =="*.e*":
                if len(lines)>0 :
                    print "\t", f
                    print "\t", lines
                continue
            for l in lines:
                if "Memory (MB)" in l:
                    memory[sample].append(int(l.split()[3]))
                if "'Root2Numpy'" in l or "'Preprocessing'" in l:
                    time[sample].append(float(l.split()[1]))
    print "\t memory:", np.min(memory[sample]), np.max(memory[sample])
    print "\t time:", len(time[sample]), np.max(time[sample]), np.sum(time[sample])/len(time[sample]), np.mean(time[sample]), np.std(time[sample])
