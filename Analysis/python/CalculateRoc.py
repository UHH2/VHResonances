from ROOT import *
from glob import glob
import numpy as np
import os, sys, matplotlib, math
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, auc
import ROOT

ROOT.gInterpreter.ProcessLine('#include "'+os.environ["CMSSW_BASE"]+'/src/UHH2/VHResonances/include/constants.hpp"')

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/macros/")
from ModuleRunnerBase import *

def SetCanvas():
    plt.cla()
    plt.xticks( np.arange(0.1,1.1,0.1) )
    plt.grid(True, which="both")
    x= np.linspace(1e-05, 1,1000)
    plt.semilogy(x, x, label = "Random classifier: auc = 0.5 ", linestyle = "--")
    plt.xlim([-0.05, 1.05])
    plt.ylim([min_bkg_reg, 5*1e00])
    plt.xlabel("Signal efficiency")
    plt.ylabel("Background mistag rate")

def sensitivity(s,b):
    # return s/math.sqrt(b)
    return math.sqrt(2*((s+b)*math.log(1+s/b)-s)) if b>0 else 0

def PlotRoc(varname,Sample, name=""):
    fpr = rocs[varname]["fpr"][Sample]
    tpr = rocs[varname]["tpr"][Sample]
    thr = rocs[varname]["thr"][Sample]
    sen = rocs[varname]["sen"][Sample]
    i_ = np.argwhere(thr > (0.8 if not "tau" in varname else 0.5) )
    i_ = i_[-1][0] if len(i_)>0 else -1
    i_sen = np.argwhere([x[0] for x in sen] > (0.8 if not "tau" in varname else 0.5) )
    i_sen = i_sen[-1][0] if len(i_sen)>0 else i_sen
    i_maxsen = np.argmax(np.array([x[1] for x in sen]))
    maxsen = sen[i_maxsen,1]
    i_maxsen = np.argwhere(thr > sen[i_maxsen,0] )
    i_maxsen = i_maxsen[-1][0] if len(i_maxsen)>0 else 0
    plt.plot(tpr[i_], fpr[i_], 'bo')
    if maxsen>0 and tpr[i_maxsen]<0.9: plt.plot(tpr[i_maxsen], fpr[i_maxsen], 'ro')
    name = varname.replace("btag_","").replace("DeepBoosted_","").replace("MassDecorrelated","MD_") if name=="" else name
    name = name + (" "*(15-len(name)))
    plt.semilogy(tpr, fpr, label=name+"auc=%0.2f sen=%0.2f@cut maxsen=%0.2f@%0.2ftpr" % (auc(fpr, tpr, reorder=applyweights),round(sen[i_sen][1],1), round(maxsen,1), round(tpr[i_maxsen],1)))
    # plt.semilogy(tpr, fpr, label=name+"auc=%0.2f cut=%0.2f sen=%0.2f" % (auc(fpr, tpr, reorder=applyweights),tpr[i_], tpr[i_maxsen]))


min_bkg_reg = 1e-03
min_sig_sen = 0.3

outputdir = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/python/Roc/"
os.system("mkdir -p "+outputdir+"samples/")

masses = [x.replace("MC_ZprimeToZH_M","") for x in ModuleRunnerBase().SignalSamples]
# masses = ["1000"]

discriminators= [ "tau1", "tau2", "tau3", "tau4", "tau21", "tau31", "tau41", "tau32", "tau42", "tau43",
                 "NN_HWW", "NN_Hbb", "NN_QCD", "NN_Top", "NN_W", "NN_Z", "NN_HWW_1", "NN_Hbb_1", "NN_QCD_1", "NN_Top_1", "NN_W_1", "NN_Z_1", "NN_HWW_2", "NN_Hbb_2", "NN_QCD_2", "NN_Top_2", "NN_W_2", "NN_Z_2",
                 "btag_DeepBoosted_H4qvsQCD", "btag_MassDecorrelatedDeepBoosted_H4qvsQCD", "btag_DeepBoosted_probHqqqq", "btag_DeepBoosted_raw_score_h", "btag_MassDecorrelatedDeepBoosted_probHqqqq",
                 "CNN_HWW", "CNN_QCD", "CNN_Top", "CNN_W", "CNN_Z", "DCL_HWW", "DCL_QCD", "DCL_Top", "DCL_W", "DCL_Z"]
best = [ "tau42", "btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_raw_score_h", "btag_DeepBoosted_probHqqqq", "btag_MassDecorrelatedDeepBoosted_H4qvsQCD", ]

discriminators= ["tau42", "btag_DeepBoosted_H4qvsQCD"]

best = [ "tau42", "btag_DeepBoosted_H4qvsQCD"]

mode = "pdf"
# mode = "png"

isfast = True
isfast = False

isLoad = False
# isLoad = True

applyweights = True
# applyweights = False

applyMatch = True
# applyMatch = False

extratext = "_weights" if applyweights else ""
extratext += "_match" if applyMatch else ""

years = ["2016", "2017", "2018"]
collections = ["Puppi"]
channels = ["muonchannel", "electronchannel"]
# years = ["2016"]
# collections = ["Puppi"]
# channels = ["muonchannel"]

for year in years:
    for collection in collections:
        for channel in channels:
            path = GenericPath().Path_STORAGE+year+"/Selection/"+collection+"/"+channel+"/nominal/"

            filename = {}
            filename["DY"] = path+"workdir_Selection_MC_DY_"+year+"/*.root"
            for mass in masses:
                filename["HWW_"+mass] = path+"workdir_Selection_MC_ZprimeToZH_M"+mass+"_"+year+"/*.root"
            if isLoad:
                for Sample in filename.keys():
                    print year, collection, channel, Sample
                    h_ = {}
                    for varname in discriminators+["weights", "Match", "MatchingStatus"]:
                        h_[varname] = []
                    for file in glob(filename[Sample]):
                        ntuple = TFile(str(file))
                        AnalysisTree = ntuple.Get("AnalysisTree")
                        for i_ev, ev in enumerate(AnalysisTree):
                            if isfast and i_ev>=100: break
                            h_["weights"].append(ev.weight_GLP)
                            for i_zp, zp in enumerate(ev.ZprimeCandidate):
                                if i_zp!=0:
                                    print "Number of Zprime > 1"
                                    continue
                                for varname in discriminators+["Match", "MatchingStatus"]:
                                    var = float(zp.discriminator(varname))
                                    h_[varname].append(var)
                    for varname in discriminators+["weights", "Match", "MatchingStatus"]:
                        h_[varname] = np.array(h_[varname])
                        np.save(outputdir+"samples/"+Sample+"_"+varname+"_"+year+"_"+collection+"_"+channel+".npy",h_[varname].astype(np.float32))
                continue

            rocs = {}
            for varname in discriminators:
                for Sample in sorted(list(set(filename.keys())-set(["DY"]))+["Allsignals"]):
                    rocs.setdefault(varname, {}).setdefault("fpr",{}).setdefault(Sample,[])
                    rocs.setdefault(varname, {}).setdefault("tpr",{}).setdefault(Sample,[])
                    rocs.setdefault(varname, {}).setdefault("thr",{}).setdefault(Sample,[])
                    rocs.setdefault(varname, {}).setdefault("sen",{}).setdefault(Sample,[])

            for Sample in sorted(list(set(filename.keys())-set(["DY"]))):
                print year, collection, channel, Sample
                SetCanvas()
                y_weights_Allsignals = {}
                y_score_Allsignals = {}
                y_true_Allsignals = {}
                for varname in discriminators:
                    temp_weights = np.load(outputdir+"samples/"+Sample+"_weights_"+year+"_"+collection+"_"+channel+".npy").astype(np.float32)
                    temp_score = np.load(outputdir+"samples/"+Sample+"_"+varname+"_"+year+"_"+collection+"_"+channel+".npy").astype(np.float32)
                    if applyMatch:
                        temp_Match = np.load(outputdir+"samples/"+Sample+"_Match_"+year+"_"+collection+"_"+channel+".npy").astype(np.float32)
                        temp_MatchingStatus = np.load(outputdir+"samples/"+Sample+"_MatchingStatus_"+year+"_"+collection+"_"+channel+".npy").astype(np.float32)
                        temp_weights = temp_weights[(temp_Match==ROOT.HWWMatch)*(temp_MatchingStatus==ROOT.Hadronic)]
                        temp_score = temp_score[(temp_Match==ROOT.HWWMatch)*(temp_MatchingStatus==ROOT.Hadronic)]
                    temp_true = np.ones(len(temp_score))
                    if not varname+"DY" in y_score_Allsignals:
                        y_weights_Allsignals[varname+"DY"]  = np.load(outputdir+"samples/DY_weights.npy").astype(np.float32)
                        y_score_Allsignals[varname+"DY"]    = np.load(outputdir+"samples/"+"DY"+"_"+varname+".npy").astype(np.float32)
                        y_true_Allsignals[varname+"DY"]     = np.zeros(len(y_score_Allsignals[varname+"DY"]))
                    DY_weights = y_weights_Allsignals[varname+"DY"]
                    DY_score = y_score_Allsignals[varname+"DY"]
                    DY_true = y_true_Allsignals[varname+"DY"]
                    y_weights = np.concatenate([temp_weights,DY_weights])
                    y_score = np.concatenate([temp_score,DY_score])
                    y_true = np.concatenate([temp_true,DY_true])
                    y_weights_Allsignals[varname+"AllSignals"]  = np.concatenate([y_weights_Allsignals[varname+"AllSignals"], y_weights]) if varname+"AllSignals" in y_weights_Allsignals else y_weights
                    y_score_Allsignals[varname+"AllSignals"]    = np.concatenate([y_score_Allsignals[varname+"AllSignals"], y_score]) if varname+"AllSignals" in y_score_Allsignals else y_score
                    y_true_Allsignals[varname+"AllSignals"]     = np.concatenate([y_true_Allsignals[varname+"AllSignals"], y_true]) if varname+"AllSignals" in y_true_Allsignals else y_true
                    fpr, tpr, thr = roc_curve(y_true, y_score, sample_weight=(y_weights if applyweights else None))
                    if auc(fpr, tpr, reorder=applyweights)<0.5:
                        fpr = 1 - fpr
                        tpr = 1 - tpr
                    mask  = fpr>min_bkg_reg
                    mask *= tpr>min_sig_sen
                    fpr = fpr[mask]
                    tpr = tpr[mask]
                    thr = thr[mask]
                    sen = np.array([(x,sensitivity(np.sum(temp_score[temp_score>x]*temp_weights[temp_score>x]),np.sum(DY_score[DY_score>x]*DY_weights[DY_score>x]))) for x in (thr if len(thr)<100 and len(thr)>1 else np.arange(0.05,0.95,0.01)) ])
                    rocs[varname]["fpr"][Sample] = fpr
                    rocs[varname]["tpr"][Sample] = tpr
                    rocs[varname]["thr"][Sample] = thr
                    rocs[varname]["sen"][Sample] = sen
                    PlotRoc(varname,Sample)
                for varname in discriminators:
                    y_weights = y_weights_Allsignals[varname+"AllSignals"]
                    y_score = y_score_Allsignals[varname+"AllSignals"]
                    y_true = y_true_Allsignals[varname+"AllSignals"]
                    fpr, tpr, thr = roc_curve(y_true, y_score, sample_weight=(y_weights if applyweights else None))
                    if auc(fpr, tpr, reorder=applyweights)<0.5:
                        fpr = 1 - fpr
                        tpr = 1 - tpr
                    mask  = fpr>min_bkg_reg
                    mask *= tpr>min_sig_sen
                    fpr = fpr[mask]
                    tpr = tpr[mask]
                    thr = thr[mask]
                    sen = np.array([(x,sensitivity(np.sum(temp_score[temp_score>x]*temp_weights[temp_score>x]),np.sum(DY_score[DY_score>x]*DY_weights[DY_score>x]))) for x in (thr if len(thr)<1000 and len(thr)>1 else np.arange(0.01,0.99,0.001)) ])
                    rocs[varname]["fpr"]["AllSignals"] = fpr
                    rocs[varname]["tpr"]["AllSignals"] = tpr
                    rocs[varname]["thr"]["AllSignals"] = thr
                    rocs[varname]["sen"]["AllSignals"] = sen


                plt.legend(loc="upper left", shadow=True, prop={'size': 6})
                plt.savefig(outputdir+"Roc_"+Sample+"_"+year+"_"+collection+"_"+channel+extratext+"."+mode)

                for case in ["tau", "NN","HWW","btag"]:
                    print "\t", case
                    SetCanvas()
                    for varname in discriminators:
                        if not case in varname: continue
                        PlotRoc(varname,Sample)
                    plt.legend(loc="upper left", shadow=True, prop={'size': 6})
                    plt.savefig(outputdir+"Roc_"+Sample+"_"+case+"_"+year+"_"+collection+"_"+channel+extratext+"."+mode)

                print "\t", "best"
                SetCanvas()
                for varname in best:
                    PlotRoc(varname,Sample)
                plt.legend(loc="upper left", shadow=True, prop={'size': 6})
                plt.savefig(outputdir+"Roc_"+Sample+"_best_"+year+"_"+collection+"_"+channel+extratext+"."+mode)


            for varname in best:
                print "\t", "allmpoints", varname
                SetCanvas()
                for Sample in sorted(list(set(filename.keys())-set(["DY"]))):
                    PlotRoc(varname,Sample,Sample)
                plt.legend(loc="upper left", shadow=True, prop={'size': 6})
                plt.savefig(outputdir+"Roc_"+varname+"_allmpoints_"+year+"_"+collection+"_"+channel+extratext+"."+mode)
            SetCanvas()
            for varname in best:
                PlotRoc(varname,"AllSignals")
            plt.legend(loc="upper left", shadow=True, prop={'size': 6})
            plt.savefig(outputdir+"Roc_"+varname+"_AllSignals_"+year+"_"+collection+"_"+channel+extratext+"."+mode)
