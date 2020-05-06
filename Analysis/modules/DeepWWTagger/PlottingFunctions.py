import numpy as np
from math import *
import sys
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from keras.utils import to_categorical
from keras.callbacks import History, ModelCheckpoint, ReduceLROnPlateau

from sklearn import preprocessing
from sklearn.metrics import roc_curve, auc

hbins = 100

colorsample = {"Other Bkg" : "k", "Random" : "k",
               "Higgs": "tab:red", "QCD": "b", "Top": "tab:green", "DY": "tab:olive", "WJets": "m", "WZ": "tab:purple", "ZZ": "tab:cyan",
               "MC_HWW": "tab:red", "MC_Hbb": "tab:orange", "MC_Top": "b", "MC_QCD": "tab:green", "MC_W": "m", "MC_Z": "tab:cyan",
               "MC_DY1JetsToLL": "m", "MC_DY2JetsToLL": "tab:purple", "MC_DY3JetsToLL": "tab:olive", "MC_DY4JetsToLL": "yellow",
               "MC_TTTo2L2Nu": "tab:purple", "MC_TTbarSemiLeptonic": "lightblue", "MC_TTToHadronic": "turquoise"}

if "nfs" in os.getcwd():
    colorsample = {"Other Bkg" : "k", "Random" : "k",
                   "Higgs": "r", "QCD": "b", "Top": "g", "DY": "tab:olive", "WJets": "m", "WZ": "tab:purple", "ZZ": "tab:cyan",
                   "MC_HWW": "r", "MC_Hbb": "tab:orange", "MC_Top": "b", "MC_QCD": "g", "MC_W": "m", "MC_Z": "c",
                   "MC_DY1JetsToLL": "m", "MC_DY2JetsToLL": "tab:purple", "MC_DY3JetsToLL": "tab:olive", "MC_DY4JetsToLL": "yellow",
                   "MC_TTTo2L2Nu": "tab:purple", "MC_TTbarSemiLeptonic": "lightblue", "MC_TTToHadronic": "turquoise"}

# print colorsample

def PlotRoc1vsAll(TB, isLogy=True, show_figure=False, save_figure=True, name="ROC", mode="validation"):
    print TB.modelPath+name+".pdf"
    print TB.modelPath+name+".png"
    classes = to_categorical(np.arange(len(TB.SamplesNames)))
    plt.cla()
    plt.xticks( np.arange(0.1,1.1,0.1) )
    plt.grid(True, which="both")
    for Sample, i_sample in TB.SamplesNames.iteritems():
        fpr, tpr, thr = roc_curve(TB.GetArray(obj="labels",mode=mode)[:,i_sample], TB.GetArray(obj="pred",mode=mode)[:,i_sample])
        label = Sample+": auc = "+str(round(auc(fpr,tpr),3))
        if isLogy:
            plt.semilogy(tpr, fpr, label=label, color=colorsample[Sample])
        else:
            plt.plot(tpr, fpr, label=label, color=colorsample[Sample])
    x= np.linspace(1e-06, 1,1000)
    if isLogy:
        plt.semilogy(x, x, label = "Random classifier: auc = 0.5 ", color=colorsample["Random"], linestyle = "--")
    else:
        plt.plot(x, x, label = "Random classifier: auc = 0.5 ", color=colorsample["Random"], linestyle = "--")
    plt.xlim([-0.05, 1.05])
    plt.ylim([1e-06, 1e00])
    plt.xlabel("Signal efficiency")
    plt.ylabel("Background mistag rate")
    plt.legend(loc="lower right", shadow=True)
    if save_figure:
        if isinstance(save_figure, bool):
            plt.savefig(TB.modelPath+name+".pdf")
            plt.savefig(TB.modelPath+name+".png")
    if show_figure:
        plt.show()


def PlotRoc1vs1(TB, isLogy=True, show_figure=False, save_figure=True, name="ROC", mode="validation"):
    for Sample, i_sample in TB.SamplesNames.iteritems():
        plt.cla()
        plt.xticks( np.arange(0.1,1.1,0.1) )
        plt.grid(True, which="both")
        for bkg, i_bkg in TB.SamplesNames.iteritems():
            if bkg == Sample:
                continue
            mask = np.any(TB.GetArray(obj="labels",mode=mode)[:,tuple(sorted((i_sample, i_bkg)))], axis=1)
            if len(TB.GetArray(obj="labels",mode=mode)[mask][:,i_sample])==0:continue
            fpr, tpr, thr = roc_curve(TB.GetArray(obj="labels",mode=mode)[mask][:,i_sample], TB.GetArray(obj="pred",mode=mode)[mask][:,i_sample])
            label = bkg+": auc = "+str(round(auc(fpr,tpr),3))
            if isLogy:
                plt.semilogy(tpr, fpr, label=label, color=colorsample[bkg])
            else:
                plt.plot(tpr, fpr, label=label, color=colorsample[bkg])
        x= np.linspace(1e-06, 1,1000)
        if isLogy:
            plt.semilogy(x, x, label = "Random classifier: auc = 0.5 ", color=colorsample["Random"], linestyle = "--")
        else:
            plt.plot(x, x, label = "Random classifier: auc = 0.5 ", color=colorsample["Random"], linestyle = "--")
        plt.xlim([-0.05, 1.05])
        plt.ylim([1e-06, 1e00])
        plt.xlabel("Signal efficiency")
        plt.ylabel("Background mistag rate")
        plt.legend(loc="lower right", shadow=True)
        if save_figure:
            if isinstance(save_figure, bool):
                plt.savefig(TB.modelPath+name+"_"+Sample+".pdf")
                plt.savefig(TB.modelPath+name+"_"+Sample+".png")
        if show_figure:
            plt.show()


def PlotOutputs1D(TB, isLogy=True, show_figure=False, save_figure=True, name = "Outputs", norm=True):
    classes = to_categorical(np.arange(len(TB.SamplesNames)))
    for cl, i_cl in TB.SamplesNames.iteritems():
        plt.cla()
        plt.xticks( np.arange(0.1,1.1,0.1) )
        plt.grid(True, which="both")
        plt.xlabel("NN response")
        plt.ylabel("A.U.")
        for Sample, i_sample in TB.SamplesNames.iteritems():
            mask_train = np.all(TB.GetArray(obj="labels",mode="train")[:]==classes[i_sample], axis=1)
            mask_val = np.all(TB.GetArray(obj="labels",mode="validation")[:]==classes[i_sample], axis=1)
            train = TB.GetArray(obj="pred",mode="train")[mask_train][:,i_cl]
            val = TB.GetArray(obj="pred",mode="validation")[mask_val][:,i_cl]
            if len(train)==0: continue
            if len(val)==0: continue
            y_val, _ = np.histogram(val, bins=hbins)
            y_train, bins_train = np.histogram(train, bins=hbins)
            integral = ((bins_train[1] - bins_train[0])*sum(y_train))
            weights_train = np.ones_like(train)/float(len(train)) if norm else np.ones_like(train)
            y_train2, bins_train2, _ = plt.hist(train, weights=weights_train, alpha=0.5, bins=hbins, histtype="step", log=True, label="Training Sample,"+Sample, color = colorsample[Sample])
            integral_norm = ((bins_train2[1] - bins_train2[0])*sum(y_train2))
            plt.errorbar(0.5*(bins_train[1:] + bins_train[:-1]), y_val*(integral_norm/integral)*(len(train)/len(val)), yerr= (y_val**0.5)*(integral_norm/integral)*(len(train)/len(val)), fmt=".", label="Validation Sample,"+Sample, color=colorsample[Sample])
        if norm:
            plt.ylim([1e-06, 5e00])
        else:
            plt.ylim([1e-02, 1e7])
        plt.legend(loc="lower right", shadow=True, ncol=2, fontsize="small")
        if save_figure:
            if isinstance(save_figure, bool):
                plt.savefig(TB.modelPath+name+"_"+cl+".png")
                plt.savefig(TB.modelPath+name+"_"+cl+".pdf")
        if show_figure:
            plt.show()



def plot_losses(hist, show_figure=False, save_figure=True, losses="loss", min_epoch=0, name = "history" ):
    plt.clf()
    plt.plot(range(1,1+len(hist.history[losses][min_epoch:])), hist.history[losses][min_epoch:], label="Training "+losses)
    plt.plot(range(1,1+len(hist.history["val_"+losses][min_epoch:])), hist.history["val_"+losses][min_epoch:], label="Validation "+losses)
    plt.grid()
    plt.legend()
    # plt.yscale("log")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.xlim(left=0)
    # plt.ylim([np.array(hist.history[losses]+hist.history["val_"+losses]).min()*0.9,np.array(hist.history[losses]+hist.history["val_"+losses]).max()*1.1])
    if save_figure:
        if isinstance(save_figure, bool):
            plt.savefig(name+".png")
            plt.savefig(name+".pdf")
    if show_figure:
        plt.show()
