from ROOT import *
from glob import glob
import numpy as np
import os, sys, matplotlib, math
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# mode = "GaussFits"
mode = "CBFits"

dic_values = {}
dic_chi = {}

for fname in glob("../Limits/*/*/*/btag_DeepBoosted_H4qvsQCD/datacards/OutputFit_btag_DeepBoosted_H4qvsQCD.txt"):
    print fname
    with open(fname) as file_:
        lines = file_.readlines()
        for line in lines:
            if not "M1000" in line: continue
            print line
            if "param" in line:
                dic_values[line.split()[0]+(line.split()[1] if len(line.split())==5 else "")] = (float(line.split()[-1]),float(line.split()[-2]))
            if "chi2-pvalue" in line:
                dic_chi[line.split()[0]+(line.split()[1] if len(line.split())==5 else "")] = (float(line.split()[-1]),float(line.split()[-2]))

bad_values = []
values = []
for x in dic_values:
    var = dic_values[x][0]/dic_values[x][1]
    if var<10: values.append(var)
    else: bad_values.append((x,var))

values = np.array(values)

for x in bad_values:
    print x

print len(values), np.max(values), np.min(values)

plt.cla()
plt.yscale("log")
plt.hist(values, bins=100)
plt.savefig("FitValues_"+mode+".pdf")

chi2 = []
pvalue = []
bad_chi2 = []
bad_pvalue = []
for x in dic_chi:
    # print x, dic_chi[x]
    pvalue_ = dic_chi[x][0]
    chi2_ = dic_chi[x][1]
    if 0.2<chi2_ and chi2_<3: chi2.append(chi2_)
    else: bad_chi2.append((x,chi2_))
    if 5<pvalue_ and pvalue_<95: pvalue.append(pvalue_)
    else: bad_pvalue.append((x,pvalue_))

chi2 = np.array(chi2)
pvalue = np.array(pvalue)

print len(dic_chi), len(chi2), len(pvalue), len(bad_chi2), len(bad_pvalue)

for x in bad_chi2:
    print x

for x in bad_pvalue:
    print x

plt.cla()
# plt.yscale("log")
plt.hist(chi2, bins=100)
plt.savefig("chi2_"+mode+".pdf")

plt.cla()
# plt.yscale("log")
plt.hist(pvalue, bins=100)
plt.savefig("pvalue_"+mode+".pdf")
