import os, glob

def ImportFunctions():
    for func in glob.glob("*.hpp"):
        # print func
        os.system("cp "+func+" "+os.environ["CMSSW_BASE"]+"/src/HiggsAnalysis/CombinedLimit/interface/")
    for func in glob.glob("*.cxx"):
        os.system("cp "+func+" "+os.environ["CMSSW_BASE"]+"/src/HiggsAnalysis/CombinedLimit/src/")

def ModifyClass():
    classes_name = os.environ["CMSSW_BASE"]+"/src/HiggsAnalysis/CombinedLimit/src/classes.h"
    classes_def_name = os.environ["CMSSW_BASE"]+"/src/HiggsAnalysis/CombinedLimit/src/classes_def.xml"
    with open(classes_name, "U") as classes:
        classes_lines = classes.readlines()
    with open(classes_def_name, "U") as classes_def:
        classes_def_lines = classes_def.readlines()
    classes_newlines = []
    SaveInClasses = {}
    SaveInClasses_def = {}
    for func in glob.glob("*.hpp"):
        name = func.replace(".hpp","")
        SaveInClasses[func] = True
        SaveInClasses_def[name] = True
    for line in classes_lines:
        for func in SaveInClasses:
            if func in line : SaveInClasses[func] = False
    for line in classes_def_lines:
        for name in SaveInClasses_def:
            if name in line : SaveInClasses_def[name] = False
    for func in SaveInClasses:
        if SaveInClasses[func]: os.system("sed -i '2i\#include \"HiggsAnalysis/CombinedLimit/interface/"+func+"\"\' "+classes_name)
    for name in SaveInClasses_def:
        if SaveInClasses_def[name]: os.system("sed -i '2i\        <class name=\""+name+"\" />\' "+classes_def_name)



def ImportToCombine():
    ImportFunctions()
    ModifyClass()
    os.system("source ./recompile.sh")

if __name__ == '__main__':
    ImportToCombine()
