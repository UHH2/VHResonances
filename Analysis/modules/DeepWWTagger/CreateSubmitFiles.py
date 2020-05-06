import os, argparse
from NtuplesHandler import *
from math import ceil

def CountFiles(sample,**kwargs):
    if len(kwargs)== 0:
        return len(glob("/nfs/dust/cms/user/amalara/WorkingArea/File/Analysis/2017/GenericCleaning/All/leptonchannel/workdir_GenericCleaning_"+sample+"/uhh2.AnalysisModuleRunner.MC."+sample+"*root"))
    else :
        # print "/nfs/dust/cms/user/amalara/WorkingArea/File/NeuralNetwork/input_varariables/VHResonances/2017/"+kwargs["Collection"]+"/leptonchannel/"+sample+"/Vars/TopJet"+kwargs["others"]+"_"+sample+"_*.npy"
        return len(glob("/nfs/dust/cms/user/amalara/WorkingArea/File/NeuralNetwork/input_varariables/VHResonances/2017/"+kwargs["Collection"]+"/leptonchannel/"+sample+"/Vars/TopJet"+kwargs["others"]+"_"+sample+"_*.npy"))

def CreateFile(path,filename,JobName="Test", executable="submit_Condor.sh", dic={}, nHours = 3, Memory=2, isGPU=False, isMaxwell=False):
    # print dic["outdir"]
    os.system("mkdir -p "+dic["outdir"]+"/log/")
    if not isMaxwell:
        lines = []
        lines.append('#HTC Submission File')
        if isGPU :
            # lines.append('requirements  = OpSysAndVer == "SL6"')
            lines.append('requirements  = ( OpSysAndVer == "CentOS7" || OpSysAndVer == "SL6")')
            lines.append('Request_GPUs  = 1')
        else:
            lines.append('requirements  = ( OpSysAndVer == "CentOS7" || OpSysAndVer == "SL6")')
        lines.append('universe          = vanilla')
        lines.append('notification      = Error')
        lines.append('notify_user       = andrea.malara@desy.de')
        lines.append('outdir            = '+dic["outdir"])
        lines.append('output            = $(outdir)/log/log_'+str(JobName)+'.o$(ClusterId).$(Process)')
        lines.append('error             = $(outdir)/log/log_'+str(JobName)+'.e$(ClusterId).$(Process)')
        lines.append('log               = $(outdir)/log/log_'+str(JobName)+'.l$(ClusterId).$(Process)')
        lines.append('RequestMemory     = '+str(int(Memory)*1024))
        lines.append('RequestDisk       = 1048576')
        lines.append('+RequestRuntime   = '+str(int(nHours*60*60)))
        #You need to set up sframe
        lines.append('getenv            = True')
        lines.append('stream_output     = True')
        lines.append('stream_error      = True')
        lines.append('environment       = "LD_LIBRARY_PATH_STORED='+os.environ.get('LD_LIBRARY_PATH')+'"')
        lines.append('JobBatchName      = '+str(JobName))
        lines.append('executable        = '+path+executable)
        lines.append('Queue')

        with open(path+filename, "w") as outputfile:
            for line in lines:
                outputfile.write(line+"\n")

        lines = []
        lines.append('#!/bin/bash')
        lines.append('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_STORED\n')
        lines.append(dic["command"])
        lines.append('echo "done"')

        with open(path+executable, "w") as outputfile:
            for line in lines:
                outputfile.write(line+"\n")
        os.system("chmod 777 "+path+executable)

    if isGPU:
        lines = []
        lines.append('#!/bin/bash')
        lines.append('#SBATCH --partition=cms-uhh') ##### partitions=cms-uhh,maxgpu,allgpu
        lines.append('#SBATCH --time=1-00:00:00')
        lines.append('#SBATCH --nodes=1')
        lines.append('#SBATCH --constraint=GPU')
        lines.append('#SBATCH --job-name  '+str(JobName))
        lines.append('#SBATCH --output    '+dic["outdir"]+'/log/sbatch-%N-%j.out')
        lines.append('#SBATCH --error     '+dic["outdir"]+'/log/sbatch-%N-%j.err')
        lines.append('#SBATCH --mail-type ALL')
        lines.append('#SBATCH --mail-user andrea.malara@desy.de')
        lines.append('\nsource ~/.setpaths\n')
        lines.append(dic["command"])
        lines.append('echo "done"')
        with open(path+executable.replace("Condor","Maxwell"), "w") as outputfile:
            for line in lines:
                outputfile.write(line+"\n")
        os.system("chmod 777 "+path+executable)


"""
USAGE
python CreateSubmitFiles.py --mode=Vars  -s
python CreateSubmitFiles.py --mode=Inputs -s
python CreateSubmitFiles.py --mode=Trainings -s

python CreateSubmitFiles.py --mode=Vars      --Collection=HOTVR -s
python CreateSubmitFiles.py --mode=Inputs    --Collection=HOTVR -s
python CreateSubmitFiles.py --mode=Trainings --Collection=HOTVR -s

python CreateSubmitFiles.py --mode=Vars      --Collection=Puppi -s
python CreateSubmitFiles.py --mode=Inputs    --Collection=Puppi -s
python CreateSubmitFiles.py --mode=Trainings --Collection=Puppi -s
"""


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', action='store_true', dest="submit")
    parser.add_argument('--mode', action='store', dest="mode")
    parser.add_argument('--Collection', action='store', dest="Collection",default="Puppi")
    args = parser.parse_args()
    if args.mode!="Vars" and args.mode!="Inputs" and args.mode!="Trainings":
        raise RuntimeError("Mode is not correct")
    print args
    MRB = ModuleRunnerBase()
    path = os.getcwd()+"/submitFiles/"+args.mode+"/"
    common_path = "/nfs/dust/cms/user/amalara/WorkingArea/File/NeuralNetwork/" if not MRB.isMaxwell else "/beegfs/desy/user/amalara/WorkingArea/File/NeuralNetwork/"
    os.system("mkdir -p "+path)
    # Max time per file in second
    Times = {
        "MC_ZZ":                    20,
        "MC_WZ":                    20,
        "MC_WZ_Zmatch":             20,
        "MC_QCD_Pt15to30":          2,
        "MC_QCD_Pt30to50":          2,
        "MC_QCD_Pt50to80":          2,
        "MC_QCD_Pt80to120":         4,
        "MC_QCD_Pt120to170":        20,
        "MC_WJetsToQQ_HT600to800":  800,
        "MC_WJetsToQQ_HT400to600":  150,
        "MC_ZJetsToQQ_HT800toInf":  2100,
        "MC_HZ_HiggsToWWZToLL":     700,
        "MC_TTToHadronic":          1000,
        "MC_TTToSemiLeptonic":      800,
        "MC_QCD_Pt170to300":        2000,
        "MC_QCD_Pt300to470":        2000,
        "MC_QCD_Pt470to600":        2100,
        "MC_QCD_Pt600to800":        2200,
        "MC_QCD_Pt800to1000":       1900,
        "MC_QCD_Pt1000to1400":      1900,
        "MC_QCD_Pt1400to1800":      1600,
        "MC_QCD_Pt1800to2400":      800,
        "MC_QCD_Pt2400to3200":      2000,
        "MC_QCD_Pt3200toInf":       2000,
        }
    # nHours = 4 if args.mode!="Trainings" else 12
    # Memory = 4 if args.mode!="Trainings" else 12
    nHours = 1 if args.mode!="Trainings" else 100
    Memory = 1 if args.mode!="Trainings" else 4
    year = "2017"
    list_submit = []
    list_local = []
    list_logfiles = []
    if args.mode!="Trainings":
        for Sample in MRB.Samples_matching:
            print Sample
            SubSamples = MRB.Samples_dict[Sample] if args.mode=="Vars" else [(x,y) for x in [["200", "4000"]] for y in ["", "_others"]]
            for SubSample in SubSamples:
                if args.mode=="Vars":
                    max = CountFiles(SubSample)
                    step = int(max/ceil(max*Times[SubSample]/(60.*60.*float(nHours/1.5))))
                else:
                    ptmin = SubSample[0][0]
                    ptmax = SubSample[0][1]
                    others = SubSample[1]
                    max = CountFiles(Sample,Collection=args.Collection,others=others)
                    step = int(60.*60.*float(nHours)/100) #TODO 100 is sec per event(it depends on size of file Check every time) Put more just to stay on time
                print " "*5, SubSample, " "*(30-len(str(SubSample))), "max:"+str(max)+" "*(3-len(str(max))), "\tstep:", step, "\tnfiles:", len(range(0,max, step)), "" if not args.mode=="Vars" else "\teta: "+str(step*Times[SubSample]/60./60.)
                for job in range(0,max, step):
                    JobName = year+"_"+args.mode+"_"+Sample
                    JobName += "_"+SubSample if args.mode == "Vars" else others+"_pt"+ptmin+"_"+ptmax
                    JobName += "_job"+str(job)
                    filename = "Condor_"+args.Collection+"_"+JobName+".submit"
                    executable="submit_Condor_"+args.Collection+"_"+JobName+".sh"
                    dic = {}
                    dic["outdir"] = common_path+"input_varariables/VHResonances/"+year+"/"+args.Collection+"/leptonchannel/"+Sample+"/"+args.mode+"/"
                    dic["command"] = "python "+os.getcwd()+"/NtuplesProduction.py --do"+args.mode+" --Samples="+Sample+" --year="+year
                    dic["command"] += " --SubSample="+SubSample if args.mode == "Vars" else " --ptmin="+ptmin+" --ptmax="+ptmax+" --others="+others
                    dic["command"] += " --firstfile="+str(job)+" --lastfile="+str(job+step)
                    dic["command"] += " --Collection="+args.Collection
                    CreateFile(path=path,filename=filename,JobName=JobName, executable=executable, dic=dic, nHours=nHours, Memory=Memory)
                    list_submit.append(["condor_submit",path+filename])
    else:
        os.system("python MakeModels.py --Collection="+args.Collection)
        for modelType in ["DCLModel"]:
        # for modelType in ["DCLModel","ConvModel","SequentialModel"]:
            mypath = path+args.Collection+"/"+modelType+"/"
            for JobName in glob(mypath+"*.json"):
                JobName = JobName.replace(".json","").replace(mypath,"")
                filename = "Condor_"+JobName+".submit"
                executable="submit_Condor_"+JobName+".sh"
                dic = {}
                dic["outdir"] = common_path+"trainings/VHResonances/"+year+"/"+args.Collection+"/leptonchannel/"+modelType+"/"+JobName+"/"
                dic["command"] = "python Training.py "+mypath+"/"+JobName+".json"
                CreateFile(path=mypath,filename=filename,JobName=JobName, executable=executable, dic=dic, nHours=nHours, Memory=Memory, isGPU=True, isMaxwell=MRB.isMaxwell)
                doFast = "10epoch" in JobName
                if MRB.isMaxwell and not doFast: continue
                if not MRB.isMaxwell and doFast: continue
                if "full" in executable: continue
                if "300epochs" in executable: continue
                # if "300kevents" in executable: continue
                list_submit.append(["condor_submit" if not MRB.isMaxwell else "sbatch", mypath+"/"+(filename if not MRB.isMaxwell else executable.replace("Condor","Maxwell") ) ])


    print "list_submit:", len(list_submit)
    print "list_local: ", len(list_local)
    for x in list_submit:
        print x
    if args.submit:
        print list_submit[-1]
        parallelise(list_submit, 20)
        parallelise(list_local, 20,list_logfiles)
