#!/usr/bin/env python
import argparse
# from NtuplesHandler import *
from PrepareInputs import *


# ##     ##  ######     ###     ######   ########
# ##     ## ##    ##   ## ##   ##    ##  ##
# ##     ## ##        ##   ##  ##        ##
# ##     ##  ######  ##     ## ##   #### ######
# ##     ##       ## ######### ##    ##  ##
# ##     ## ##    ## ##     ## ##    ##  ##
#  #######   ######  ##     ##  ######   ########

"""
USAGE
python CreateSubmitFiles.py --mode=Vars --Collection=HOTVR -s
python CreateSubmitFiles.py --mode=Vars --Collection=Puppi -s

# ./NtuplesProduction.py --doMergeVars --Collection=HOTVR

./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_HWW --others= &
./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_QCD --others= &
./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_Top --others= &
./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_W   --others= &
./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_Z   --others= &

./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_HWW --others=_others &
./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_QCD --others=_others &
./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_Top --others=_others &
./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_W   --others=_others &
./NtuplesProduction.py --doMergeVars --Collection=HOTVR --Sample=MC_Z   --others=_others &

./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_HWW --others= &
./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_QCD --others= &
./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_Top --others= &
./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_W   --others= &
./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_Z   --others= &

./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_HWW --others=_others &
./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_QCD --others=_others &
./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_Top --others=_others &
./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_W   --others=_others &
./NtuplesProduction.py --doMergeVars --Collection=Puppi --Sample=MC_Z   --others=_others &

python CreateSubmitFiles.py --mode=Inputs -s
python CreateSubmitFiles.py --mode=Inputs --Collection=HOTVR -s
python CreateSubmitFiles.py --mode=Inputs --Collection=Puppi -s

./NtuplesProduction.py --doMergeInputs

./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_HWW  --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_QCD  --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_Top  --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_W    --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_Z    --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_HWW  --others=_others &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_QCD  --others=_others &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_Top  --others=_others &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_W    --others=_others &
./NtuplesProduction.py --doMergeInputs --Collection=HOTVR --Sample=MC_Z    --others=_others &


./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_HWW  --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_QCD  --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_Top  --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_W    --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_Z    --others=  &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_HWW  --others=_others &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_QCD  --others=_others &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_Top  --others=_others &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_W    --others=_others &
./NtuplesProduction.py --doMergeInputs --Collection=Puppi --Sample=MC_Z    --others=_others &

./NtuplesProduction.py --doHistos     --Collection=HOTVR
./NtuplesProduction.py --doPlotImages --Collection=HOTVR

./NtuplesProduction.py --doHistos     --Collection=Puppi
./NtuplesProduction.py --doPlotImages --Collection=Puppi
"""



def NtuplesProductionMain(input_options):
    parser = argparse.ArgumentParser()
    for option in ["doVars","doMergeVars","doInputs","doMergeInputs","prepareNorm","doHistos","doPlotImages","isTest"]:
        parser.add_argument('--'+option,action='store_true',dest=option, help="Boolean option: simply type --"+option)
    parser.add_argument('--year',       action='store',     dest="year",        default="2017")
    parser.add_argument('--Samples',    action='append',    dest="Samples",     default=[])
    parser.add_argument('--SubSamples', action='store',     dest="SubSamples",  default="")
    parser.add_argument('--Channel',    action='store',     dest="Channel",     default="lepton")
    parser.add_argument('--Collection', action='store',     dest="Collection",  default="Puppi")
    parser.add_argument('--ptmin',      action='store',     dest="ptmin",       default="")
    parser.add_argument('--ptmax',      action='store',     dest="ptmax",       default="")
    parser.add_argument('--firstfile',  action='store',     dest="firstfile",   default="None")
    parser.add_argument('--lastfile',   action='store',     dest="lastfile",    default="None")
    parser.add_argument('--others',     action='store',     dest="others",      default="")

    args = parser.parse_args()

    print "\n****************************************"
    for arg in args.__dict__:
        print arg, " "*(15-len(arg)), args.__dict__[arg]
    print "****************************************\n"

    MRB = ModuleRunnerBase()

    if not args.Samples: args.Samples = MRB.Samples_matching
    print "Samples", args.Samples, args.SubSamples

    pts = [[args.ptmin,args.ptmax]] if args.ptmin!=args.ptmax else [["200", "4000"]]
    plotting_pts = pts+[["200", "300"],["300", "500"],["500", "4000"],["4000", "10000"]]
    print "pts", pts

    if args.doVars:
        NH = NtuplesHandler(Channel=args.Channel, Collection=args.Collection, Samples=args.Samples, year=args.year, extraText="", isTest=args.isTest)
        for Sample in args.Samples:
            print Sample
            for SubSample in MRB.Samples_dict[Sample]:
                if SubSample!=args.SubSamples and args.SubSamples!="": continue
                print "CreateVars:", SubSample
                NH.NHBDict[Sample][SubSample].ConvertRoot2Numpy(Sample=SubSample, firstfile=eval(args.firstfile),lastfile=eval(args.lastfile))
    if args.doMergeVars:
        NH = NtuplesHandler(Channel=args.Channel, Collection=args.Collection, Samples=args.Samples, year=args.year, extraText="", isTest=args.isTest)
        NH.MergeVars(others=args.others)
    if args.doInputs:
        PreIn = PrepareInputs(Channel=args.Channel, Collection=args.Collection, Samples=args.Samples, year=args.year, pt_min=int(args.ptmin),pt_max=int(args.ptmax), extraText="", isTest=args.isTest)
        PreIn.Preprocessing(firstfile=eval(args.firstfile),lastfile=eval(args.lastfile), others=args.others)
    if args.doMergeInputs:
        print args.Samples
        for Sample in args.Samples:
            for pt in pts:
                print Sample, pt
                PreIn = PrepareInputs(Channel=args.Channel, Collection=args.Collection, Samples=[Sample], year=args.year, pt_min=int(pt[0]),pt_max=int(pt[1]), extraText="", isTest=args.isTest)
                PreIn.MergeInputs(others=args.others)
    if args.doHistos:
        print "doHistos"
        for Sample in args.Samples:
            for pt in pts:
                print Sample, pt
                PreIn = PrepareInputs(Channel=args.Channel, Collection=args.Collection, Samples=[Sample], year=args.year, pt_min=int(pt[0]),pt_max=int(pt[1]), extraText="", isTest=args.isTest)
                for pt_cuts in plotting_pts:
                    PreIn.FileToHist(pt_cuts=pt_cuts)
    if args.doPlotImages:
        for pt in pts:
            print "doPlotImages", pt
            PreIn = PrepareInputs(Channel=args.Channel, Collection=args.Collection, Samples=args.Samples, year=args.year, pt_min=int(pt[0]),pt_max=int(pt[1]), extraText="", isTest=args.isTest)
            for pt_cuts in plotting_pts:
                PreIn.PlotInputVars(pt_cuts=pt_cuts)

    if args.prepareNorm:
        PreIn = PrepareInputs(Channel=args.Channel, Collection=args.Collection, Samples=args.Samples, extraText="",year=args.year, controls=[])
        PreIn.FitNorm()

            # PreIn.LoadNormalization()
            # PreIn.PlotInputVars(mode="PF", extraText="")
            # PreIn.PlotInputVars(mode="PF", extraText="_trasl")
            # PreIn = PrepareInputs(Channel=args.Channel, Collection=args.Collection, Samples=args.Samples, extraText="",year=args.year, controls=match)
            # PreIn.InputShape(Mode="Jet")
            # PreIn.InputShape(Mode="PF")
            # PreIn.InputShape(Mode="Images")
            # PreIn.PlotInputVars(mode=mode, extraText="")
            # PreIn.InputShape(Mode="Images")
            # PreIn.PlotInputVars(mode="Images", extraText="")
            # PreIn.PlotInputVars(mode="Images", extraText="_trasl")
            # if mode=="Matching":
            #     PreIn.ApplyNormalization(doObj=True)
            # if mode=="Inputs":
            #     PreIn.ApplyNormalization(doInput=True)
            # PreIn.PlotInputVars(mode=mode,extraText="_norm")

if __name__ == "__main__":
    #print 'Arguments',sys.argv[1:]
    status = NtuplesProductionMain(sys.argv[1:])
    exit(status)
