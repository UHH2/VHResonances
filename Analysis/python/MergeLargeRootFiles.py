#!/usr/bin/env python

import ROOT, sys, glob

def MergeLargeRootFiles_(list_, file_output,size_max,noTree=False):
    ROOT.TTree.SetMaxTreeSize(size_max)
    # print "Updated tree size",ROOT.TTree.GetMaxTreeSize()/1E06
    rm = ROOT.TFileMerger(False)
    rm.SetFastMethod(True)
    rm.SetNotrees(noTree) # True for skipping tree
    rm.OutputFile(file_output)
    count = 0
    for F in list_:
        rm.AddFile(F)
        count += 1; sys.stdout.write( 'File added '+str(count)+' out of '+str(len(list_))+' \r'); sys.stdout.flush()
    sys.stdout.write( '\nMerge in progress\n')
    rm.Merge()
    sys.stdout.write( 'Merge complete\n')

def MergeLargeRootFiles(path, file_output,noTree=False):
    ext_ = ".root"
    if file_output[-5:] != ext_: file_output += ext_
    print "Merge files in:\t", path
    print "Output file:\t",file_output
    size_max = 200*(1024**3) # 200 Gb
    file_list = list(map(lambda x: x, glob.glob(path+"/*"+ext_) ))
    print len(file_list)
    size_ = 0
    file_list_list = []
    list_ = []
    for F in file_list:
        file_size = ROOT.TFile.Open(F).GetSize()
        if file_size > size_max: continue
        size_ += file_size
        list_.append(F)
        if size_>size_max and not noTree:
            if len(list_)>0: file_list_list.append(list_[:-1])
            list_= [F]
            size_ = 0
    if len(list_)>0: file_list_list.append(list_)
    for i, list_ in enumerate(file_list_list):
        extratext = "_"+str(i)+ext_ if len(file_list_list)!=1 else ext_
        MergeLargeRootFiles_(list_, file_output[:-5]+extratext,size_max,noTree=noTree)




if __name__ == "__main__":
    noTree = eval(sys.argv[1]) # True for skipping tree
    print noTree
    if len(sys.argv)>4 :
        print sys.argv[4:]
        # MergeLargeRootFiles_(sys.argv[3:], file_output[:-5]+extratext,size_max,noTree=noTree)
    else:
        MergeLargeRootFiles(sys.argv[2],sys.argv[3],noTree=noTree)
