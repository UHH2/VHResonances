from NtuplesHandler import *
from numpy import linalg as LA
from root_numpy import array2hist
from fast_histogram import histogram2d as fastHist2d

from sklearn import preprocessing
from sklearn.model_selection import train_test_split


class PrepareInputs(NtuplesHandler):
    def __init__(self, Channel="", Collection="", Samples=[], year="2017", pt_min=200,pt_max=4000, SubSetVars = [],extraText="", isTest=False):
        NtuplesHandler.__init__(self,Channel=Channel, Collection=Collection, Samples=Samples, year=year, extraText=extraText, isTest=isTest)
        if len(SubSetVars)==0:
            SubSetVars += ["jetpt", "jeteta", "jetphi", "jetenergy", "jettau1", "jettau2", "jettau3", "jettau4"]
            SubSetVars += ["PFpt", "PFeta", "PFphi", "PFenergy", "PFcharge","PFparticleID", "PFpuppi_weight"]
        self.filesize = 1000
        self.Radius = 1.0
        self.n_eta = 40
        self.n_phi = 40
        self.images_dict = {0:"PFpt", 1:"CH", 2:"NH"}
        self.pt_min=pt_min
        self.pt_max=pt_max
        self.ptrange = str(self.pt_min)+"_"+str(self.pt_max)
        self.doTraining = False
        self.doWeights = False
        self.isSubset = False
        self.Modes = ["Jet","PF"]
        self.SubSetVars = {}
        # self.SubSetVars["Jet"] = [var for var in SubSetVars if "jet" in var]
        # self.SubSetVars["PF"] = [var for var in SubSetVars if "PF" in var]
        # print len(self.SubSetVars["Jet"])
        # print len(self.SubSetVars["PF"])
        self.SubSetVars["Jet"] = dict([(var,(index, self.VarNames.index(var.replace("jet","TopJet")))) for index, var in enumerate(SubSetVars) if "jet" in var])
        self.SubSetVars["PF"] = dict([(var, (index,self.RootVarNameListsDict["PFParticles"].index(var.replace("PF",""))))  for index, var in enumerate(SubSetVars) if "PF" in var])
        self.SubSetVars["Image"] = {"PFpt":(0,0), "CH":(1,1), "NH": (2,2)}
        self.ObjNames = self.SubSetVars.keys()
        print "SubSetVars"
        prettydic(self.SubSetVars)
        self.Means = np.zeros(len(self.SubSetVars["Jet"]))
        self.Sigmas = np.ones(len(self.SubSetVars["Jet"]))
        self.isNormalized = False
    @timeit
    def Preprocessing(self, firstfile=None,lastfile=None, others=""):
        self.InputVars = {}
        Objects = ["TopJet","PFParticles"]
        for Sample in self.SamplesNames:
            a = os.system("mkdir -p "+self.FileStorageOutput+Sample+"/Inputs/")
            print Sample, others
            self.InputVars[Sample] = {}
            SubSample= self.Samples_dict[Sample][0]
            listFiles= sorted(glob(self.NHBDict[Sample][SubSample].outdir+"TopJet"+others+self.NHBDict[Sample][SubSample].FileExt.replace("index","*").replace(SubSample,Sample)))
            if self.isTest: listFiles= listFiles[:2]
            print len(listFiles), firstfile, lastfile
            listFiles= listFiles[firstfile:lastfile]
            print "Running on nfiles:", len(listFiles)
            Vars = {}
            for fname in listFiles:
                print fname
                index = fname[fname.rfind('_')+1:-4]
                savename = fname.replace("/Vars/","/Inputs/").replace("TopJet"+others,"Object").replace(Sample+"_",Sample+"_"+self.ptrange+"_")
                for obj in Objects:
                    Object = obj+others
                    Vars[Object] = np.load(fname.replace("TopJet"+others,Object)).astype(self.FS)
                    print "\tShape pre\t", Object, Vars[Object].shape
                # PT SELECTION
                PT_mask = (Vars["TopJet"+others][:,self.SubSetVars["Jet"]["jetpt"][1]]>self.pt_min)*(Vars["TopJet"+others][:,self.SubSetVars["Jet"]["jetpt"][1]]<self.pt_max)
                PT_mask*= ~(Vars["PFParticles"+others]==0).all(axis=(1,2))
                for obj in Objects:
                    Object = obj+others
                    Vars[Object] = Vars[Object][PT_mask]
                # CREATE IMAGE
                Object = "PFParticles"+others
                ImageName = "Image"+others
                # JET PREPROCESSING
                eta = Vars[Object][:, :, self.SubSetVars["PF"]["PFeta"][1]]
                phi = Vars[Object][:, :, self.SubSetVars["PF"]["PFphi"][1]]
                pt  = Vars[Object][:, :, self.SubSetVars["PF"]["PFpt"][1]]
                if len(eta)==0: raise RuntimeError("Check me! Unexpected input!")
                Vars[ImageName] = np.zeros((eta.shape[0],len(self.SubSetVars["Image"]),self.n_eta,self.n_phi),dtype=self.FS)
                # translation wrt the centre of "mass"
                eta -= np.expand_dims(np.sum(eta*pt,axis=1)/np.sum(pt,axis=1),axis=1)
                phi -= np.expand_dims(np.sum(phi*pt,axis=1)/np.sum(pt,axis=1),axis=1)
                # CALCULATE MATRIX OF INERTIA
                I = np.array([[np.sum(phi*phi*pt,axis=1), np.sum(-phi*eta*pt,axis=1)], [np.sum(-phi*eta*pt,axis=1), np.sum(eta*eta*pt,axis=1)]]).swapaxes(0,1).swapaxes(0,2).astype(self.FS)
                # condition of being diagonal
                mask = np.abs(I[:,0,1])>9*1e-03
                # CREATE ROTATION MATRIX
                RR = LA.eigh(I)
                # f1 corresponds to the highest EV
                f1 = np.zeros(len(I))
                f2 = np.zeros(len(I))
                H = np.zeros(len(I))
                f1[mask] = (RR[0][:,1][mask]-I[:,1,1][mask])/I[:,0,1][mask]
                f2[mask] = (RR[0][:,0][mask]-I[:,1,1][mask])/I[:,0,1][mask]
                f1[~mask] = np.amax(RR[0][~mask],axis=1)
                f2[~mask] = np.amin(RR[0][~mask],axis=1)
                H[mask] = 1./np.sqrt(f1[mask]*(f1[mask]-f2[mask]))
                Rot = np.zeros((I.shape[0],I.shape[1],I.shape[2]))
                Rot[:,0,0] = +H*f1
                Rot[:,0,1] = -H*f1*f2
                Rot[:,1,0] = +H
                Rot[:,1,1] = -H*f1
                Rot[~mask*(I[:,0,0]>=I[:,1,1])] = np.identity(2)
                Rot[~mask*(I[:,0,0]< I[:,1,1])] = np.identity(2)[:,(1,0)]
                # JET ROTATION
                angles = np.array((eta,phi)).swapaxes(0,1).swapaxes(1,2)
                angles = np.matmul(Rot, angles.swapaxes(1,2)).swapaxes(1,2)
                # DEFINE 1 QUADRANT TO BE THE ONE WITH HIGHER SUMPTself
                # THIS OPERATION IS NEEDED TO DEFINE 1 UNIQUE DIRECTION
                quad = []
                for el in range(angles.shape[0]):
                    quad1 = quad2 = quad3 = quad4 = 0
                    for pf in range(angles.shape[1]):
                        if angles[el,pf,0]>=0:
                            if angles[el,pf,1]>=0: quad1 += pt[el,pf]
                            else: quad4 += pt[el,pf]
                        else:
                            if angles[el,pf,1]>=0: quad2 += pt[el,pf]
                            else: quad3 += pt[el,pf]
                    quad.append(np.expand_dims(np.array([quad1,quad2,quad3,quad4]),axis=1))
                quad = np.concatenate(quad,axis=1).swapaxes(0,1)
                quad = np.argmax(quad,axis=1)+1
                angles[:,:,0][quad==2] *= -1.
                angles[:,:,1][quad==4] *= -1.
                angles[:,:,:][quad==3] *= -1.
                # STORE ETA AND PHI IN THE VECTORS
                Vars[Object][:, :, self.SubSetVars["PF"]["PFeta"][1]] = angles[:,:,0]
                Vars[Object][:, :, self.SubSetVars["PF"]["PFphi"][1]] = angles[:,:,1]
                # CREATE IMAGES
                PFpt  = Vars[Object][:,:,self.SubSetVars["PF"]["PFpt"][1]]
                PFeta = Vars[Object][:,:,self.SubSetVars["PF"]["PFeta"][1]]
                PFphi = Vars[Object][:,:,self.SubSetVars["PF"]["PFphi"][1]]
                PFparticleID   = Vars[Object][:,:,self.SubSetVars["PF"]["PFparticleID"][1]]
                PFpuppi_weight = Vars[Object][:,:,self.SubSetVars["PF"]["PFpuppi_weight"][1]]
                # (nevent, 3,40,40)
                for x in range(PFeta.shape[0]):
                    # PFparticleID[x]==1 charged hadron; PFparticleID[x]==5 neutral hadron;
                    H = fastHist2d(PFeta[x], PFphi[x], range = [[-self.Radius,self.Radius],[-self.Radius,self.Radius]], bins=(self.n_eta, self.n_phi), weights=PFpuppi_weight[x]*(PFparticleID[x]==1).astype(np.int)*PFpt[x]/np.cosh(PFeta[x]))
                    Vars[ImageName][x,0,:,:] = H
                    H = fastHist2d(PFeta[x], PFphi[x], range = [[-self.Radius,self.Radius],[-self.Radius,self.Radius]], bins=(self.n_eta, self.n_phi), weights=PFpuppi_weight[x]*(PFparticleID[x]==1).astype(np.int))
                    Vars[ImageName][x,1,:,:] = H
                    H = fastHist2d(PFeta[x], PFphi[x], range = [[-self.Radius,self.Radius],[-self.Radius,self.Radius]], bins=(self.n_eta, self.n_phi), weights=PFpuppi_weight[x]*((PFparticleID[x]==4)*(PFparticleID[x]==5)).astype(np.int))
                    Vars[ImageName][x,2,:,:] = H
                toSave = True
                for obj in Objects:
                    Object = obj+others
                    toSave *= (len(Vars[Object]) == len(Vars["Image"+others]))
                if not toSave:
                    print "Not saving"
                    continue
                for obj in Objects+["Image"]:
                    Object = obj+others
                    print savename.replace("Object",Object), Vars[Object].shape
                    np.save(savename.replace("Object",Object), Vars[Object].astype(self.FS))
    @timeit
    def MergeInputs(self, others=""):
        debug = False
        # debug = True
        for Sample in self.SamplesNames:
            SubSample= self.Samples_dict[Sample][0]
            index =0
            print "START MERGE", Sample, others
            os.system("mkdir -p "+self.NHBDict[Sample][SubSample].outdir.replace("/Vars/","/Inputs/temp/"))
            savename = self.NHBDict[Sample][SubSample].outdir.replace("/Vars/","/Inputs/")+"TopJet"+others+self.NHBDict[Sample][SubSample].FileExt.replace("index",self.ptrange+"_*").replace(SubSample,Sample)
            print savename
            filenameTemplate = sorted(glob(savename))
            savename = savename.replace("TopJet","/temp/TopJet")
            print savename
            nElements = 0
            TJs = []
            PFs = []
            Ims = []
            print "REMOVING", len(filenameTemplate), Sample, others
            for fn in filenameTemplate:
                TopJet = np.load(fn)
                if len(TopJet)>0:
                    print "\t", fn[fn.find("TopJet"):], len(filenameTemplate)
                    PF = np.load(fn.replace("TopJet","PFParticles"))
                    Image = np.load(fn.replace("TopJet","Image"))
                    nElements += len(TopJet)
                    TJs.append(TopJet)
                    PFs.append(PF)
                    Ims.append(Image)
                    if debug: print "nElements:", nElements
                    while nElements >= self.filesize:
                        if debug: print "nElements cicle:", nElements
                        TJs = np.concatenate(TJs)
                        PFs = np.concatenate(PFs)
                        Ims = np.concatenate(Ims)
                        if not debug: np.save(savename.replace("*",str(index)),TJs[:self.filesize].astype(self.FS))
                        if not debug: np.save(savename.replace("*",str(index)).replace("TopJet","PFParticles"),PFs[:self.filesize].astype(self.FS))
                        if not debug: np.save(savename.replace("*",str(index)).replace("TopJet","Image"),Ims[:self.filesize].astype(self.FS))
                        TJs = [TJs[self.filesize:]]
                        PFs = [PFs[self.filesize:]]
                        Ims = [Ims[self.filesize:]]
                        nElements = len(TJs[0])
                        if debug: print "nElements cicle post:", nElements
                        if debug: print "\tSaving index:", index
                        index += 1
                if not debug:
                    if os.path.exists(fn): os.remove(fn)
                    if os.path.exists(fn.replace("TopJet","PFParticles")): os.remove(fn.replace("TopJet","PFParticles"))
                    if os.path.exists(fn.replace("TopJet","Image")): os.remove(fn.replace("TopJet","Image"))
                if debug: print "nElements post:", nElements
            if debug: print "nElements last:", nElements
            if nElements!=0:
                TJs = np.concatenate(TJs)
                PFs = np.concatenate(PFs)
                Ims = np.concatenate(Ims)
                if not debug: np.save(savename.replace("*",str(index)),TJs.astype(self.FS))
                if not debug: np.save(savename.replace("*",str(index)).replace("TopJet","PFParticles"),PFs.astype(self.FS))
                if not debug: np.save(savename.replace("*",str(index)).replace("TopJet","Image"),Ims.astype(self.FS))
            print "MERGING", Sample, others
            for Object in ["TopJet","PFParticles", "Image"]:
                for el in glob(savename.replace("TopJet",Object)):
                    if not debug: os.system("mv "+el+" "+el.replace("/temp/","/"))
    def LoadNormalization(self):
        with open(self.FileStorageOutput+"NormInfo.txt", "r") as f:
            lines = f.readlines()
            Vars = self.SubSetVars["Jet"].keys()
            for line in lines[1:]:
                if not line.split()[0] in Vars:
                    continue
                self.Means[self.SubSetVars["Jet"][line.split()[0]][0]] = float(line.split()[2])
                self.Sigmas[self.SubSetVars["Jet"][line.split()[0]][0]] = float(line.split()[3])
    def ApplyNormalization1(self, vec):
        vec -= self.Means
        vec /= self.Sigmas
    def ApplyNormalization(self, vec,Var):
        index = self.SubSetVars["Jet"][Var][0]
        vec -= self.Means[index]
        vec /= self.Sigmas[index]
    def ApplyNormalization2(self, vec,Var):
        index = self.SubSetVars["Jet"][Var][0]
        vec -= self.Means[index]
        vec /= self.Sigmas[index]
        return vec
    def SetRanges(self,var):
        min, max, bin = (0,100,100)
        if "tags" in var:                                           min, max, bin = (-0.5,5.5,6)
        if "charge" in var:                                         min, max, bin = (-2,2,5)
        if "Area" in var:                                           min, max, bin = (0,10,100)
        if "Fraction" in var or "tau" in var or "weight" in var:    min, max, bin = (-0.05,1.05,100)
        if "btag" in var:                                           min, max, bin = (-1,1,100)
        if "eta" in var or "phi" in var:                            min, max, bin = (-math.pi,math.pi,100)
        if "pt" in var or "energy" in var:                          min, max, bin = (250,1250,100)
        if "PFpt" in var or "PFenergy" in var:                      min, max, bin = (0,1000,100)
        if "mass" in var:                                           min, max, bin = (0,200,100)
        if "particleID" in var:                                     min, max, bin = (-0.5,9.5,10)
        if "PFpt" in var or "PFenergy" in var:                      min, max, bin = (0,1000,1000)
        return min, max, bin
    @timeit
    def FileToHist(self, maxfilesize=100000, pt_cuts=[]):
        pt_cut_min =float(pt_cuts[0])
        pt_cut_max =float(pt_cuts[1])
        ptcut = str(int(pt_cut_min))+"_"+str(int(pt_cut_max))
        self.LoadNormalization()
        os.system("mkdir -p "+self.pathPlots+"/histos/")
        for Sample in self.SamplesNames:
            SubSample= self.Samples_dict[Sample][0]
            for others in ["", "_others"]:
                f = ROOT.TFile.Open(self.pathPlots+"/histos/"+Sample+others+"_"+self.ptrange+"_ptcut_"+ptcut+".root","RECREATE")
                savename = self.NHBDict[Sample][SubSample].outdir.replace("/Vars/","/Inputs/")+"TopJet"+others+self.NHBDict[Sample][SubSample].FileExt.replace("index",self.ptrange+"_*").replace(SubSample,Sample)
                filenameTemplate = sorted(glob(savename))
                savename = savename.replace("TopJet","/histos/TopJet")
                nElements = 0
                Objs = {"Jet":[], "PF":[], "Image":[]}
                for fn in filenameTemplate:
                    if nElements >= maxfilesize: break
                    TopJet = np.load(fn)
                    if len(TopJet)<=0: continue
                    nElements += len(TopJet)
                    PT_mask = (TopJet[:,self.SubSetVars["Jet"]["jetpt"][1]]>pt_cut_min)*(TopJet[:,self.SubSetVars["Jet"]["jetpt"][1]]<pt_cut_max)
                    Objs["Jet"].append(TopJet[PT_mask])
                    Objs["PF"].append(np.load(fn.replace("TopJet","PFParticles"))[PT_mask])
                    Objs["Image"].append(np.load(fn.replace("TopJet","Image"))[PT_mask])
                for Obj in self.ObjNames:
                    if len(Objs[Obj])<=0: continue
                    isImage = Obj=="Image"
                    Objs[Obj] = np.concatenate(Objs[Obj])
                    if len(Objs[Obj])<=0: continue
                    if Obj=="PF" : Objs[Obj] = np.concatenate(Objs[Obj],axis=0)
                    for VarName in self.SubSetVars[Obj]:
                        index = self.SubSetVars[Obj][VarName][1]
                        min, max, bin = self.SetRanges(VarName)
                        h_ = ROOT.TH1F(VarName, "; "+VarName+"; A.U.", bin, min, max) if not isImage else ROOT.TH2F("Image_"+VarName, "; #eta; #phi", self.n_eta, -self.Radius, self.Radius, self.n_phi, -self.Radius, self.Radius)
                        fill_hist(h_, Objs[Obj][:,index]) if not isImage else array2hist(np.sum(Objs[Obj][:,index,:,:], axis=0),h_)
                        h_.Write()
                        if Obj!="Jet": continue
                        min = self.ApplyNormalization2(min,VarName)
                        max = self.ApplyNormalization2(max,VarName)
                        self.ApplyNormalization(Objs[Obj][:,index],VarName)
                        h_ = ROOT.TH1F(VarName+"_norm", "; "+VarName+"; A.U.", bin, min, max) if not isImage else ROOT.TH2F("Image_"+VarName+"_norm", "; #eta; #phi", self.n_eta, -self.Radius, self.Radius, self.n_phi, -self.Radius, self.Radius)
                        fill_hist(h_, Objs[Obj][:,index]) if not isImage else array2hist(np.sum(Objs[Obj][:,index,:,:], axis=0),h_)
                        h_.Write()
                f.Close()
    @timeit
    def PlotInputVars(self, pt_cuts=[]):
        ptcut = pt_cuts[0]+"_"+pt_cuts[1]
        self.LoadNormalization()
        ROOT.gROOT.SetBatch(ROOT.kFALSE)
        ROOT.gROOT.SetBatch(ROOT.kTRUE)
        ROOT.gStyle.SetOptStat(1)
        ROOT.gStyle.SetOptStat(0)
        for Obj in self.ObjNames:
            isImage = Obj=="Image"
            for others in ["", "_others"]:
                isothers = "_others"==others
                os.system("mkdir -p "+self.pathPlots+Obj+others)
                for VarName in self.SubSetVars[Obj]:
                    for norm in ["", "_norm"]:
                        isNorm = "_norm"==norm
                        if Obj!="Jet" and isNorm: continue
                        min, max, bin = self.SetRanges(VarName)
                        if "jetpt"==VarName and "300" in self.ptrange : max/2
                        if "PFpt"==VarName or "PFenergy"==VarName : max = 300
                        if isNorm:
                            min = self.ApplyNormalization2(min,VarName)
                            max = self.ApplyNormalization2(max,VarName)
                        c_ = tdrCanvas(VarName+norm, min, max, 1.e-7, 1.e00, VarName, "A.U.")
                        c_.SetLogy(1)
                        # leg = tdrLeg(0.50,0.60,0.9,0.9, 0.03)
                        leg = tdrLeg(0.20,0.80,0.9,0.9, 0.03)
                        leg.SetNColumns(5)
                        # leg.SetHeader(VarName,"L")
                        for Sample in self.SamplesNames:
                            f = ROOT.TFile.Open(self.pathPlots+"/histos/"+Sample+others+"_"+self.ptrange+"_ptcut_"+ptcut+".root","r")
                            h_ = f.Get("Image_"+VarName+norm if isImage else VarName+norm)
                            if not h_: continue
                            if (h_.Integral()!=0):
                                h_.Scale(1./h_.Integral())
                            if isImage:
                                c_ = tdrCanvas("Image_"+Sample+others+VarName+norm, -self.Radius, self.Radius, -self.Radius, self.Radius, "#eta '", "#phi '", square=kSquare)
                                c_.SetLogz(1)
                                c_.SetRightMargin(0.15)
                                g_ = ROOT.TGraph2D(h_)
                                g_.SetNpx(200)
                                g_.SetNpy(200)
                                g_.GetHistogram().SetContour(2000)
                                g_.GetHistogram().SetMinimum(1e-04 if "NH"==VarName else 5*1e-06)
                                g_.GetHistogram().SetMaximum(5*1e-03 if "NH"==VarName else 5*1e-02)
                                g_.Draw("colz")
                                g_.SetDirectory(0)
                                c_name = self.pathPlots+Obj+others+"/"+"Image_"+Sample+others+"_"+VarName+norm+"_pt"+self.ptrange+"_ptcut_"+ptcut+".pdf"
                                c_.SaveAs(c_name, "pdf")
                            else:
                                tdrDraw(h_, "HIST", lcolor=self.colors[Sample], lstyle=1, fstyle=0)
                                h_.SetDirectory(0)
                                leg.AddEntry(h_, Sample ,"l")
                            f.Close()
                        if not isImage:
                            leg.SetLineColorAlpha(1, 0.7)
                            leg.Draw("same")
                            c_name = self.pathPlots+Obj+others+"/"+VarName+norm+others+"_pt"+self.ptrange+"_ptcut_"+ptcut+".pdf"
                            c_.SaveAs(c_name, "pdf")
