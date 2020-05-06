from TaggerBase import *
from numpy import exp, sqrt

def InvariantMass(E,pt,eta):
    t= exp(-eta)
    return sqrt(E**2-(pt*(1+t*t)/(2*t))**2)


class PostTrainingEvaluation(TaggerBase):
    def __init__(self, Channel="lepton", Collection="Puppi", extraText="", ModelInfoDict={}, doTraining=False, Object=None):
        TaggerBase.__init__(self,Channel=Channel,Collection=Collection, extraText=extraText,ModelInfoDict=ModelInfoDict, doTraining=doTraining)
        self.modelFolder = ModelInfoDict["modelName"]
        with open(self.OutputFolder+self.modelFolder+"/mymodelinfo.json") as f:
            ModelInfoDict = json.load(f)
        ModelInfoDict["doTraining"]=False
        if Object:
            self.Object = Object
        self.MakeModel()
        # self.LoadModel(modelName="mymodel")
        # self.LoadModel(modelName="models/bestmodel_epoch007_acc0.90")
        self.LoadModel(modelName="models/bestmodel_epoch002_acc0.94")
    def MakePredictions(self, load=False, step=20):
        # self.LoadAllData(load=load, step=step)
        self.Predict(step=step)
        # for x in ["train", "validation", "test"]:
        #     print "inputs", x, self.GetArray("inputs",x).shape
        #     print "pred", x, self.GetArray("pred",x).shape
        self.Plots(extraName="_post")
    def MaskSample(self,sample,mode="validation"):
        label = np.zeros((len(self.SamplesNames)),dtype=np.int)
        label[self.SamplesNames[sample]]=1
        label = [label]*len(data)
        mask = np.all(self.GetArray("labels",mode=mode)==label,axis=1)
        return mask
    def PlotsPostTraining(self):
        for cl, i_cl in self.Tagger.SamplesNames.iteritems():
            for Object in [""]:
                plt.cla()
                plt.xticks( np.arange(-0.1,1.1,0.1) )
                plt.grid(True, which="both")
                plt.xlabel("NN response")
                plt.ylabel("A.U.")
                for Sample, i_sample in self.SamplesNames.iteritems():
                    if Object == "SemiMatched":
                        sample = []
                        for obj in self.Tagger.Objects:
                            if obj == self.MatchingDict[Sample] or len(self.InputVars[Sample][obj])==0 or obj=="NotMatched" or len(self.Predictions[Sample][obj])==0: continue
                            sample.append(self.Predictions[Sample][obj])
                        if len(sample)==0: continue
                        sample = np.concatenate(sample)
                    else:
                        sample = self.Predictions[Sample][Object]
                    print "Check1", Object, cl, i_cl, Sample, sample.shape
                    if len(sample)==0: continue
                    # _, _, _ = plt.hist(sample[:,i_cl], bins=hbins, weights=weights_train, histtype="step", log=True, label="Test Sample,"+Sample+" "+Object, color = colorsample[Sample])
                    _, _, _ = plt.hist(sample[:], bins=hbins, histtype="step", log=True, label="Test Sample,"+Sample+" "+Object, color = colorsample[Sample])
                plt.xlim([-0.1, 1.1])
                plt.ylim([1e-01, 1e7])
                plt.legend(loc="upper center", shadow=True, fontsize="small")
                print "save", cl, Object, self.OutputFolder+self.modelFolder
                plt.savefig(self.OutputFolder+self.modelFolder+"/"+Object+"_"+cl+"_PostFit.png")
                plt.savefig(self.OutputFolder+self.modelFolder+"/"+Object+"_"+cl+"_PostFit.pdf")
