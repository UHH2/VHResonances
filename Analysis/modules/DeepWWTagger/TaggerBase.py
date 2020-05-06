from PrepareInputs import *
from Models import *
import json

from keras.models import model_from_json, load_model
from keras.utils import to_categorical, plot_model

from ExportModel import *

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


class TaggerBase(PrepareInputs):
    def __init__(self, Channel="lepton", Collection="Puppi", extraText="", ModelInfoDict={}, doTraining=True):
        PrepareInputs.__init__(self,Channel=Channel,Collection=Collection, SubSetVars = ModelInfoDict["SubSetVars"], extraText=extraText)
        self.OutputFolder =  self.FileStorageOutput.replace("input_varariables","trainings")+ModelInfoDict["modelType"]+"/"
        self.UpdateVars(ModelInfoDict, doTraining=doTraining)
        self.LoadNormalization()
        np.random.seed(4444)
    def UpdateVars(self, ModelInfoDict,doTraining=True):
        self.doTraining = ModelInfoDict["doTraining"]
        if not doTraining: self.doTraining = False
        self.isFitGenerator = ModelInfoDict["isFitGenerator"]
        # self.isFitGenerator = True
        # self.isFitGenerator = False
        self.pt_cut_min = ModelInfoDict["pt_min"]
        self.pt_cut_max = ModelInfoDict["pt_max"]
        self.isSubset = ModelInfoDict["isSubset"]
        self.isSubset = ModelInfoDict["isSubset"]
        self.datasetsize = ModelInfoDict["max_size"]
        self.modelPath = self.OutputFolder+ModelInfoDict["modelName"]+"/"
        self.modelName = ModelInfoDict["modelName"]
        self.modelType = ModelInfoDict["modelType"]
        if self.modelType == "SequentialModel" :
            self.Object = "TopJet"
        if self.modelType == "ConvModel" :
            self.Object = "Image"
        if self.modelType == "DCLModel" :
            self.Object = ["TopJet","PFParticles","Image"]
        if type(self.Object)==type(list()):
            self.inputcolumns = {}
            for var in self.Object:
                self.inputcolumns[var] = tuple([y for x,y in sorted(dict((x[0], x[1]) for x in self.SubSetVars[var.replace("Top","").replace("Particles","")].values()).iteritems())])
        else:
            self.inputcolumns = tuple([y for x,y in sorted(dict((x[0], x[1]) for x in self.SubSetVars[self.Object.replace("Top","").replace("Particles","")].values()).iteritems())])
        self.params = ModelInfoDict["params"]
        self.params["ConvLayer"] = ModelInfoDict["ConvLayer"]
        self.params["kernel_size"] = ModelInfoDict["ConvLayer"]
        self.params["kernel_size"] = ModelInfoDict["kernel_size"]
        self.params["strides"] = ModelInfoDict["strides"]
        self.params["DenseLayer"] = ModelInfoDict["DenseLayer"]
        self.params["LSTMLayer"] = ModelInfoDict["LSTMLayer"]
        if len(ModelInfoDict["SamplesNames"])!=0:
            self.SamplesNames = dict((name, k) for k,name in enumerate(ModelInfoDict["SamplesNames"]))
        print "SamplesNames"
        prettydic(self.SamplesNames)
        if not os.path.exists(self.modelPath):
            os.makedirs(self.modelPath)
        if self.doTraining:
            for file in glob(self.modelPath+"*"):
                if not os.path.isdir(file):
                    os.remove(file)
            with open(self.modelPath+"mymodelinfo.json","w") as f:
                f.write(json.dumps(ModelInfoDict))
            with open(self.modelPath+"mymodelinfo.txt","w") as f:
                for el in ModelInfoDict:
                    f.write(str(el)+"\t"+str(ModelInfoDict[el]).replace("u'","").replace("'","")+"\n")
    @timeit
    def MakeModel(self):
        self.datasets = {"train": {}, "validation": {}, "test": {} }
        for mode in self.datasets :
            for sample in self.SamplesNames :
                self.datasets[mode][sample] = []
        for sample in self.SamplesNames :
            subsample = self.Samples_dict[sample][0]
            listFiles = self.NHBDict[sample][subsample].outdir+"TopJet"+self.NHBDict[sample][subsample].FileExt.replace("index","*").replace(subsample,sample)
            listFiles = listFiles.replace("/Vars/","/Inputs/").replace(sample+"_",sample+"_"+self.ptrange+"_")
            if type(self.Object)!=type(list()):
                listFiles = listFiles.replace("/Vars/","/Inputs/").replace("TopJet",self.Object)
            listFiles = sorted(glob(listFiles))
            if self.datasetsize>0:
                listFiles = listFiles[:int(floor(self.datasetsize/self.filesize))]
            print sample, len(listFiles), listFiles[0]
            train_percentage = int(floor(len(listFiles)*90./100))
            self.datasets["train"][sample] += listFiles[:train_percentage]
            listFiles = listFiles[train_percentage:]
            val_percentage = int(floor(len(listFiles)*50./100))
            self.datasets["validation"][sample] += listFiles[:val_percentage]
            self.datasets["test"][sample] += listFiles[val_percentage:]
        for mode in self.datasets :
            min_ = np.min(np.array([len(self.datasets[mode][sample]) for sample in self.SamplesNames]))
            for sample in self.SamplesNames :
                self.datasets[mode][sample] = self.datasets[mode][sample][:min_]
                dataset = len(self.datasets[mode][sample])
                print mode, " "*(12-len(str(mode))), sample, " "*(12-len(str(sample)))+"nFiles:", dataset, " "*(5-len(str(dataset)))+"nEvents:", dataset*self.filesize
        if self.doTraining:
            # TODO shuffle
            if self.modelType == "SequentialModel" :
                self.model, self.callbacks = SequentialModel(input_shape = (len(self.SubSetVars["Jet"]),), output_shape=len(self.SamplesNames), params=self.params, modelPath=self.modelPath)
            if self.modelType == "ConvModel" :
                self.model, self.callbacks = JetImageModel(input_shape = (len(self.images_dict),self.n_eta,self.n_phi), output_shape=len(self.SamplesNames), params=self.params,modelPath=self.modelPath)
            if self.modelType == "DCLModel" :
                self.model, self.callbacks = DCLModel(input_shape = {"TopJet":(len(self.SubSetVars["Jet"]),), "PFParticles":(self.n_PF,len(self.SubSetVars["PF"])), "Image":(len(self.images_dict),self.n_eta,self.n_phi)}, output_shape=len(self.SamplesNames), params=self.params,modelPath=self.modelPath)
            self.LoadAllData()
    @timeit
    def LoadAllData(self, load=False, step=50):
        for mode in ["train", "validation", "test"]:
            self.LoadData(mode=mode)
        if load and self.isFitGenerator:
            self.inputs_train, self.labels_train = self.training_gen.GetAll(step=step)
            self.inputs_val, self.labels_val = self.validation_gen.GetAll(step=step)
            self.inputs_test, self.labels_test = self.test_gen.GetAll(step=step)
    @timeit
    def LoadData(self, mode ="train"):
        if self.isFitGenerator:
            self.SetArray("generator",mode, DataGenerator(TB= self,  mode=mode))
            # print self.GetArray("generator",mode)
        else:
            X = []
            Y = []
            for sample in self.SamplesNames:
                for ind, ID in enumerate(self.datasets[mode][sample]):
                    if "PF" in self.Object:
                        raise Exception("Selection for PFParticles is not implemented. you need to use 1 more column")
                    data = np.load(ID)[:,self.inputcolumns].astype(np.float32)
                    TopJet = np.load(ID.replace(self.Object,"TopJet")).astype(np.float32) if not "Jet" in self.Object else data
                    PT_mask = (TopJet[:,self.SubSetVars["Jet"]["jetpt"][1]]>self.pt_cut_min)*(TopJet[:,self.SubSetVars["Jet"]["jetpt"][1]]<self.pt_cut_max)
                    if len(data)==0: continue
                    if len(data.shape)==2:
                        axis_ = (1,)
                    elif len(data.shape)==4:
                        axis_ = (1,2,3)
                    mask = np.isinf(data).any(axis=axis_)
                    mask += np.isnan(data).any(axis=axis_)
                    mask += np.all(data==0,axis=axis_)
                    data = data[(~mask)*PT_mask]
                    if "Jet" in self.Object:
                        self.ApplyNormalization1(data) #TODO why 1?
                    X.append(data)
                    label = np.zeros((len(self.SamplesNames)),dtype=np.int)
                    label[self.SamplesNames[sample]]=1
                    label = [label]*len(data)
                    Y.append(label)
            X = np.concatenate(X)
            Y = np.concatenate(Y)
            self.SetArray("inputs",mode,X)
            self.SetArray("labels",mode,Y)
            # print self.GetArray("inputs",mode).shape,self. GetArray("inputs",mode).shape
    @timeit
    def FitModel(self):
        print "Training"
        a= os.system("cp "+self.FileStorageOutput+"NormInfo.txt "+ self.modelPath+"NormInfo.txt")
        if self.isFitGenerator:
            self.model.fit_generator(generator=self.training_gen, validation_data=self.validation_gen, max_queue_size= 10, use_multiprocessing=True, workers=10, epochs=self.params["epochs"], verbose=1, callbacks=self.callbacks)
            # self.model.fit_generator(generator=self.training_gen, validation_data=self.validation_gen, epochs=self.params["epochs"], verbose=1, callbacks=self.callbacks)
        else:
            # self.model.fit(self.inputs_train, self.labels_train, sample_weight=self.weights_train, batch_size=self.params["batch_size"], epochs=self.params["epochs"], verbose=1, validation_data=(self.inputs_val, self.labels_val), callbacks=self.callbacks)
            self.model.fit(self.inputs_train, self.labels_train, batch_size=self.params["batch_size"], epochs=self.params["epochs"], verbose=1, validation_data=(self.inputs_val, self.labels_val), callbacks=self.callbacks)
    @timeit
    def Predict(self, step=500):
        if self.isFitGenerator:
            self.inputs_train, self.labels_train = self.training_gen.GetAll(step=step)
            self.inputs_val, self.labels_val = self.validation_gen.GetAll(step=step)
            self.inputs_test, self.labels_test = self.test_gen.GetAll(step=step)
        self.predictions_train = self.model.predict(self.inputs_train)
        self.predictions_val = self.model.predict(self.inputs_val)
        self.predictions_test = self.model.predict(self.inputs_test)
    def SetArray(self,obj,mode,arr):
        if obj.lower()=="generator".lower():
            if mode.lower()=="train".lower():       self.training_gen = arr
            if mode.lower()=="validation".lower():  self.validation_gen = arr
            if mode.lower()=="test".lower():        self.test_gen = arr
        if obj.lower()=="inputs".lower():
            if mode.lower()=="train".lower():       self.inputs_train = arr
            if mode.lower()=="validation".lower():  self.inputs_val = arr
            if mode.lower()=="test".lower():        self.inputs_test = arr
        if obj.lower()=="labels".lower():
            if mode.lower()=="train".lower():       self.labels_train = arr
            if mode.lower()=="validation".lower():  self.labels_val = arr
            if mode.lower()=="test".lower():        self.labels_test = arr
        if obj.lower()=="pred".lower():
            if mode.lower()=="train".lower():       self.predictions_train = arr
            if mode.lower()=="validation".lower():  self.predictions_val = arr
            if mode.lower()=="test".lower():        self.predictions_test = arr
    def GetArray(self,obj,mode):
        if obj.lower()=="generator".lower():
            if mode.lower()=="train".lower():       return self.training_gen
            if mode.lower()=="validation".lower():  return self.validation_gen
            if mode.lower()=="test".lower():        return self.test_gen
        if obj.lower()=="inputs".lower():
            if mode.lower()=="train".lower():       return self.inputs_train
            if mode.lower()=="validation".lower():  return self.inputs_val
            if mode.lower()=="test".lower():        return self.inputs_test
        if obj.lower()=="labels".lower():
            if mode.lower()=="train".lower():       return self.labels_train
            if mode.lower()=="validation".lower():  return self.labels_val
            if mode.lower()=="test".lower():        return self.labels_test
        if obj.lower()=="pred".lower():
            if mode.lower()=="train".lower():       return self.predictions_train
            if mode.lower()=="validation".lower():  return self.predictions_val
            if mode.lower()=="test".lower():        return self.predictions_test
    @timeit
    def Plots(self, extraName=""):
        if self.doTraining:
            plot_losses(self.callbacks[0], min_epoch=0,  losses="loss", name=self.modelPath+"loss"+extraName)
            plot_losses(self.callbacks[0], min_epoch=0,  losses="acc",  name=self.modelPath+"acc"+extraName)
            plot_losses(self.callbacks[0], min_epoch=10, losses="loss", name=self.modelPath+"loss_10"+extraName)
            plot_losses(self.callbacks[0], min_epoch=10, losses="acc",  name=self.modelPath+"acc_10"+extraName)
        PlotRoc1vsAll(self, name = "ROC"+extraName)
        PlotRoc1vs1(self, name = "ROC"+extraName)
        PlotOutputs1D(self, name = "Outputs"+extraName)
    @timeit
    def SaveModel(self):
            with open(self.modelPath+"mymodeljson.json", "w") as f:
                f.write(str(self.model.to_json()) + "\n")
            self.model.save_weights(self.modelPath+"myweights.h5")
            self.model.save(self.modelPath+"mymodel.h5")
            ExportModel(self.model, self.modelPath+"mymodel.txt")
    @timeit
    def LoadModel(self,modelName="mymodel"):
        self.modelPath = str(self.modelPath)
        self.model = load_model(self.modelPath+modelName+".h5")
        self.model.summary()
