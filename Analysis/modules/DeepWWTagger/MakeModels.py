import os, argparse
import json

def SetDefaultModel():
    ModelInfoDict = {"doTraining"   : True,
                     "isSubset"     : False,
                     "SubSetVars"   : ["jetpt", "jeteta", "jetphi", "jetenergy", "jettau1", "jettau2", "jettau3", "jettau4", "PFpt", "PFeta", "PFphi", "PFenergy", "PFcharge","PFparticleID", "PFpuppi_weight"],
                     "pt_min"       : 300,
                     "pt_max"       : 500,
                     "max_size"     : -1,
                     "SamplesNames" : ["MC_HWW", "MC_Top", "MC_QCD", "MC_W", "MC_Z"],
                     "modelName"    : "model",
                     "modelType"    : "ConvModel",
                     "ConvLayer"    : [48,24,12,2],
                     "kernel_size"  : [(3,3),(3,3),(3,3),(3,3)],
                     "strides"      : [(1,1),(1,1),(1,1),(1,1)],
                     "DenseLayer"   : [50,50,10],
                     "LSTMLayer"    : 64,
                     "params"       : {"epochs" : 500, "batch_size" : 512, "activation" : "relu", "kernel_initializer": "glorot_normal", "bias_initializer": "ones", "activation_last": "softmax", "padding" : "valid", "data_format" : "channels_first", "optimizer": "adam", "metrics":["accuracy"], "dropoutRate": 0.01}
                     }
    return ModelInfoDict



def WriteModeljson(ModelInfoDict):
    path = "submitFiles/Trainings/"+ModelInfoDict["Collection"]+"/"+ModelInfoDict["modelType"]
    a = os.system("mkdir -p "+path)
    # print path+"/"+ModelInfoDict["modelName"]
    with open(path+"/"+ModelInfoDict["modelName"]+".json","w") as f:
        f.write(json.dumps(ModelInfoDict))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--Collection', action='store', dest="Collection",default="Puppi")
    args = parser.parse_args()

    DenseLayers = [[50,50,10], [100,100,10],[100,100,100,100,30,10],[50,50,50,50,30,10]]
    ConvLayers  = [[48,24,8],[128,64,24,8],[48,128,256]]
    LSTMLayers  = [32,64,128,256,512]

    # pt_mins=[300,500,300]
    # pt_maxs=[500,10000,10000]
    pt_mins     = [200]
    pt_maxs     = [4000]
    NEpochs     = [20,300]
    NEvents     = [300,1000, "all"]
    NSamples    = ["reduced"] #["reduced","full"]
    modelTypes  = ["DCLModel","ConvModel","SequentialModel"]

    for nSamples in NSamples:
        for modelType in modelTypes:
            for index_pt,pt_min in enumerate(pt_mins):
                pt_max=pt_maxs[index_pt]
                for nEpochs in NEpochs:
                    for nEvents in NEvents:
                        for DenseLayer in DenseLayers:
                            for ConvLayer in ConvLayers if modelType!="SequentialModel" else [[0]]:
                                for LSTMLayer in LSTMLayers if modelType=="DCLModel" else [[0]]:
                                    extraText = ""
                                    if modelType!="SequentialModel":
                                        for x in ConvLayer:
                                            extraText += "_"+str(x)
                                        extraText += "_Conv"
                                        if modelType!="ConvModel":
                                            extraText += "_"+str(LSTMLayer)+"_LSTM"
                                    for x in DenseLayer:
                                        extraText += "_"+str(x)
                                    extraText += "_Dense"
                                    extraText += "_pt"+str(pt_min)+"_"+str(pt_max)
                                    extraText += "_"+nSamples
                                    ModelInfoDict = SetDefaultModel()
                                    ModelInfoDict["max_size"]   = nEvents*1000 if nEvents!="all" else -1
                                    ModelInfoDict["isFitGenerator"] = True
                                    ModelInfoDict["params"]["epochs"] = nEpochs
                                    ModelInfoDict["modelType"]  = modelType
                                    ModelInfoDict["modelName"]  = "model_"+str(nEvents)+"kevents_"+str(nEpochs)+"epochs"+extraText
                                    ModelInfoDict["DenseLayer"] = DenseLayer
                                    ModelInfoDict["ConvLayer"]  = ConvLayer
                                    ModelInfoDict["LSTMLayer"]  = LSTMLayer
                                    ModelInfoDict["pt_min"]     = pt_min
                                    ModelInfoDict["pt_max"]     = pt_max
                                    ModelInfoDict["SamplesNames"] = ["MC_HWW", "MC_QCD"] if "reduced"==nSamples else ["MC_HWW", "MC_Top", "MC_QCD", "MC_W", "MC_Z"]
                                    ModelInfoDict["Collection"] = args.Collection
                                    WriteModeljson(ModelInfoDict)
