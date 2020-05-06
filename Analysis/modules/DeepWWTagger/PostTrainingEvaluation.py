from PostTrainingPlots import *
# from TaggerBase import *

fileName = {}
# fileName["Dense"] = "model_300kevents_10epochs_100_100_100_100_30_10_Dense_pt300_10000_reduced.json"
# fileName["Conv"] = "model_1000kevents_10epochs_48_24_8_Conv_50_50_10_Dense_pt300_10000_reduced.json"
# fileName["DCLModel"] = "DCLmodel.json"
fileName["DCLModel"] = "/nfs/dust/cms/user/amalara//WorkingArea/File//NeuralNetwork/trainings/HiggsToWWTagger/2017/HOTVR/leptonchannel/DCLModel/model_300kevents_20epochs_48_24_8_Conv_50_50_10_Dense_pt200_4000_reduced/mymodelinfo.json"
fileName["DCLModel"] = "/nfs/dust/cms/user/amalara//WorkingArea/File//NeuralNetwork/trainings/HiggsToWWTagger/2017/Puppi/leptonchannel/DCLModel/model_300kevents_300epochs_48_128_256_Conv_256_LSTM_50_50_10_Dense_pt200_4000_reduced/mymodelinfo.json"

PTE = {}

for model in fileName.keys():
    with open(fileName[model]) as f:
        ModelInfoDict = json.load(f)

    print "PostTrainingEvaluation on:"
    print fileName[model]
    prettydic(ModelInfoDict)
    PTE[model] = PostTrainingEvaluation(Collection=ModelInfoDict["Collection"] , ModelInfoDict=ModelInfoDict)
    PTE[model].LoadAllData(load=True, step=2)
    print PTE[model].inputs_train[0].shape, PTE[model].inputs_train[1].shape, PTE[model].inputs_train[2].shape
    PTE[model].MakePredictions()
    PTE[model].PlotsPostTraining()

# PTE = PostTrainingEvaluation(ModelInfoDict=ModelInfoDict)
# PTE.MakePredictions()
# PTE_JET = PostTrainingEvaluation(ModelInfoDict=ModelInfoDict, Object="TopJet")
# PTE_JET.LoadAllData(load=True)
#
# PTE_JET.GetArray("inputs","train").shape

# print "PlotsPostTraining"
# PTE.PlotsPostTraining()
