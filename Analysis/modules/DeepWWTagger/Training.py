from PostTrainingPlots import *
# from TaggerBase import *

if len(sys.argv) > 1:
    fileName = sys.argv[1]
else:
    fileName = "trainings/model.json"

with open(fileName) as f:
    ModelInfoDict = json.load(f)

print "Training on:"
print fileName
prettydic(ModelInfoDict)


TB=TaggerBase(ModelInfoDict=ModelInfoDict, Collection=ModelInfoDict["Collection"])
TB.MakeModel()
TB.FitModel()
TB.Predict()
TB.Plots()
TB.SaveModel()

# print "Start PostTrainingEvaluation"


# PTE = PostTrainingEvaluation(ModelInfoDict=ModelInfoDict)
# print "MakePredictions"
# PTE.MakePredictions()
# print "PlotsPostTraining"
# PTE.PlotsPostTraining()
