import os
import json

path = "/nfs/dust/cms/user/amalara/WorkingArea/File/NeuralNetwork/model_6SampleSelected_reduced_300epochs_1000kevents/"

file_ = "mymodelinfo"

with open(path+file_+".json") as f:
    ModelInfoDict = json.load(f)

with open(path+file_+".txt","w") as f:
    for el in ModelInfoDict:
        f.write(str(el)+"\t"+str(ModelInfoDict[el]).replace("u'","").replace("'","")+"\n")


a = os.system("cp /nfs/dust/cms/user/amalara/WorkingArea/File/NeuralNetwork/input_varariables/VHResonances/Puppi/leptonchannel/NormInfo.txt " + path)
