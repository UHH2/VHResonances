import numpy as np


def WriteArray(fout, array):
    array = array.tolist()
    str_ = ""
    for el in array:
        str_ += " " if el>0 else ""
        str_ += str(round(el,6)) + ",\t"
    fout.write(str_+"\n")


def WriteWeights(fout, weights):
    if len(weights.shape) == 1:
        WriteArray(fout, weights)
    elif len(weights.shape) == 2:
        for w in weights:
            WriteArray(fout, w)

def ExportModel(model, filename):
    with open(filename, 'w') as fout:
        for ind, l in enumerate(model.get_config()):
            if l["class_name"] == "Dropout":
                continue
            fout.write("New Layer\n")
            fout.write(l["class_name"]+"\n")
            if l['class_name'] == "Dense":
                fout.write("weights"+"\n")
                WriteWeights(fout, model.layers[ind].get_weights()[0])
                fout.write("bias"+"\n")
                WriteWeights(fout, np.expand_dims(model.layers[ind].get_weights()[1],axis=0) )
                fout.write("activation\t"+l["config"]["activation"]+"\n")
            if l['class_name'] == "BatchNormalization":
                fout.write("epsilon"+"\n")
                fout.write(str(model.layers[ind].epsilon)+"\n")
                fout.write("gamma"+"\n")
                WriteWeights(fout, np.expand_dims(model.layers[ind].get_weights()[0], axis=0) )
                fout.write("beta"+"\n")
                WriteWeights(fout, np.expand_dims(model.layers[ind].get_weights()[1], axis=0) )
                fout.write("moving_mean"+"\n")
                WriteWeights(fout, np.expand_dims(model.layers[ind].get_weights()[2], axis=0) )
                fout.write("moving_variance"+"\n")
                WriteWeights(fout, np.expand_dims(model.layers[ind].get_weights()[3], axis=0) )
