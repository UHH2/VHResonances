import os
from keras.models import Sequential, Model
from keras.layers import Input, Dense, Dropout, BatchNormalization, Flatten, Conv2D, MaxPooling2D, LSTM, merge
from keras.optimizers import Adam
from keras import metrics, regularizers
from keras.callbacks import History, ModelCheckpoint, ReduceLROnPlateau, Callback


import random
import numpy as np
import keras

from PlottingFunctions import *

class DataGenerator(keras.utils.Sequence):
    'Generates data for Keras'
    def __init__(self, TB, mode="", shuffle=True):
        'Initialization'
        self.TB = TB
        self.shuffle = shuffle
        self.mode = mode
        self.batch_size = self.TB.params["batch_size"]
        self.labels = TB.SamplesNames
        self.datasets = TB.datasets[self.mode]
        self.n_classes = len(TB.SamplesNames)
        self.input_shape = TB.inputcolumns
        self.SampleRef = self.TB.SampleRef
        self.step = int(self.batch_size/self.n_classes)
        if self.step<1:
            raise Exception("batch_size<n_classes! Not supported!")
        if any(len(self.datasets[sample]) != len(self.datasets[self.SampleRef]) for sample in self.labels):
            raise Exception("Sample of different length! Not supported!", str([(x, len(self.datasets[x])) for x in self.labels]))
        self.on_epoch_end()
        self.debug = False
        # self.debug = True


    def __len__(self):
        'Denotes the number of batches per epoch. Used internally by Keras'
        #TODO check if tot number of files is respected (keras should add 1 internally?)
        if self.debug:
            print "@ __len__"+" "+str(len(self.datasets)*self.TB.filesize)+" "+str(self.batch_size)+" "+str(int(np.floor((len(self.datasets)*self.TB.filesize) / self.batch_size)))+" "+str(int(np.floor((len(self.datasets)))))
        self.on_epoch_end()
        return len(self.indexes)

    def __getitem__(self, index):
        if self.debug: print "@ __getitem__"+" "+str("start ")+str(self.input_shape)
        index = self.indexes[index]["index"]
        batch_min = self.indexes[index]["batch_min"]
        batch_max = self.indexes[index]["batch_max"]
        if self.debug: print "START LOOP"+str(index)+" "+str(batch_min)+" "+str(batch_max)
        if type(self.TB.Object)==type(list()):
            X = {}
            for obj in self.TB.Object:
                X[obj] = []
        else:
            X = []
        Y = []
        for sample in self.labels:
            if "PF" in self.TB.Object:
                raise Exception("Selection for PFParticles is not implemented. you need to use 1 more column")
            ID = self.datasets[sample][index]
            if self.debug: print ID
            if type(self.TB.Object)==type(list()):
                data = {}
                for obj in self.TB.Object:
                    if "PF" in obj:
                        data[obj] = np.load(ID.replace("TopJet",obj))[batch_min:batch_max,:,self.input_shape[obj]].astype(np.float32)
                        #TODO Select here NParticles wanted
                    else:
                        data[obj] = np.load(ID.replace("TopJet",obj))[batch_min:batch_max,self.input_shape[obj]].astype(np.float32)
                TopJet = np.load(ID.replace(self.TB.Object[0],"TopJet"))[batch_min:batch_max].astype(np.float32)
                Images = np.load(ID.replace(self.TB.Object[0],"Image"))[batch_min:batch_max].astype(np.float32)
            else:
                data = np.load(ID)[batch_min:batch_max,self.input_shape].astype(np.float32)
                TopJet = np.load(ID.replace(self.TB.Object,"TopJet"))[batch_min:batch_max].astype(np.float32)
                Images = np.load(ID.replace(self.TB.Object,"Image"))[batch_min:batch_max].astype(np.float32)
            PT_mask = (TopJet[:,self.TB.SubSetVars["Jet"]["jetpt"][1]]>self.TB.pt_cut_min)*(TopJet[:,self.TB.SubSetVars["Jet"]["jetpt"][1]]<self.TB.pt_cut_max)
            if self.debug: print "check1", np.load(ID).shape, (data.shape if type(self.TB.Object)!=type(list()) else [data[x].shape for x in data]), batch_min, batch_max, ID
            if self.debug: print "shape check", Images[PT_mask].shape, TopJet[PT_mask].shape, " exp",
            mask = np.isinf(TopJet).any(axis=(1,))
            mask += np.isinf(Images).any(axis=(1,2,3))
            mask += np.isnan(TopJet).any(axis=(1,))
            mask += np.isnan(Images).any(axis=(1,2,3))
            mask += np.all(TopJet==0,axis=(1,))
            mask += np.all(Images==0,axis=(1,2,3))
            if self.debug: print Images[(~mask)*PT_mask].shape, TopJet[(~mask)*PT_mask].shape
            if type(self.TB.Object)==type(list()):
                for obj in self.TB.Object:
                    data[obj] = data[obj][(~mask)*PT_mask]
            else:
                data = data[(~mask)*PT_mask]
            if self.debug:
                if type(self.TB.Object)==type(list()):
                    print "post", str([data[x].shape for x in data])
                else:
                    print "post", data.shape
            if type(self.TB.Object)==type(list()):
                for obj in self.TB.Object:
                    if len(data[obj])==0: continue
            else:
                if len(data)==0: continue
            if "Jet" in self.TB.Object:
                self.TB.ApplyNormalization1(data)
            if type(self.TB.Object)==type(list()):
                for obj in self.TB.Object:
                    if "Jet" in obj: self.TB.ApplyNormalization1(data[obj])
                    X[obj].append(data[obj])
            else:
                X.append(data)
            label = np.zeros((self.n_classes),dtype=np.int)
            label[self.labels[sample]]=1 #TODO PASS THE CORRECT INDEX
            if self.debug:
                print "load"+" "+str(ID[115:])+" Sample:"+" "+str(sample)+" shape "
                if type(self.TB.Object)==type(list()):
                    for obj in self.TB.Object:
                        print str(data[obj].shape)+" ",
                else:
                    print str(data.shape)+" ",
                print "n_classes "+str(self.n_classes)+" label "+str(self.labels[sample])+" "+str(label)
            if type(self.TB.Object)==type(list()):
                label = [label]*len(data[self.TB.Object[0]])
            else:
                label = [label]*len(data)
            Y.append(label)
        if type(self.TB.Object)==type(list()):
            for obj in self.TB.Object:
                X[obj] = np.concatenate(X[obj])
        else:
            X = np.concatenate(X)
        Y = np.concatenate(Y)
        if self.debug:
            print "@ __getitem__"+" "+"END ",
            if type(self.TB.Object)==type(list()):
                for obj in self.TB.Object:
                    print str(X[obj].shape)+" "+str(X[obj].dtype)+" ",
            else:
                print str(X.shape)+" "+str(X.dtype),
            print str(Y.shape)+" "+str(Y.dtype)
            print "IT's OK, ready to return"
        if self.debug:
            import time
            print "\n"
            time.sleep(10)
            print "\n"
        if type(self.TB.Object)==type(list()):
            X = [X[obj] for obj in self.TB.Object]
        return X, Y


    def on_epoch_end(self):
        #TODO make sure that nevents are the same for each sample
        'Updates indexes after each epoch'
        batches = range(0,self.TB.filesize,self.step)
        if self.shuffle == True:
            np.random.shuffle(batches)
        indexes = {}
        for sample in self.labels:
            indexes[sample] = list(range(len(self.datasets[sample])))
            np.random.shuffle(indexes[sample])
        self.indexes = []
        for index in range(len(indexes[self.SampleRef])):
            for x in batches:
                self.indexes.append({"index": index, "batch_min":x, "batch_max":x+self.step})

    def GetAll(self, step=None):
        if type(self.TB.Object)==type(list()):
            X = {}
            for obj in self.TB.Object:
                X[obj] = []
        else:
            X = []
        Y = []
        for i in range( min(step,self.__len__()) if step!=None else self.__len__()):
            x,y = self.__getitem__(i)
            if type(self.TB.Object)==type(list()):
                for index,obj in enumerate(self.TB.Object):
                    X[obj].append(x[index])
            else:
                X.append(x)
            Y.append(y)
        if type(self.TB.Object)==type(list()):
            for obj in self.TB.Object:
                X[obj] = np.concatenate(X[obj])
        else:
            X = np.concatenate(X)
        Y = np.concatenate(Y)
        if type(self.TB.Object)==type(list()):
            X = [X[obj] for obj in self.TB.Object]
        return X,Y


def SequentialModel(input_shape,output_shape,params,modelPath):
    model = Sequential()
    # Define layers
    model.add(Dense(params["DenseLayer"][0], input_shape=input_shape, activation=params["activation"],kernel_initializer=params["kernel_initializer"]))
    for i in range(1,len(params["DenseLayer"])):
        model.add(Dense(params["DenseLayer"][i], activation="tanh",kernel_initializer=params["kernel_initializer"],bias_initializer=params["bias_initializer"]))
        model.add(BatchNormalization())
        model.add(Dropout(params["dropoutRate"]))
    model.add(Dense(output_shape, activation=params["activation_last"]))
    # Compile
    myloss = "categorical_crossentropy"
    if output_shape == 1:
        myloss = "binary_crossentropy"
    # TODO "tanh" "Adamax"
    model.compile(loss=myloss, optimizer="adamax", metrics=params["metrics"])
    model.summary()
    callbacks = DefineCallbacks(modelPath)
    # TODO
    # save model params
    # plot_model(model, to_file="model.png", show_shapes=True, show_layer_names=False, rankdir="LR")
    return (model, callbacks)

def JetImageModel(input_shape,output_shape,params,modelPath):
    model = Sequential()
    # Define ConvLayer
    model.add(Conv2D(params["ConvLayer"][0], input_shape=input_shape, kernel_size=params["kernel_size"][0], strides=params["strides"][0], padding=params["padding"], activation=params["activation"], data_format=params["data_format"]))
    model.add(MaxPooling2D())
    model.add(BatchNormalization())
    for i in range(1,len(params["ConvLayer"])):
        model.add(Conv2D(params["ConvLayer"][i], kernel_size=params["kernel_size"][i], strides=params["strides"][i], padding=params["padding"], activation=params["activation"], data_format=params["data_format"]))
        if i<1: # TODO understand algebra
            model.add(MaxPooling2D())
        model.add(BatchNormalization())
        # model.add(Dropout(params["dropoutRate"]))
    model.add(Flatten())
    for i in range(0,len(params["DenseLayer"])):
        model.add(Dense(params["DenseLayer"][i], activation=params["activation"],kernel_initializer=params["kernel_initializer"],bias_initializer=params["bias_initializer"]))
        model.add(BatchNormalization())
        # model.add(Dropout(params["dropoutRate"]))
    model.add(Dense(output_shape, activation=params["activation_last"]))
    # Compile
    myloss = "categorical_crossentropy"
    if output_shape == 1:
        myloss = "binary_crossentropy"
    model.compile(loss=myloss, optimizer=params["optimizer"], metrics=params["metrics"])
    model.summary()
    callbacks = DefineCallbacks(modelPath)
    return (model, callbacks)


def DCLModel(input_shape,output_shape,params,modelPath):
    # Define Inputs
    TJs= Input(shape=input_shape["TopJet"])
    PFs = Input(shape=input_shape["PFParticles"])
    Ims = Input(shape=input_shape["Image"])

    # Define Recursive
    rec = LSTM(params["LSTMLayer"])(PFs)
    rec = Dropout(params["dropoutRate"])(rec)

    # Define ConvLayer
    conv = Conv2D(params["ConvLayer"][0], kernel_size=params["kernel_size"][0], strides=params["strides"][0], padding=params["padding"], activation=params["activation"], data_format=params["data_format"])(Ims)
    conv = MaxPooling2D()(conv)
    conv = BatchNormalization()(conv)
    conv = Dropout(params["dropoutRate"])(conv)
    for i in range(1,len(params["ConvLayer"])):
        conv = Conv2D(params["ConvLayer"][i], kernel_size=params["kernel_size"][i], strides=params["strides"][i], padding=params["padding"], activation=params["activation"], data_format=params["data_format"])(conv)
        if i<1: # TODO understand algebra
            conv = MaxPooling2D()(conv)
        conv = BatchNormalization()(conv)
        conv = Dropout(params["dropoutRate"])(conv)

    conv = Flatten()(conv)

    # Combine NN
    preds = merge( [TJs, rec, conv] , mode='concat')

    for i in range(0,len(params["DenseLayer"])):
        preds = Dense(params["DenseLayer"][i], activation=params["activation"],kernel_initializer=params["kernel_initializer"],bias_initializer=params["bias_initializer"])(preds)
        preds = BatchNormalization()(preds)
        preds = Dropout(params["dropoutRate"])(preds)

    preds = Dense(output_shape, activation=params["activation_last"])(preds)

    model = Model(inputs=[TJs, PFs, Ims], outputs=preds)
    # Compile
    myloss = "categorical_crossentropy"
    if output_shape == 1:
        myloss = "binary_crossentropy"

    model.compile(loss=myloss, optimizer=params["optimizer"], metrics=params["metrics"])
    model.summary()
    callbacks = DefineCallbacks(modelPath)
    return (model, callbacks)

class Plot(Callback):
    def __init__(self, modelPath):
        self.modelPath = modelPath
        self.history = {}
    def on_epoch_end(self, epoch, logs=None):
        for k, v in logs.items():
            self.history.setdefault(k, []).append(v)
        plot_losses(self, losses="loss", name=self.modelPath+"loss_training")
        plot_losses(self, losses="acc",  name=self.modelPath+"acc_training")


def DefineCallbacks(modelpath):
    if not os.path.exists(modelpath+"/models/"):
        os.makedirs(modelpath+"/models/")
    callbacks = []
    history = History()
    callbacks.append(history)
    modelCheckpoint_loss_best   = ModelCheckpoint(modelpath+"models/bestmodel_epoch{epoch:03d}_loss{val_loss:.2f}.h5", monitor="val_loss", save_best_only=True)
    callbacks.append(modelCheckpoint_loss_best)
    modelCheckpoint_acc_best    = ModelCheckpoint(modelpath+"models/bestmodel_epoch{epoch:03d}_acc{val_acc:.2f}.h5",   monitor="val_acc",  save_best_only=True)
    callbacks.append(modelCheckpoint_acc_best)
    reduceLROnPlateau           = ReduceLROnPlateau(monitor="val_loss", factor=0.5, patience=1, min_lr=0.001, cooldown=10)
    callbacks.append(reduceLROnPlateau)
    PlotTraining = Plot(modelPath=modelpath)
    callbacks.append(PlotTraining)
    return callbacks


# class DataGenerator_Save(keras.utils.Sequence):
#     'Generates data for Keras'
#     def __init__(self, TB, datasets, labels, input_shape=None, filesize = 10000, batch_size=64, n_classes=10, shuffle=True, mode=""):
#         'Initialization'
#         self.TB = TB
#         self.TB.filesize = filesize
#         self.batch_size = batch_size
#         self.labels = labels
#         self.datasets = datasets
#         self.n_classes = n_classes
#         self.shuffle = shuffle
#         self.input_shape = input_shape
#         self.mode = mode
#         self.on_epoch_end()
#         self.debug = False
#
#     def __len__(self):
#         'Denotes the number of batches per epoch. Used internally by Keras'
#         #TODO check if tot number of files is respected (keras should add 1 internally?)
#         if self.debug: print "@ __len__"+" "+str(len(self.datasets)*self.TB.filesize)+" "+str(self.batch_size)+" "+str(int(np.floor((len(self.datasets)*self.TB.filesize) / self.batch_size)))+" "+str(int(np.floor((len(self.datasets)))))
#         return int(np.floor((len(self.datasets))))
#
#     def __getitem__(self, index):
#         if self.debug: print "@ __getitem__"+" "+"start"+" "+str(index)
#         'Generate one batch of data'
#         # Generate indexes of the batch
#         # indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]
#         indexes = self.indexes[index:(index+1)]
#
#         if self.debug:
#             text = "@ __getitem__"+" "+str("indexes")+" "+str("[")
#             for x in indexes:
#                 text += str(x)+","
#             text += str("]")
#             print text
#
#         # Find list of IDs
#         datasets_temp = [self.datasets[k] for k in indexes]
#         if self.debug: print "@ __getitem__"+" "+str("datasets_temp")+" "+str(len(datasets_temp))
#
#         # Generate data
#         X, Y = self.__data_generation(datasets_temp)
#         if self.debug: print "@ __getitem__"+" "+str(X.shape)+" "+str(Y.shape)+" "+str(X.dtype)+" "+str(Y.dtype)
#         # print "@ __getitem__"+" "+str(X.shape)+" "+str(Y.shape)+" "+str(X.dtype)+" "+str(Y.dtype)
#
#         return X, Y
#
#     def on_epoch_end(self):
#         'Updates indexes after each epoch'
#         self.indexes = np.arange(len(self.datasets))
#         if self.shuffle == True:
#             np.random.shuffle(self.indexes)
#
#     def __data_generation(self, datasets_temp):
#         if self.debug: print "@ __data_generation"+" "+str("start ")+str(self.input_shape)
#         # X = np.zeros((len(datasets_temp),self.input_shape[0]))
#         # X[ind] = np.load(ID)[0,:self.input_shape[0]]
#         # Y = np.zeros((len(datasets_temp),self.n_classes),dtype=np.int)
#         X = []
#         Y = []
#         if self.debug: print "START LOOP", len(datasets_temp)
#         for ind, ID in enumerate(datasets_temp):
#             for x in self.labels:
#                 if x in ID : Sample = x
#             data = np.load(ID)[:,self.input_shape].astype(np.float32)
#             if self.debug: print "check1", np.load(ID).shape, data.shape, ID
#             if len(data.shape)==2:
#                 axis_ = (1,)
#             elif len(data.shape)==4:
#                 axis_ = (1,2,3)
#             mask = np.isinf(data).any(axis=axis_)
#             mask += np.isnan(data).any(axis=axis_)
#             data = data[~mask]
#             self.TB.ApplyNormalization1(data)
#             X.append(data)
#             label = np.zeros((self.n_classes),dtype=np.int)
#             label[self.labels[Sample]]=1 #TODO PASS THE CORRECT INDEX
#             if self.debug: print "load"+" "+str(ID[115:])+" Sample:"+" "+str(Sample)+" shape "+str(data.shape)+" n_classes "+str(self.n_classes)+" label "+str(self.labels[Sample])+" "+str(label)
#             # print "\nload "+str(self.mode)+" "+str(ID[115:])+" Sample:"+" "+str(Sample)+" shape "+str(data.shape)+" n_classes "+str(self.n_classes)+" label "+str(self.labels[Sample])+" "+str(label)
#             label = [label]*len(data)
#             with open("test_iter_"+self.mode+".txt", "a") as f:
#                 f.write(ID+"\t"+str(data[0,0])+"\t"+str(label[0])+"\n")
#             # Y[ind:self.labels[Sample]] = 1
#             Y.append(label)
#         X = np.concatenate(X)
#         Y = np.concatenate(Y)
#         if self.debug: print "@ __data_generation"+" "+"END"+" "+str(X.shape)+" "+str(Y.shape)+" "+str(X.dtype)+" "+str(Y.dtype)
#         # import time
#         # print "\n"
#         # time.sleep(2)
#         # print "\n"
#         return X, Y
