import pandas as pd
import numpy as np
from sequence import *
  
from sklearn.preprocessing import LabelEncoder
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestClassifier 
from sklearn.metrics import log_loss
from nolearn.lasagne import NeuralNet
from lasagne.layers import DenseLayer
from lasagne.layers import InputLayer
from lasagne.layers import DropoutLayer
from lasagne.updates import adagrad
from lasagne.nonlinearities import softmax

def LasagneNN1(num_classes, num_features, n_units0, p_dropout0, learning_rate, n_epochs):    
    layers = [ ('input', InputLayer),
               ('dense0', DenseLayer),('dropout0', DropoutLayer),
               ('output', DenseLayer)]
    net1 = NeuralNet(layers=layers,
                 input_shape=(None, num_features),
                 dense0_num_units=n_units0, dropout0_p=p_dropout0,
                 output_num_units=num_classes,
                 output_nonlinearity=softmax,
                 update=adagrad,
                 update_learning_rate=learning_rate,
                 eval_size=0.2,
                 verbose=1,
                 max_epochs=n_epochs)
    return net1


class RGENPrediction:
    def __init__(self, model):
        self.model=model
    def __enter__(self):
        return self
    def __exit__(self, *err):
        pass
    def __del__(self):
        #print "deling", self
        pass

    def Train(self, data):
        index=list(data.index)
        np.random.shuffle(index)
        shuffled_data = data.ix[index]
        
        X=shuffled_data.drop('y',axis=1).values
        encoder=LabelEncoder()
        y = encoder.fit_transform(shuffled_data.y.values).astype(np.int32)
        
        self.model.fit(X,y)    
        
    def TrainWithBaggingAndTest(self, datalist, n_iter, test_1, test_0):
        for i in range(n_iter):
            idx_data = np.random.choice(len(datalist))
            data=datalist[idx_data]
            self.Train(data)

            status=self.model.train_history_[-1:][0]
            if status['valid_loss'] < 0.00001:
                break
            if i%10==1:
                print '%d-th train, %d-th dataset, train=%.5r valid=%.5r ratio=%.5r accuracy=%.5r'%(\
                i, idx_data, status['train_loss'], status['valid_loss'],\
                status['train_loss']/status['valid_loss'] ,status['valid_accuracy'])

        self.UpdatePrediction(test_1, test_0)