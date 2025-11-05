# ============================================================
# Author: Dr. Mori (Morteza) Chalabi
# Created: October 2025
# Description: Deep learning model trained on PDAC TME (LIV-MR signature)
# ============================================================

# Importing libraries

import fnmatch                                  # lisiting files with specific patterns as in list.files in R
import os
import numpy as np
import scipy
import pandas as pd

# R-related imports
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects.packages import importr      # for being able to import R pacakges

# loading R script

r.source("preprocess.R")                      # loads R script
r_preprocess = ro.globalenv['preprocess']     # handle to preprocess() function

# Preprocessing method
def load_examples(min_features = 200, max_features = float('inf'), max_percent_mt = 10):

  # reading in and preprocessing each file

  examples_ = fnmatch.filter(os.listdir('data/train/'), '*.mtx')                      # listing example files
  examples_ = [s_.split('_')[0] for s_ in examples_]                                  # extracting example IDs/names
  labels_ = pd.read_table('data/train/labels.tsv', header = None, index_col = 0)      # reading in labels file
  exmp_vects = []                                                                     # list to hold examples' vectors: n by m array where m is number of examples and n is number of features
  exmp_labs = []                                                                      # list to hold examples' labels: 1 by m array

  for s_ in examples_:

    # reading in matrix-related files

    mat_fl = str('data/train/'+s_+'_matrix.mtx')
    features_fl = str('data/train/'+s_+'_features.tsv')
    barcodes_fl = str('data/train/'+s_+'_barcodes.tsv')

    # preprocessing count matrix

    exmp_vect = r_preprocess(mat_ = mat_fl, features_ = features_fl, barcodes_ = barcodes_fl, sample_ = s_, min_features = min_features, max_features = max_features, max_percent_mt = max_percent_mt)

    # updating lists

    exmp_vect = np.array(exmp_vect)         # convert to numpy array (must be 1 x n)
    exmp_vects.append(exmp_vect)            # updating train list

    exmp_labs.append(labels_.loc[s_,1])     # updating train labels list

  # stack into (m, n) array

  exmp_vects = np.hstack(exmp_vects)                              # hstack to get n by m array to comply with deep learning libraries
  exmp_labs = np.array(exmp_labs).reshape(1, len(examples_))      # conversion to 1 by m numpy array

  # saving as numpy array

  np.savez('data/train/examples_labels.npz', exmp_vects = exmp_vects, exmp_labs = exmp_labs)
  
  return exmp_vects, exmp_labs

# Train method
# The input shapes must be:
# X: n x m numpy array (m is number of samples)
# Y: 1 x m numpy array
# W: n x 1 numpy array (n is number of features)
# b = 0
def train(min_features = 200, max_features = float('inf'), max_percent_mt = 10, alpha_ = 0.5, num_iterations = 2000, print_cost = True):
  
  # load examples
  
  if not os.path.exists('data/train/examples_labels.npz'):
    X_, Y_ = load_examples(min_features = min_features, max_features = max_features, max_percent_mt = max_percent_mt)
  else:
    objects_ = np.load('data/train/examples_labels.npz')
    X_ = objects_['exmp_vects']
    Y_ = objects_['exmp_labs']
  
  # training
  
  n_ = X_.shape[0]                                                                  # number of features
  m_ = X_.shape[1]                                                                  # number of examples
  Y_ = np.array([1 if(l_ == 'liver') else 0 for l_ in Y_[0,:]]).reshape(1, m_)      # encoding labels and reshaping Y to be a 1 by m array
  W_ = np.zeros((n_, 1))                                                            # weights to be updated
  b_ = 0                                                                            # bias to be updated
  Js_ = []                                                                          # list to hold costs
  for i_ in range(num_iterations):
    
    # forward propagation
    
    Z_ = np.dot(W_.T,X_) + b_                               # linear function
    Y_hat = scipy.special.expit(Z_)                         # sigmoid activation function: 1/(1 + e^-z)
    L_ = (Y_*np.log(Y_hat) + (1-Y_)*np.log(1 - Y_hat))      # binary cross-entropy loss (aka log loss or log likelihood)
    J_ = (-1/m_) * np.sum(L_)                               # cost function
    
    # backward probagation
    
    nabla_L_Z = (Y_hat - Y_).T                    # derivative of loss w.r.t Z; shape: m x 1
    nabla_J_W = (1/m_)*np.dot(X_, nabla_L_Z)      # derivative of cost w.r.t W; shape: n x 1
    nabla_J_b = (1/m_)*np.sum(nabla_L_Z)          # derivative of cost w.r.t b; shape: scalar
    
    # updating parameters using gradient descent
    
    W_ = W_ - alpha_*nabla_J_W
    b_ = b_ - alpha_*nabla_J_b
    
    # printing cost every 100 iterations
    
    if print_cost and i_ % 100 == 0:
      Js_.append(J_)      # appending cost to costs list
      print('cost at iteration %i is: %.5f' % (i_, J_))

  # performance of model on training set

  Z_ = np.dot(W_.T,X_) + b_                     # linear function with learned W and b
  Y_hat = scipy.special.expit(Z_)               # sigmoid activation function
  Y_hat[Y_hat < 0.5] = 0
  Y_hat[0.5 <= Y_hat] = 1
  rslt_ = np.mean(Y_ == Y_hat.astype(int))      # elementwise XNOR and mean to get accuracy
  print('\n[ Accuracy on training set is: %.2f%% ]\n' % (rslt_*100))

  return W_, b_, np.array(Js_)

# Predict method
def predict(mat_ = None, features_ = None, barcodes_ = None, label_ = None, example_name = None,  min_features = 200, max_features = float('inf'), max_percent_mt = 10, W_ = None, b_ = None):
  
  x_ = r_preprocess(mat_ = mat_, features_ = features_, barcodes_ = barcodes_, sample_ = example_name, min_features = min_features, max_features = max_features, max_percent_mt = max_percent_mt)
  x_ = np.array(x_)
  
  # forward propagation
  
  y_ = label_                         # true label
  z_ = np.dot(W_.T,x_) + b_           # linear function
  y_hat = scipy.special.expit(z_)     # sigmoid activation function
  
  # converting probabilities to class labels
  
  class_ = 'liver' if(0.5 <= y_hat) else 'non-liver site'
  
  return y_hat, class_


# Main body

if __name__ == "__main__":

  W_, b_, costs_ = train(alpha_ = 0.5)      # calling train function to get learned parameters
  # prob_, class_ = predict(mat_ = 'data/train/PN13_matrix.mtx', features_ = 'data/train/PN13_features.tsv', barcodes_ = 'data/train/PN13_barcodes.tsv', label_ = 'liver', example_name = 'PN13', W_ = W_, b_ = b_)
  # print("\n[ Probability of liver recurence is %.2f%%. Therefore, the sample may recur in %s ]\n" % (prob_[0,0], class_))
