#
# mse_check.py - To check the result of libmfcc.c and compare with the output of python_speech_features 
# Author - Akarsh Prabhakara	
#

import numpy as np
from python_speech_features import mfcc

raw = [line.rstrip('\n') for line in open('recorded_data.dat')]
raw = np.array(raw)
rawnp = raw.astype(np.int)

# If you change the MFCC / audio parameters in real_time_mfcc.c, supplement the right arguments for mfcc() 
mf = mfcc(rawnp)

result = [line.rstrip('\n') for line in open('mfcc_result.dat')]
result = np.array(result)
resultnp = result.astype(np.float64)
mf1 = np.zeros((mf.shape[0],mf.shape[1]))
for i in range(mf.shape[0]):
	for j in range(mf.shape[1]):
		mf1[i,j] = resultnp[i*mf.shape[1]+j]

error = np.mean(np.square(mf - mf1))	

print "Mean Squared Error between MFCCs produced by libmfcc and python_speech_features is : ", error
