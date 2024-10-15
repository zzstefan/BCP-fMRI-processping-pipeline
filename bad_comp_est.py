import numpy as np
import time
import sys
import h5py

import my_data_processing as dp
import my_data_rw as drw
import my_3DCNN as cnn

dataFileName=sys.argv[1]
outputFileName=sys.argv[2]
print(dataFileName)

timevalue = np.zeros([1,10], dtype=np.float16)

data=dp.read_nii_data(dataFileName)
[dimrX, dimrY, dimrZ, numComp] = data.shape

numModel=1


for k in range(0,10):
    start_time = time.time()
    estLabel=np.zeros([numModel, numComp], dtype=np.int)

    for idxModel in range(0, numModel):
        #load a CNN model
        f = h5py.File('noise_comp_detection/model/cnn_cross_validation_%02d.hdf5' %(idxModel + 1), 'r')

        for i in range(0, numComp):
            single_data = np.zeros([dimrX, dimrY, dimrZ], dtype=np.float64)
            single_data = data[:,:,:,i]
            boundary_removed_data = dp.remove_predifined_boundary(single_data)
            normalized_data = dp.data_normalization(boundary_removed_data)
            estLabel[idxModel, i]=cnn.feedforward(normalized_data,f)

    votingLabel=dp.majority_voting_label(estLabel)
    drw.export_noise_index(votingLabel, outputFileName)

    timevalue[0,k]=time.time() - start_time

    print(timevalue)
print(np.mean(timevalue))
print(np.std(timevalue))








