import h5py
import numpy as np
import nibabel as nib
import scipy.io
from scipy import ndimage
import time

def padding_CNN(data, pad_size):
    [numF, dimX, dimY, dimZ] = data.shape
    # print(data.shape)
    pData = np.zeros([numF, dimX + 2 * pad_size, dimY + 2 * pad_size, dimZ + 2 * pad_size], dtype=np.float64)
    pData[:, pad_size: dimX + pad_size, pad_size: dimY + pad_size, pad_size: dimZ + pad_size] = data
    return pData

def padding_CNN_diff_size(data, padX_l, padX_r, padY_l, padY_r, padZ_l, padZ_r):
    [numF, dimX, dimY, dimZ] = data.shape
    # print(data.shape)
    pData = np.zeros([numF, dimX + padX_l + padX_r, dimY + padY_l + padY_r, dimZ + padZ_l + padZ_r], dtype=np.float64)
    pData[:, padX_l: dimX + padX_l, padY_l: dimY + padY_l, padZ_l: dimZ + padZ_l] = data
    return pData

def conv_forward(data, w, b, stride_size, pad_size):

    [numF, dimX, dimY, dimZ] = data.shape
    pData=padding_CNN(data, pad_size)

    [numF, dimX, dimY, dimZ] = pData.shape
    [numCurrentF, numPrevF, dimkX, dimkY, dimkZ] = w.shape
    #print(w.shape)


    convData=np.zeros([numCurrentF, (dimX-dimkX)//stride_size+1, (dimY-dimkY)//stride_size+1, (dimZ-dimkZ)//stride_size+1], dtype=np.float64)

    for i in range(0, numCurrentF):

        centerX = 0

        startX = 0
        endX = startX + dimkX
        while endX<=dimX:
            startY = 0
            endY = startY + dimkY
            centerY = 0
            while endY<=dimY:
                startZ = 0
                endZ = startZ + dimkZ
                centerZ = 0
                while endZ<=dimZ:
                    for k in range(0, numF):
                        tempData = np.squeeze(pData[k, startX : endX, startY : endY, startZ : endZ])
                        convData[i, centerX, centerY, centerZ] = convData[i, centerX, centerY, centerZ] + np.sum(np.multiply(tempData, np.squeeze(w[i, k, :, :, :])))

                    convData[i, centerX, centerY, centerZ] = convData[i, centerX, centerY, centerZ] + b[i]

                    startZ = startZ + stride_size
                    endZ = endZ + stride_size
                    centerZ = centerZ + 1
                        #print(endZ)

                startY = startY + stride_size
                endY = endY + stride_size
                centerY = centerY + 1

            startX = startX + stride_size
            endX = endX + stride_size
            centerX = centerX + 1

    return convData

def ReLU(x):
    return abs(x) * (x > 0)


def maxPooling(data, kernel_size, stride_size, pad_size):

    [numF, dimX, dimY, dimZ] = data.shape

    pData=padding_CNN(data, pad_size)

    [numF, dimX, dimY, dimZ] = pData.shape

    if (dimX - kernel_size) % stride_size == 0:
        padX=0
    else:
        padX=kernel_size-(dimX - kernel_size) % stride_size

    if (dimY - kernel_size) % stride_size == 0:
        padY = 0
    else:
        padY=kernel_size-(dimY - kernel_size) % stride_size

    if (dimZ - kernel_size) % stride_size == 0:
        padZ=0
    else:
        padZ=kernel_size-(dimZ - kernel_size) % stride_size

    pData2=padding_CNN_diff_size(pData, 0, padX, 0, padY, 0, padZ)
    [numF, dimX, dimY, dimZ] = pData2.shape

    poolData=np.zeros([numF, (dimX-kernel_size)//stride_size+1, (dimY-kernel_size)//stride_size+1, (dimZ-kernel_size)//stride_size+1], dtype=np.float64)
    for k in range(0, numF):
        centerX = 0
        startX = 0
        endX = startX + kernel_size
        while endX<=dimX:
            startY = 0
            endY = startY + kernel_size
            centerY = 0
            while endY<=dimY:
                startZ = 0
                endZ = startZ + kernel_size
                centerZ = 0
                while endZ<=dimZ:
                    tempData = np.squeeze(pData2[k, startX : endX, startY : endY, startZ : endZ])
                    poolData[k, centerX, centerY, centerZ] = tempData.max()

                    startZ = startZ + stride_size
                    endZ = endZ + stride_size
                    centerZ = centerZ + 1

                startY = startY + stride_size
                endY = endY + stride_size
                centerY = centerY + 1

            startX = startX + stride_size
            endX = endX + stride_size
            centerX = centerX + 1

    return poolData


def FC_forward(data, w, b):

    [numF, dimX, dimY, dimZ] = data.shape
    [numbF] = b.shape
    flat_Data = data.reshape(numF * dimX * dimY * dimZ, 1)
    weighted_Data = np.matmul(w, flat_Data)
    weighted_Data = np.add(weighted_Data, b.reshape(numbF,1))
    [dimwX, dimwY]=weighted_Data.shape

    fcData = np.zeros([ 1, 1, 1, dimwX], dtype = np.float64)
    fcData [0, 0, 0, :] = weighted_Data[:,0]
    return fcData

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    return np.exp(x) / np.sum(np.exp(x), axis=0)



def feedforward(data, f):
    [dimX, dimY, dimZ] = data.shape
    netData = np.zeros([1, dimX, dimY, dimZ], dtype=np.float64)
    netData[0, :, :, :] = data

    numConv = 1
    numFC = 1

    while 'conv%01d_w' % numConv in f.keys():
        w = np.array(f.get('conv%01d_w' % numConv))
        b = np.array(f.get('conv%01d_b' % numConv))
        kernel_size = np.array(f.get('conv%01d_kernel' % numConv))
        stride_size = np.array(f.get('conv%01d_stride' % numConv))
        pad_size = np.array(f.get('conv%01d_pad' % numConv))

        netData = conv_forward(netData, w, b, stride_size, pad_size)
        netData = ReLU(netData)

        kernel_size = np.array(f.get('pool%01d_kernel' % numConv))
        stride_size = np.array(f.get('pool%01d_stride' % numConv))
        pad_size = np.array(f.get('pool%01d_pad' % numConv))

        netData = maxPooling(netData, kernel_size, stride_size, pad_size)

        numConv = numConv + 1

    while 'f%01d_w' % numFC in f.keys():
        w = np.array(f.get('f%01d_w' % numFC))
        b = np.array(f.get('f%01d_b' % numFC))

        netData = FC_forward(netData, w, b)

        if 'f%01d_w' % (numFC + 1) in f.keys():
            netData = ReLU(netData)

        numFC = numFC + 1

    prob = softmax(np.squeeze(netData))

    estLabel = prob.argmax(axis=0)
    return estLabel