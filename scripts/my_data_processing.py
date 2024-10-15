import nibabel as nib
import numpy as np

def read_nii_data(filePath):
    if not filePath.endswith('nii.gz') and not filePath.endswith('nii'):
        raise Exception('The Data file should be nii or nii.gz format')

    img = nib.load(filePath)
    data = img.get_data()
    return data

def remove_predifined_boundary(data):
    minX = 8
    maxX = 90
    minY = 9
    maxY = 107
    minZ = 1
    maxZ = 83
    boundary_removed_data = np.zeros([1, 1, maxX - minX, maxY - minY, maxZ - minZ], dtype=np.float64)
    boundary_removed_data = np.array(data[minX : maxX, minY : maxY, minZ : maxZ])
    return boundary_removed_data

def data_normalization(data):
    data[np.isinf(data)] = 0
    data[np.isneginf(data)] = 0
    data[np.isnan(data)] = 0

    voxInd = [data != 0]
    backInd = [data == 0]
    vox = data[voxInd]

    # data normalization
    meanVal = vox.mean()
    stdVal = vox.std()

    data = data - meanVal
    data = data / stdVal
    data[backInd] = 0

    minval = data.min()
    maxval = data.max()

    posInd = (data > 0)
    data[posInd] = data[posInd] / float(maxval)

    negInd = (data < 0)
    data[negInd] = data[negInd] / abs(float(minval))
    return data

def padding_3D_diff_size(data, padX_l, padX_r, padY_l, padY_r, padZ_l, padZ_r):
    [dimX, dimY, dimZ] = data.shape

    pData = np.zeros([dimX + padX_l + padX_r, dimY + padY_l + padY_r, dimZ + padZ_l + padZ_r], dtype=np.float64)
    pData[padX_l: dimX + padX_l, padY_l: dimY + padY_l, padZ_l: dimZ + padZ_l] = data
    return pData


def majority_voting_label(estLabel):
    [numModel, numLabel]=estLabel.shape
    votingLabel = estLabel.sum(axis=0)

    Th=(numModel + 1) // 2

    noise=votingLabel < Th
    signal=votingLabel >= Th

    votingLabel[noise] = 0
    votingLabel[signal] = 1

    return votingLabel

