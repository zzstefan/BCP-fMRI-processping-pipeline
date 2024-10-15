import nibabel as nib
import numpy as np

def load_nii(filePath):
    if not filePath.endswith('nii.gz') and not filePath.endswith('nii'):
        raise Exception('The Data file should be nii or nii.gz format')

    img = nib.load(filePath)
    data = img.get_data()
    return data


def export_noise_index(label, filePath):
    [numRow,] = label.shape

    num_noise=np.sum(label == 0)
    t = open(filePath, 'w')

    idx_noise=0
    for i in range(0, numRow):
        if label[i] == 0:
            if idx_noise == (num_noise - 1):
                t.write('%01d' % (i + 1))
            else:
                t.write('%01d, ' % (i + 1))
            idx_noise = idx_noise + 1
    t.close()
    return 0