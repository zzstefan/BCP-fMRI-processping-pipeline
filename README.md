# BCP-fMRI-processping-pipeline

The BCP fMRI data processing pipeline is in the form of a bash script. It is easy to run for someone with knowledge of the BASH scripting language. It is an infant-dedicated pipeline[1-2] which include the following steps: head-motion and distortion correction, infant-dedicated spatial registration, one-time resampling of fMRI data and individual independent component analysis (ICA)-based extensive denoising[3], similar to the state-of-the-art HCP pipeline.

The pipeline can be run on a Linux system computer (My own PC is Ubuntu system and the pipeline runs well), including the high-performance-computing clusters. One must install ANTS (ver. 2.3.1, http://stnava.github.io/ANTs/), FSL (ver. 5.0.10/5.0.11, ver. 6 may causes some mistakes, https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation), AFNI (ver. 17.2.09, https://afni.nimh.nih.gov) and python2 to run the script. If you run it on a cluster using slurm to submit jobs, simply just add these modules each time you run the code (using ‘module add’ command). Note that by slurm, you can use command ‘SLURM_ARRAY_TASK_ID’ to run hundreds of subjects simultaneously. If run it on your own PC, use the for loop to preprocess data one by one. The estimated time for running one sample is approximately 6 hours in my own PC.


References:
1.	Zhou, Z., et al., Multi-layer Temporal Network Analysis Reveals Increasing Temporal Reachability and Spreadability in the First Two Years of Life, in Medical Image Computing and Computer Assisted Intervention – MICCAI 2019. 2019. p. 665-672.
2.	Hu, D., Wang, F., Zhang, H., Wu, Z., Zhou, Z., Li, G., Wang, L., Lin, W., Li, G. and UNC/UMN Baby Connectome Project Consortium, 2022. Existence of Functional Connectome Fingerprint during Infancy and Its Stability over Months. Journal of Neuroscience, 42(3), pp.377-389.
3.	Kam, T.-E., et al. A Deep Learning Framework for Noise Component Detection from Resting-State Functional MRI. 2019. Cham: Springer International Publishing.
