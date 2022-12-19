#!/bin/bash


########################################################################################################################################################
# If you have questions in this part, please ask Xuyun Wen for help
# Head Motion Correction, Distortion Correction, EPI to T1 registration, 
# One-step resampling, high-pass filtering and ICA decomposition, 
# Done By Xuyun Wen

i=1
DATAPATH="/media/zz/data/test_BCPpipeline/test_data/All/After_Sorted"
OUTPUTPATH="/media/zz/data/test_BCPpipeline/test_data/All/Functional_Preprocessed"
FILEPATH="/media/zz/data/test_BCPpipeline/test_data/All/MNI_Template"
STRUCPATH="/media/zz/data/test_BCPpipeline/test_data/All/Structural_Preprocessed"
ATLASPATH="/media/zz/data/test_BCPpipeline/test_data/All/NativeSpace_atlas"
OUTPUTRegPic="/media/zz/data/test_BCPpipeline/test_data/All/Structural_Preprocessed/PICforRegCheckMNItoT2"
ICANORMPATH="/media/zz/data/test_BCPpipeline/test_data/All/ICA_Normalized"
TIMECOURSEPATH="/media/zz/data/test_BCPpipeline/test_data/All/TIMECOURSE"
NETWORK="/media/zz/data/test_BCPpipeline/test_data/All/NETWORK"
CLEANDATAPATH="/media/zz/data/test_BCPpipeline/test_data/All/AFTERDENOISE"
MINIOUTPUTPATH="/media/zz/data/test_BCPpipeline/test_data/All/Minimal_Preprocessed"
mkdir -p ${OUTPUTPATH}
mkdir -p ${STRUCPATH}
mkdir -p ${ATLASPATH}
mkdir -p ${OUTPUTRegPic}
mkdir -p ${ICANORMPATH}
mkdir -p ${TIMECOURSEPATH}
mkdir -p ${CLEANDATAPATH}


echo DataPreprocessing Sample $i

# unzip 4D fMRI image into multiple 3D images
OUTPUTFMRI="${MINIOUTPUTPATH}/Sample${i}"
mkdir -p ${OUTPUTFMRI}
mkdir -p ${OUTPUTFMRI}/prevols    # fmri before preprocessing
fMRI="${DATAPATH}/Sample${i}/REST.nii.gz"
fslsplit ${fMRI} ${OUTPUTFMRI}/prevols/vol -t


#Head Motion Correction using mcflirt
echo Head Motion Correction Sample $i
#if SBref image exists,  SBref.nii is used as the reference image; 
#otherwise, the 10th volume of fmri image is used as the reference image
SBREFNAME="${DATAPATH}/Sample${i}/sbref.txt"
SBRefFlag=$(cat $SBREFNAME)
if [[ $SBRefFlag == "1" ]]; then
SBRef="${DATAPATH}/Sample${i}/SBRef.nii.gz"
elif [[ $SBRefFlag == "0" ]]; then
SBRef="${OUTPUTFMRI}/prevols/vol0009.nii.gz"
fi

OUTPUTMC="${OUTPUTPATH}/Sample${i}/MC"
mkdir -p $OUTPUTMC
mcflirt -in ${fMRI} -r ${SBRef} -mats -plots -o ${OUTPUTMC}/mc

#EPI Distortion Correction using topup function
echo EPI Distortion Correction Subject ${i}
OUTPUTDC="${OUTPUTPATH}/Sample${i}/DC"
mkdir -p $OUTPUTDC
# Read FieldMAP AP and PA and combine them together
PhaseAP="${DATAPATH}/Sample${i}/FieldMapAP.nii.gz"
PhasePA="${DATAPATH}/Sample${i}/FieldMapPA.nii.gz"
fslmerge -t ${OUTPUTDC}/BothPhases ${PhaseAP} ${PhasePA}
fslmaths ${PhaseAP} -mul 0 -add 1 ${OUTPUTDC}/Mask
# Read Topup configuration files
TopupConfig="${DATAPATH}/b02b0.cnf"
txtfname="${DATAPATH}/acqparams.txt"
# Use topup function to do distortion correction
numslice=`fslval ${OUTPUTDC}/BothPhases dim3`
if [ ! $(($numslice % 2)) -eq "0" ] ; then
  log_Msg "Padding Z by one slice"
  echo Padding Z by one slice
  for Image in ${OUTPUTDC}/BothPhases ${OUTPUTDC}/Mask ; do
    fslroi ${Image} ${OUTPUTDC}/slice.nii.gz 0 -1 0 -1 0 1 0 -1
    fslmaths ${OUTPUTDC}/slice.nii.gz -mul 0 ${OUTPUTDC}/slice.nii.gz
    fslmerge -z ${Image} ${Image} ${OUTPUTDC}/slice.nii.gz
    rm ${OUTPUTDC}/slice.nii.gz
  done
fi
fslmaths ${OUTPUTDC}/BothPhases -abs -add 1 -mas ${OUTPUTDC}/Mask -dilM -dilM -dilM -dilM -dilM ${OUTPUTDC}/BothPhases
topup --imain=${OUTPUTDC}/BothPhases --datain=$txtfname --config=$TopupConfig --out=${OUTPUTDC}/Coefficents --iout=${OUTPUTDC}/Magnitudes --fout=${OUTPUTDC}/TopupField --dfout=${OUTPUTDC}/WarpField --rbmout=${OUTPUTDC}/MotionMatrix --jacout=${OUTPUTDC}/Jacobian -v

#  register SBRef to FieldMap
#  use Phase to distinguish PA and AP; if x: AP, if y: PA
dimtOne=`${FSLDIR}/bin/fslval ${PhaseAP} dim4`
dimtTwo=`${FSLDIR}/bin/fslval ${PhasePA} dim4`
PHASENAME="${DATAPATH}/Sample${i}/Phase.txt"  # store a label to indicate the AP (x) or PA (y)
PhaseD=$(cat $PHASENAME)
if [[ $PhaseD == "x" ]]; then   # the phase direction is AP
   VolumeNumber=$((0 + 1))
   vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
   flirt -dof 6 -interp spline -in $SBRef -ref $PhaseAP -omat ${OUTPUTDC}/SBRef2PhaseAP.mat -out ${OUTPUTDC}/SBRef2PhaseAP
   convert_xfm -omat ${OUTPUTDC}/SBRef2WarpField.mat -concat ${OUTPUTDC}/MotionMatrix_${vnum}.mat ${OUTPUTDC}/SBRef2PhaseAP.mat
   convertwarp --relout --rel -r $PhaseAP --premat=${OUTPUTDC}/SBRef2WarpField.mat --warp1=${OUTPUTDC}/WarpField_${vnum} --out=${OUTPUTDC}/WarpField.nii.gz
   imcp ${OUTPUTDC}/Jacobian_${vnum}.nii.gz ${OUTPUTDC}/Jacobian.nii.gz

elif [[ $PhaseD == "y" ]]; then # the phase direction is PA
   VolumeNumber=$(($dimtOne + 1))
   vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
   flirt -dof 6 -interp spline -in $SBRef -ref $PhasePA -omat ${OUTPUTDC}/SBRef2PhasePA.mat -out ${OUTPUTDC}/SBRef2PhasePA
   convert_xfm -omat ${OUTPUTDC}/SBRef2WarpField.mat -concat ${OUTPUTDC}/MotionMatrix_${vnum}.mat ${OUTPUTDC}/SBRef2PhasePA.mat
   convertwarp --relout --rel -r $PhasePA --premat=${OUTPUTDC}/SBRef2WarpField.mat --warp1=${OUTPUTDC}/WarpField_${vnum} --out=${OUTPUTDC}/WarpField.nii.gz 
   imcp ${OUTPUTDC}/Jacobian_${vnum}.nii.gz ${OUTPUTDC}/Jacobian.nii.gz
fi

# PhaseTwo (first vol) - warp and Jacobian modulate to get distortion corrected output
echo warp and Jacobian modulate to get distortion corrected output sample ${i}
echo warp and Jacobian modulate to get distortion corrected PhasePA
VolumeNumber=$(($dimtOne + 1))
vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
applywarp --rel --interp=spline -i ${PhasePA} -r ${PhasePA} --premat=${OUTPUTDC}/MotionMatrix_${vnum}.mat -w ${OUTPUTDC}/WarpField_${vnum} -o ${OUTPUTDC}/PhasePA_dc
fslmaths ${OUTPUTDC}/PhasePA_dc -mul ${OUTPUTDC}/Jacobian_${vnum} ${OUTPUTDC}/PhasePA_dc_jac
# PhaseOne (first vol) - warp and Jacobian modulate to get distortion corrected output
echo warp and Jacobian modulate to get distortion corrected PhaseAP
VolumeNumber=$((0 + 1))
vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
applywarp --rel --interp=spline -i ${PhaseAP} -r ${PhaseAP} --premat=${OUTPUTDC}/MotionMatrix_${vnum}.mat -w ${OUTPUTDC}/WarpField_${vnum} -o ${OUTPUTDC}/PhaseAP_dc
fslmaths ${OUTPUTDC}/PhaseAP_dc -mul ${OUTPUTDC}/Jacobian_${vnum} ${OUTPUTDC}/PhaseAP_dc_jac

# Scout - warp and Jacobian modulate to get distortion corrected output
echo warp and Jacobian modulate to get distortion corrected SBRef
applywarp --rel --interp=spline -i ${SBRef} -r ${SBRef} -w ${OUTPUTDC}/WarpField.nii.gz -o ${OUTPUTDC}/SBRef_dc.nii.gz
fslmaths ${OUTPUTDC}/SBRef_dc.nii.gz -mul ${OUTPUTDC}/Jacobian.nii.gz ${OUTPUTDC}/SBRef_dc_jac.nii.gz

# Calculate Equivalent Field Map
fslmaths ${OUTPUTDC}/TopupField -mul 6.283 ${OUTPUTDC}/TopupField
fslmaths ${OUTPUTDC}/Magnitudes.nii.gz -Tmean ${OUTPUTDC}/Magnitude.nii.gz
bet ${OUTPUTDC}/Magnitude ${OUTPUTDC}/Magnitude_brain -f 0.35 -m   #Brain extract the magnitude image


# Register EPI to T1
echo register EPI to T1 sample ${i}
OUTPUTEPI2STR="${OUTPUTPATH}/Sample${i}/EPI2STR"
mkdir -p $OUTPUTEPI2STR
StructureIntensityImage=${DATAPATH}/Sample${i}/T1stripped*   # Please check the structural data type in raw data file. The data type from Zhengwang may be .nii or .nii.gz  
# create the white matter segetation from tissue map
echo create the white matter issue sample ${i}
T2wseg=${DATAPATH}/Sample${i}/tissuelabel*                   # Please check the structural data type in raw data file. The data type from Zhengwang may be .nii or .nii.gz 
fslmaths ${T2wseg} -thr 2.5 -div 3 ${OUTPUTEPI2STR}/T2wwmseg.nii.gz

# use flirt do registration from SBRef to tissuelabel image
vout=Epi2Str
dof=6
# do a standard flirt pre-alignment
echo "FLIRT pre-alignment"
flirt -ref ${StructureIntensityImage} -in ${OUTPUTDC}/SBRef_dc_jac.nii.gz -dof ${dof} -omat ${OUTPUTEPI2STR}/${vout}_init.mat -out ${OUTPUTEPI2STR}/${vout}_init -cost mutualinfo -searchcost mutualinfo -searchrx -180 180 -searchry -180 180 -searchrz -180 180
# do the second time regisration using bbr
echo "Running BBR"
flirt -ref ${StructureIntensityImage} -in ${OUTPUTDC}/SBRef_dc_jac.nii.gz -dof ${dof} -cost bbr -wmseg ${OUTPUTEPI2STR}/T2wwmseg.nii.gz -init ${OUTPUTEPI2STR}/${vout}_init.mat -applyxfm -omat ${OUTPUTEPI2STR}/${vout}.mat -out ${OUTPUTEPI2STR}/${vout} -schedule ${FSLDIR}/etc/flirtsch/bbr.sch
applywarp -i ${OUTPUTDC}/SBRef_dc_jac.nii.gz -r ${StructureIntensityImage} -o ${OUTPUTEPI2STR}/${vout} --premat=${OUTPUTEPI2STR}/${vout}.mat --interp=spline

cp ${OUTPUTEPI2STR}/${vout}.mat ${OUTPUTEPI2STR}/fMRI2str.mat
#generate combined warpfields and spline interpolated images
convertwarp --relout --rel -r ${StructureIntensityImage} --warp1=${OUTPUTDC}/WarpField.nii.gz --postmat=${OUTPUTEPI2STR}/${vout}.mat -o ${OUTPUTEPI2STR}/${vout}_warp
applywarp --rel --interp=spline -i ${OUTPUTDC}/Jacobian.nii.gz -r ${StructureIntensityImage} --premat=${OUTPUTEPI2STR}/${vout}.mat -o ${OUTPUTEPI2STR}/Jacobian2T2w.nii.gz
# 1-step resample from input (gdc) scout - NOTE: no longer includes jacobian correction, if specified
applywarp --rel --interp=spline -i ${SBRef} -r ${StructureIntensityImage} -w ${OUTPUTEPI2STR}/${vout}_warp -o ${OUTPUTEPI2STR}/SBRef_undistorted2T1w


# one step resampling
# first: down resampling T1 image
RefImgforOneStep=${OUTPUTFMRI}/T2DownResample
flirt -in ${StructureIntensityImage} -ref ${StructureIntensityImage} -o ${RefImgforOneStep} -applyisoxfm 2

echo One Step Resampling sample ${i}
OUTPUTREAMPLE="${OUTPUTPATH}/Sample${i}/ONESTEPREAMPLING"
mkdir -p ${OUTPUTREAMPLE}
mkdir -p ${OUTPUTFMRI}/postvols   # fmri after preprocessing

OUTPUTREGQ="${OUTPUTPATH}/PICforREGCHECK/Sample${i}"
mkdir -p ${OUTPUTREGQ}

TR_vol=`${FSLDIR}/bin/fslval ${fMRI} pixdim4 | cut -d " " -f 1`
NumFrames=`${FSLDIR}/bin/fslval ${fMRI} dim4`
# Apply combined transformations to fMRI (combines gradient non-linearity distortion, motion correction, and registration to T1w space, but keeping fMRI resolution)
FrameMergeSTRING=""
FrameMergeSTRINGII=""
k=0
while [ $k -lt $NumFrames ] ; do
 vnum=`${FSLDIR}/bin/zeropad $k 4`
 echo $k
 prevmatrix="${OUTPUTMC}/mc.mat"
 convertwarp --relout --rel --ref=${RefImgforOneStep} --warp1=${OUTPUTEPI2STR}/${vout}_warp --premat=$prevmatrix/MAT_${vnum} --out=${OUTPUTREAMPLE}/vol${vnum}_all_warp.nii.gz
 fslmaths ${OUTPUTFMRI}/prevols/vol${vnum}.nii.gz -mul 0 -add 1 ${OUTPUTFMRI}/prevols/vol${vnum}_mask.nii.gz
 applywarp --rel --interp=spline --in=${OUTPUTFMRI}/prevols/vol${vnum}.nii.gz --warp=${OUTPUTREAMPLE}/vol${vnum}_all_warp.nii.gz --ref=${RefImgforOneStep} --out=${OUTPUTFMRI}/postvols/vol${k}.nii.gz
 applywarp --rel --interp=nn --in=${OUTPUTFMRI}/prevols/vol${vnum}_mask.nii.gz --warp=${OUTPUTREAMPLE}/vol${vnum}_all_warp.nii.gz --ref=${RefImgforOneStep} --out=${OUTPUTFMRI}/postvols/vol${k}_mask.nii.gz
 FrameMergeSTRING="${FrameMergeSTRING}${OUTPUTFMRI}/postvols/vol${k}.nii.gz " 
 FrameMergeSTRINGII="${FrameMergeSTRINGII}${OUTPUTFMRI}/postvols/vol${k}_mask.nii.gz " 
 # output images for registration quality check
 slicer ${OUTPUTFMRI}/postvols/vol${k}.nii.gz ${RefImgforOneStep} -S 2 3000 ${OUTPUTREGQ}/vol_${vnum}.png

 k=`echo "$k + 1" | bc`
done

echo merge multiple 3D image into one 4D image
fslmerge -tr ${OUTPUTFMRI}/fMRIAfterMinP $FrameMergeSTRING $TR_vol
fslmerge -tr ${OUTPUTFMRI}/fMRIAfterMinP_mask $FrameMergeSTRINGII $TR_vol
fslmaths ${OUTPUTFMRI}/fMRIAfterMinP_mask -Tmin ${OUTPUTFMRI}/fMRIAfterMinP_mask

# high-pass filter
echo High Pass Filter
OUTPUTFILTER="${OUTPUTPATH}/Sample${i}/FILTER"
mkdir -p ${OUTPUTFILTER}
# calculate the sigma to volume: sigma = 1/(2*tr*(1/cut in seconds))
cut=1000
SIGMA=`echo "${cut}/${TR_vol}" | bc`
echo sigma2volume =$SIGMA
fslmaths ${OUTPUTFMRI}/fMRIAfterMinP.nii.gz -bptf ${SIGMA} -1 ${OUTPUTFILTER}/fMRIAfterfilter.nii.gz

echo calculate the mean fmri image
fslmaths ${OUTPUTFMRI}/fMRIAfterMinP.nii.gz -Tmean ${OUTPUTFILTER}/fMRIMean.nii.gz
fslmaths ${OUTPUTFILTER}/fMRIAfterfilter.nii.gz -add ${OUTPUTFILTER}/fMRIMean.nii.gz ${OUTPUTFILTER}/fMRIAfterfilterAddMean.nii.gz

# ICA decomposition
 OUTPUTICA="${OUTPUTPATH}/Sample${i}/ICA"
mkdir -p ${OUTPUTICA}
MaskforICA="${OUTPUTICA}/BrainMask"
fslmaths $RefImgforOneStep -bin ${MaskforICA}
melodic -i ${OUTPUTFILTER}/fMRIAfterfilterAddMean.nii -o ${OUTPUTICA}/ICAResult -m ${MaskforICA} --bgimage=${RefImgforOneStep} --nobet -d 150 --report -v --tr=${TR_vol} --report_maps="-S 2 2000 "


###########################################################################################################################################################
# For this part, if you have questions, please ask Bin Jing for help
# Regisration from Structure images to MNI space,  DONE by Bin Jing. 
MNILabelImage="${FILEPATH}/label_template_mni_sym_addsymsubcort_masked.nii" 

LabelRunRegName="${DATAPATH}/Sample${i}/EPItoMNILabel.txt"
RunRegFlag=$(cat $LabelRunRegName)

# Registration from native space to MNI space
#if [[ $RunRegFlag == "1" ]]; 
#then
  echo Registration from native space to MNI space
  REGMNItoT2="${STRUCPATH}/Sample${i}"
  mkdir -p ${REGMNItoT2}
  antsRegistrationSyNQuick.sh -d 3 -f ${MNILabelImage} -m ${T2wseg} -o ${REGMNItoT2}/Str2MNIAnt -n 8
  echo Generate NeareatNeighbour interpolated normalized image
  antsApplyTransforms -d 3 -i ${T2wseg} -r ${MNILabelImage} -o ${REGMNItoT2}/Str2MNIAntWarped_NN.nii.gz -t ${REGMNItoT2}/Str2MNIAnt1Warp.nii.gz -t ${REGMNItoT2}/Str2MNIAnt0GenericAffine.mat  --interpolation NearestNeighbor

  echo CreateWarpedTemplate Sample $i
  OUTPUTNativeSpaceATLAS="${ATLASPATH}/Sample${i}"
  mkdir -p ${OUTPUTNativeSpaceATLAS}
  antsApplyTransforms -d 3 -i ${FILEPATH}/new_HarvardOxford_2mm_rmoverlap.nii -o ${OUTPUTNativeSpaceATLAS}/Harvard_warped.nii.gz -r ${OUTPUTFILTER}/fMRIMean.nii.gz -t [${REGMNItoT2}/Str2MNIAnt0GenericAffine.mat,1] -t ${REGMNItoT2}/Str2MNIAnt1InverseWarp.nii.gz --interpolation NearestNeighbor
  antsApplyTransforms -d 3 -i ${FILEPATH}/HOA_Contract_1024_2MM_10.nii.gz -o ${OUTPUTNativeSpaceATLAS}/HOA1024_warped.nii.gz -r ${OUTPUTFILTER}/fMRIMean.nii.gz -t [${REGMNItoT2}/Str2MNIAnt0GenericAffine.mat,1] -t ${REGMNItoT2}/Str2MNIAnt1InverseWarp.nii.gz --interpolation NearestNeighbor
  antsApplyTransforms -d 3 -i ${FILEPATH}/Zalesky_980_parcellated_compact.nii -o ${OUTPUTNativeSpaceATLAS}/Zalesky980_warped.nii.gz -r ${OUTPUTFILTER}/fMRIMean.nii.gz -t [${REGMNItoT2}/Str2MNIAnt0GenericAffine.mat,1] -t ${REGMNItoT2}/Str2MNIAnt1InverseWarp.nii.gz --interpolation NearestNeighbor
  antsApplyTransforms -d 3 -i ${FILEPATH}/Zalesky_1024_parcellated_uniform.nii -o ${OUTPUTNativeSpaceATLAS}/Zalesky1024_warped.nii.gz -r ${OUTPUTFILTER}/fMRIMean.nii.gz -t [${REGMNItoT2}/Str2MNIAnt0GenericAffine.mat,1] -t ${REGMNItoT2}/Str2MNIAnt1InverseWarp.nii.gz --interpolation NearestNeighbor
  antsApplyTransforms -d 3 -i ${FILEPATH}/Comp_total.nii -o ${OUTPUTNativeSpaceATLAS}/Comp_total_warped.nii.gz -r ${OUTPUTFILTER}/fMRIMean.nii.gz -t [${REGMNItoT2}/Str2MNIAnt0GenericAffine.mat,1] -t ${REGMNItoT2}/Str2MNIAnt1InverseWarp.nii.gz --interpolation NearestNeighbor
  
  # save registration results
  echo Save RegistrationResults $i
  slicer ${REGMNItoT2}/Str2MNIAntWarped.nii ${MNILabelImage} -S 2 3000 ${OUTPUTRegPic}/Sample$i.png
#fi

echo Normalize ICA components of each indiviual to MNI space
OUTPUTNORMIC="${ICANORMPATH}/Sample${i}" 
Resliced_MNILabelImage="${FILEPATH}/Resliced_label_template_mni_sym_addsymsubcort_masked.nii" 
mkdir -p ${OUTPUTNORMIC}
#LabelICANormName="${DATAPATH}/Sample${i}/ICANormLabel.txt"
#ICANormRefSubID=$(cat $LabelICANormName)

antsApplyTransforms -d 3 -e 3 -i ${OUTPUTICA}/ICAResult/melodic_IC.nii.gz -o ${OUTPUTNORMIC}/ICA_MNI.nii.gz -r ${Resliced_MNILabelImage} -t ${STRUCPATH}/Sample${i}/Str2MNIAnt1Warp.nii.gz -t ${STRUCPATH}/Sample${i}/Str2MNIAnt0GenericAffine.mat 

################################################################################################################################################################
    # Data Cleaning using the classified bad ICA components. Component classification is done by Tae-Eui
    # Done by Tae-Eui
   echo Detect the noise related ICA components using Tae-Eui model
   OUTPUTCLEANDATAICA="${CLEANDATAPATH}/BadLabel"
   mkdir -p ${OUTPUTCLEANDATAICA}
   LabelName="${OUTPUTCLEANDATAICA}/Sample${i}.txt"
   python2 noise_comp_detection/bad_comp_est.py "${OUTPUTNORMIC}/ICA_MNI.nii.gz" "${LabelName}"
  

   # done by xuyun wen
   echo clean fMRI data sample ${i}
     OUTPUTCLEANDATA="${CLEANDATAPATH}/CLEANEDDATA/Sample${i}"
   mkdir -p ${OUTPUTCLEANDATA}
   Label=$(cat $LabelName)
   fsl_regfilt -i ${OUTPUTFILTER}/fMRIAfterfilterAddMean.nii.gz -o ${OUTPUTCLEANDATA}/fMRIAfterDenoise -d ${OUTPUTICA}/ICAResult/melodic_mix -f "${Label}"


    # Folder compression and removal
    cd ${OUTPUTPATH}/Sample${i}
    zip -r ICA.zip ICA/
    cd -

    rm -rv ${OUTPUTFMRI}/postvols 
    rm -rv ${OUTPUTREAMPLE}

    rm -rv ${OUTPUTMC}/mc.mat
    rm -rv ${OUTPUTICA}
	
	mv ${OUTPUTREGQ}/vol_0009.png ${OUTPUTPATH}/PICforREGCHECK/Sample${i}.png
    rm -rv ${OUTPUTREGQ}


#################################################################################################################################################################
#extract time series using four altases
#Done by Zhen Zhou
#OUTPUTCLEANDATA="${OUTPUTPATH}/AFTERDENOISE/CLEANEDDATA/Sample${i}"
fslmeants -i ${OUTPUTCLEANDATA}/fMRIAfterDenoise --label=${ATLASPATH}/Sample${i}/HOA1024_warped.nii.gz     -o ${TIMECOURSEPATH}/HOA1024_Sample${i}.txt
fslmeants -i ${OUTPUTCLEANDATA}/fMRIAfterDenoise --label=${ATLASPATH}/Sample${i}/Harvard_warped.nii.gz     -o ${TIMECOURSEPATH}/Harvard112_Sample${i}.txt
fslmeants -i ${OUTPUTCLEANDATA}/fMRIAfterDenoise --label=${ATLASPATH}/Sample${i}/Zalesky980_warped.nii.gz     -o ${TIMECOURSEPATH}/Zalesky980_Sample${i}.txt
fslmeants -i ${OUTPUTCLEANDATA}/fMRIAfterDenoise --label=${ATLASPATH}/Sample${i}/Zalesky1024_warped.nii.gz     -o ${TIMECOURSEPATH}/Zalesky1024_Sample${i}.txt

################################################################################################################################################################
#seed-based correlation analysis
NETWORKPATH="${NETWORK}/Sample${i}_network"
mkdir -p ${NETWORKPATH}
3dAutomask -prefix ./Sample${i}_mask.nii ${OUTPUTCLEANDATA}/fMRIAfterDenoise.nii.gz

for comp in 1 2 3 4 5 6 7 8 9 10
do
3dmaskave -mask ${OUTPUTNativeSpaceATLAS}/Comp_total_warped.nii.gz -mrange ${comp} ${comp} -quiet ${OUTPUTCLEANDATA}/fMRIAfterDenoise.nii.gz > Sample${i}_timecourse_${comp}.1D
3dDeconvolve -input ${OUTPUTCLEANDATA}/fMRIAfterDenoise.nii.gz -mask Sample${i}_mask.nii -jobs 8 -float -num_stimts 1 -stim_file 1 Sample${i}_timecourse_${comp}.1D -stim_label 1 "seeds_${comp}" -tout -rout -bucket Sample${i}_comp${comp}_buc
3dcalc -a Sample${i}_comp${comp}_buc+orig'[4]' -b Sample${i}_comp${comp}_buc+orig'[2]' -expr "ispositive(b)*sqrt(a)-isnegative(b)*sqrt(a)" -prefix Sample${i}_comp${comp}_buc_r
3dcalc -a Sample${i}_comp${comp}_buc_r+orig -expr "0.5*log((1+a)/(1-a))" -datum float -prefix Sample${i}_comp${comp}_buc_r_z.nii
done


#cd ./Network/
mv Sample${i}_*_z.nii ${NETWORK}/Sample${i}_network/
mv Sample${i}_mask.nii ${NETWORK}/Sample${i}_network/
rm  Sample${i}_*

sleep 1
