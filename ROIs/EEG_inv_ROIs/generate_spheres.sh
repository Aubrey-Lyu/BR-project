#!/bin/bash
ROIcord=("71 1 51 1 62 1 0 1" "43 1 59 1 34 1 0 1" "59 1 64 1 28 1 0 1" \
"35 1 84 1 28 1 0 1" "16 1 50 1 58 1 0 1" "41 1 49 1 50 1 0 1" "16 1 50 1 58 1 0 1"\
 "56 1 60 1 29 1 0 1" "29 1 61 1 28 1 0 1" "51 1 81 1 44 1 0 1" "35 1 31 1 60 1 0 1" \
 "74 1 51 1 50 1 0 1" "53 1 30 1 41 1 0 1" "42 1 58 1 33 1 0 1" "32 1 76 1 32 1 0 1" "57 1 61 1 27 1 0 1" \
 "48 1 37 1 47 1 0 1" "53 1 30 1 41 1 0 1" "42 1 58 1 33 1 0 1" "32 1 76 1 32 1 0 1" "57 1 61 1 27 1 0 1" \
 "48 1 37 1 47 1 0 1" "45 1 52 1 37 1 0 1" "75 1 49 1 50 1 0 1" "40 1 26 1 62 1 0 1" "26 1 62 1 30 1 0 1" \
 "52 1 39 1 38 1 0 1" "38 1 48 1 32 1 0 1" "64 1 51 1 36 1 0 1" "59 1 62 1 28 1 0 1" \
 "47 1 35 1 58 1 0 1" "41 1 51 1 56 1 0 1" "47 1 62 1 53 1 0 1" "35 1 29 1 57 1 0 1" "26 1 58 1 23 1 0 1" \
 "65 1 46 1 59 1 0 1" "26 1 63 1 29 1 0 1" "34 1 29 1 57 1 0 1" "35 1 29 1 57 1 0 1" "60 1 32 1 56 1 0 1" "65 1 46 1 59 1 0 1" \
 "26 1 63 1 29 1 0 1" "34 1 29 1 57 1 0 1")
ROIname=("BR-Rpl_mixed_pos_IPL_inv1" "BR-Rpl_rdgn_neg_THALAMUS_inv1" "BR-Rpl_rdgn-mixed_PARAHPCAMPAL_inv1"\
 "BR-Rpl_rdgn-mixed_ACC_inv1" "mixed_BR_IPL_inv2" "mixed_BR_PCC_inv2" "BR-Rpl_mixed_pos_IPL_inv2" \
"BR-Rpl_rdgn_pos_lPARAHPCAMPAL_inv2" "BR-Rpl_rdgn_pos_rPARAHPCAMPAL_inv2" "BR-Rpl_rdgn_neg_ACC_inv2" "BR-Rpl_rdgn_neg_PCU_inv2"\
 "BR-Rpl_rdgn-mixed_IPL_inv2" "mixed_BR_PCC_inv3" "BR-Rpl_rdgn_neg_THALAMUS_inv3" "BR-Rpl_rdgn_neg_CLAUSTRUM_inv3" "BR-Rpl_rdgn-mixed_PARAHPCAMPAL_inv3" \
 "BR-Rpl_rdgn-mixed_PCC_inv3" "mixed_BR_PCC_inv4" "BR-Rpl_rdgn_neg_THALAMUS_inv4" "BR-Rpl_rdgn_neg_CLAUSTRUM_inv4" "BR-Rpl_rdgn-mixed_PARAHPCAMPAL_inv4" \
 "BR-Rpl_rdgn-mixed_PCC_inv4" "rdgn_BR_THALAMUS_inv5" "BR-Rpl_mixed_pos_IPL_inv5" "BR-Rpl_mixed_pos_PCU_inv5" "BR-Rpl_mixed_pos_rCLAUSTRUM_inv5" \
 "BR-Rpl_mixed_pos_PARAHPCAMPAL_inv5" "BR-Rpl_mixed_pos_THALAMUS_inv5" "BR-Rpl_mixed_pos_lCLAUSTRUM_inv5" "BR-Rpl_rdgn_pos_PARAHMCAMPAL_inv5" \
 "BR-Rpl_rdgn-mixed_PCU_inv5" "BR-Rpl_rdgn-mixed_PCC_inv5" "BR-Rpl_rdgn-mixed_midCC_inv5" "mixed_BR_PCU_inv6" "mixed_BR_PARAHPCAMPAL_inv6" \
 "BR-Rpl_mixed_pos_IPL_inv6" "BR-Rpl_rdgn-mixed_CLAUSTRUM_inv6" "BR-Rpl_rdgn-mixed_PCU_inv6" "mixed_BR_PCU_inv7" "mixed_BR_AG_inv7" "BR-Rpl_mixed_pos_IPL_inv7" \
 "BR-Rpl_rdgn-mixed_CLAUSTRUM_inv7" "BR-Rpl_rdgn-mixed_PCU_inv7")


for index in ${!ROIname[*]}; do 
fslmaths $FSLDIR/data/standard/avg152T1.nii.gz -mul 0 -add 1 -roi ${ROIcord[$index]} ${ROIname[$index]} -odt float
fslmaths ${ROIname[$index]} -kernel sphere 8 -fmean "${ROIname[$index]}_sphere" -odt float
done
