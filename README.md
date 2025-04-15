# MarkVCID Free Water Kit

**Overview**
In the context of the MarkVCID consortium, NINDS leadership selected 5 candidate imaging and fluid-based biomarkers to undergo clinical validation in MarkVCID2 to the point of readiness for incorporation into future VCID interventional trials. The MRI Free Water (FW) kit, developed by UCD and UNM, is one of these biomarkers. The FW kit demonstrates remarkable analytical performances with excellent inter-
rater reliability (ICC=0.997, confidence interval: [0.993; 0.999], p<0.001;), test-retest repeatability (0.995 [0.99; 0.997], p<0.001) and inter-scanner reproducibility (3 scanners: 0.96 [0.911; 0.985], p<0.001; 4 scanners: 0.964 [0.926; 0.986], p<0.001) of mFW measures (Maillard et al. Alzheimers Dement (Amst). 2022a). FW measure was also found to be strongly associated with cognition in the MarkVCID 1 cohort (Maillard et al. Alzheimers Dement (Amst). 2022b).

Maillard P, Lu H, Arfanakis K, et al. Instrumental validation of free water, peak-width of skeletonized mean diffusivity, and white matter hyperintensities: MarkVCID neuroimaging kits. Alzheimer's & Dementia:
Diagnosis, Assessment & Disease Monitoring 2022;14:e12261.

Maillard P, Lu H, Arfanakis K, Gold BT, Bauer CE, Zachariou V, et al. Instrumental validation of free water, peak-width of skeletonized mean diffusivity, and white matter hyperintensities: MarkVCID neuroimaging kits. Alzheimers Dement (Amst). 2022 Mar 28;14(1):e12261. doi: 10.1002/dad2.12261. eCollection 2022b. PMID: 35382232

# **I.	Components**
Different files are provided in the distributed directory /scripts_FW_CONSORTIUM/:  
•	MAIN_script_FW.sh: batch that launches the main shell script,  
•	mrn.py: script that computes the FW map  
•	FMRIB58_FA_1mm_thr: White matter mask (FSL FMRIB58_FA map thresholded at 3000)  

# **II.	Prerequisites**
# **a)	Hardware:** 
•	Computer with Linux or Mac OS X.  
•	For Windows, a Linux virtual machine is needed, e.g. the NeuroDebian Virtual Machine (http://neuro.debian.net/vm.html)  
# **b)	Software:** 
•	Mandatory: (1) An installation of the FMRIB Software Library (FSL, https://www.fmrib.ox.ac.uk/fsl). FSL is a free tool. See license details at: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence. (2) An installation of Python and of the following python libraries: “matplotlib” and “dipy”.  
•	Optional: dcm2niix (https://github.com/rordenlab/dcm2niix). This tool is designed to convert neuroimaging data from the DICOM format to the NIfTI format. This software is open source. The bulk of the code is covered by the BSD license. Some units are either public domain or use the MIT license.  
# **III.	Usage**
# **a)	Input needed** 
•	Full path to the 4D DTI volumes (nii.gz format). 
•	Full path to the brain mask image (nii.gz format)  
•	Full path to the bval file  
•	Full path to the bvec file  
•	Full path for the output directory  
# **b)	Modifications that must be done in both scripts before launching** 
Open MAIN_script_FW.sh file in /scripts_FW_CONSORTIUM/ directory. Three paths at the top of the script should be modified including FWMRN_PATH, BINFSL_PATH and PYTHON_EXEC (see below).  
 
•	FWMRN_PATH: must be replaced by the path leading to the /scripts_FW_CONSORTIUM/ directory  
•	BINFSL_PATH: must be replaced by the path leading to the FSL main directory  
•	PYTHON_EXEC: must be replaced by the path of the python executable  
# **c)	Before running the script for the first time**
Make sure the script is executable by running the following commands:  
cd /scripts_FW_CONSORTIUM/  
chmod +x MAIN_script_FW.sh  
# **d)	How to run the script**
Example:  
MAIN_script_FW.sh /dirsubject1/data.nii.gz /dirsubject1/brain_mask.nii.gz /dirsubject1/file.bval /dirsubject1/file.bvec /myoutputdir/  
# **IV.	Pipeline**
During the process, the main script will:  
1)	Create a new directory for the output files. The directory is defined by the user (see section III. a.)  
2)	Compute the FW  and FW-corrected FA maps.  
3)	Compute transformation parameters (named nat2std_warp) from the subject’s DTI space into the FSL FA template space   
4)	Coregister FW and FW-corrected FA maps in the FSL FA template space  
5)	Compute the mean FW and FW-corrected FA metrics within white matter voxels and store measures in the “summary.txt” file.  

**Reference paper:**
Maillard P, Lu H, Arfanakis K, et al. Instrumental validation of free water, peak-width of skeletonized mean diffusivity, and white matter hyperintensities: MarkVCID neuroimaging kits. Alzheimer's & Dementia:
Diagnosis, Assessment & Disease Monitoring 2022;14:e12261.
