import fnmatch
import subprocess
import argparse
import os, sys
import errno
from shutil import copyfile
import glob
import pandas as pd
import numpy as np
import getopt
import xml.etree.ElementTree as ET
import time

from distutils.spawn import find_executable

############ Various functions ############
def make_sure_path_exists(path):
    """
    make_sure_path_exists(path) call mkdir and fails gracefully if the directory already exists.
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def checkForDependencies():
    """
    Checks for dependencies and returns nonzero values if any are missing.
    Dependencies:
      ANTS
      MRtrix
      FSL
      NiftyReg
    """
    dependencyCheckValue = 0
    #* Check for ANTS
    checkForANTS = os.system('N4BiasFieldCorrection --version')
    if checkForANTS!=0:
        print("ANTS does not seem to be installed. Have a look here:\n\thttps://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/")
        dependencyCheckValue = 1
    #* Check for MRtrix
    checkForMRtrix = os.system('mrinfo --version')
    if checkForMRtrix!=0:
        print('Error: MRtrix not installed!')
        dependencyCheckValue = 2
    #* Check for FSL
    checkForFSL = os.system('which fslinfo')
    if checkForFSL!=0:
        print('Error: FSL not installed!')
        dependencyCheckValue = 3
    #* Check for NiftyReg
    os.system('export PATH=${PATH}:/share/apps/cmic/niftk/niftk/bin')
    checkForNiftyReg = os.system('reg_aladin --version')
    if checkForNiftyReg!=0:
        print('Error: NiftyReg not installed!')
        dependencyCheckValue = 4
    #* Other dependencies would go here
    return dependencyCheckValue

def traverseFolders(topFolder='./', wildcard='*_T1.nii.gz', depth=1):
    """Traverses the file hierarchy within topFolder looking for files named *_T1.nii.gz 
    and returns a list containing the paths, including filename, for each file."""
    if depth==1:
        print 'traverseFolders: depth==1 => globbing'
        f = os.path.join(topFolder,wildcard)
        imageFiles = glob.glob(f) 
        return imageFiles
    else:
        print 'traverseFolders: depth!=1 => os.walk'
        matches = []
        for root, dirnames, filenames in os.walk(topFolder):
            for filename in fnmatch.filter(filenames, wildcard):
                matches.append(os.path.join(root, filename))
        return matches

def getPhaseEncoding(DWI,ADNI_boolean=True):
    """Returns 'AP' (the phase encoding for ADNI images). Eventually, code this to extract the PE direction from the input DWI."""
    if ADNI_boolean:
        pe_dir = 'AP'
    else: # Extract the phase encoding direction from the DWI???
        pe_dir = 'LR'
    return pe_dir
############ End: Various functions ############


############ ACT functions ############
def denoise(dMRI,dMRI_denoised=None):
    if dMRI_denoised is None:
        dMRI_denoised = dMRI.replace('.mif','_denoised.mif')
    if not os.path.isfile(dMRI_denoised):
        dwidenoise_check = subprocess.call(['dwidenoise',dMRI,dMRI_denoised])
    else:
        dwidenoise_check = 0
    #*** Parallelised batch version using MRtrix3's bash scripts (I haven't used these for a while: not guaranteed to work)
    #foreach dMRI/* : dwidenoise IN dMRI_denoised/NAME "&&" dwipreproc AP -rpe_none dMRI_denoised/NAME IN dMRI_preprocessed/NAME
    #foreach dMRI/* : dwidenoise IN dMRI_denoised/NAME
    return (dwidenoise_check,dMRI_denoised)

def preproc(dMRI_denoised,phase_encoding_dir='AP',rpe_none=False):
    #*** Cluster version: qsub dwipreproc_qsub.sh
    dMRI_preproc = dMRI_denoised.replace('.mif','_preprocessed.mif')
    if not os.path.isfile(dMRI_preproc):
        if rpe_none is False:
            dwipreproc_check = subprocess.call(['dwipreproc','-pe_dir',phase_encoding_dir,'-rpe_none',dMRI_denoised,dMRI_preproc])
        else
            print('Not implemented: need reverse-phase-encoded image to be passed as -rpe_pair argument.\nSee: http://mrtrix.readthedocs.io/en/latest/reference/scripts/dwipreproc.html for info.')
            dwipreproc_check = 1
            dMRI_preproc = None
    else:
        dwipreproc_check = 0
    return (dwipreproc_check,dMRI_preproc)

def biascorrect(dMRI_preproc,dMRI_mask=None):
    #*** Parallelised batch version using MRtrix3's bash scripts:
    #foreach dMRI_preprocessed/* : dwibiascorrect -ants -mask dMRI_masks/UNI_mask.mif dMRI_preprocessed/NAME IN dMRI_biascorrected/NAME
    dMRI_biascorrected = dMRI_preproc.replace('.mif','_biascorrected.mif')
    
    #*** Force the dwibiascorrect script to generate its own mask using dwi2mask
    dMRI_mask = None
    
    if not os.path.isfile(dMRI_biascorrected):
        if dMRI_mask is None:
            dwibiascorrect_check = subprocess.call(['dwibiascorrect','-ants',dMRI_preproc,dMRI_biascorrected])
        else:
            dwibiascorrect_check = subprocess.call(['dwibiascorrect','-ants','-mask',dMRI_mask,dMRI_preproc,dMRI_biascorrected])
    else:
        print("biascorrect(): bias-corrected file exists, skipping - {0}".format(dMRI_biascorrected))
        dwibiascorrect_check = 0
    return (dwibiascorrect_check,dMRI_biascorrected)

def reg_sMRI_to_dMRI(sMRI_nii,dMRI_biascorrected,sMRI_ext='nii.gz'):
    # NIFTI file is needed by NiftyReg
    dMRI_biascorrected_nifti = dMRI_biascorrected.replace('mif',sMRI_ext).replace('biascorrected','biascorrected_nifti') 
    sMRI_reg = sMRI_nii.replace('.nii','_reg_aladin.nii').replace('sMRI/','sMRI_reg/')
    reg_aladin_output_transform = sMRI_reg.replace('.{0}'.format(sMRI_ext),'_transform.txt')
    if not os.path.isfile(dMRI_biascorrected_nifti):
        mrconvert_check = subprocess.call(['mrconvert',dMRI_biascorrected,dMRI_biascorrected_nifti])
    if not os.path.isfile(sMRI_reg):
        reg_aladin_check = subprocess.call(['reg_aladin','-ref',dMRI_biascorrected_nifti,'-flo',sMRI_nii,'-rigOnly','-res',sMRI_reg,'-aff',reg_aladin_output_transform])
    else:
        reg_aladin_check = 0
    return (reg_aladin_check,sMRI_reg,reg_aladin_output_transform)

def fiveTT_GIF3(sMRI,algorithm = 'gif3'):
    fiveTT_check = -1
    fiveTT_ext = '_5TT_{0}.mif'.format(algorithm.upper())
    if algorithm == 'gif3':
        print('')
    elif algorithm == 'gif2':
        print('')
    else:
        print('fiveTT_GIF3(): Error! Supports only gif3 and gif2 algorithms. Could be modified for FreeSurfer fairly easily.')
    
    # Assume GIF naming defaults for output parcellation file
    (sMRI_path,sMRI_basename) = os.path.split(sMRI)
    file_base = sMRI_basename.split('.')[0]
    file_ext = '.'.join(sMRI.split(file_base)[1:])
    sMRI_GIF = os.path.join(sMRI_path,file_base,sMRI_basename.replace(file_ext,'_NeuroMorph_Parcellation.nii.gz'))
    sMRI_5TT = os.path.join(sMRI_path,file_base + fiveTT_ext)
    # Generate 5TT image using MRtrix
    print('\n\nProcessing structural image: {0}\n'.format(sMRI))
    print('           and parcellation: {0}\n'.format(sMRI_GIF))
    try:
        if not os.path.isfile(sMRI_5TT):
            fiveTT_check = subprocess.call(['5ttgen',algorithm,sMRI_GIF,sMRI_5TT])
    except subprocess.CalledProcessError as err:
        print('5ttgen failed :-(')
        fiveTT_check = err.returncode
    
    return (fiveTT_check,sMRI_5TT,sMRI_GIF)

def normalise():
    

def convert_sMRI_labels(sMRI_GIF,algorithm='gif3'):
    """
    Uses MRtrix's labelconvert to modify the integer values in the GIF-parcellated image, 
    such that the numbers in the image no longer correspond to entries in the colour lookup 
    table, but rows and columns of the connectome.
    
    The configuration file (fs_default.txt/gif3_default.txt/gif2_default.txt) is also a handy 
    text file containing a structure name for every row / column of the connectome matrix.
    """
    mrtrix_loc = os.path.abspath(os.path.join(find_executable('mrconvert'), os.pardir,os.pardir)) 
    connectome_LUT_loc = os.path.join(mrtrix_loc, 'share/mrtrix3/labelconvert')
    connectome_LUT = os.path.join(connectome_LUT_loc,algorithm + '_default.txt')
    if not os.path.isfile(connectome_LUT):
        print("convert_sMRI_labels(): Something ain't right! Connectome LUT not found: {0}".format(connectome_LUT))
    
    if algorithm == 'gif3':
        cluster = True
        if cluster:
            colour_LUT_loc = connectome_LUT_loc
        else:
            colour_LUT_loc = os.environ.get('GIFDB_HOME') # This environment variable is defined in my .bash_profile
        colour_LUT = os.path.join(colour_LUT_loc,algorithm.upper()+'ColourLUT.txt')
    else:
        print 'ERROR: this script was written specifically for the GIF+Neuromorphometrics parcellation.\n       To use a FreeSurfer parcellation (for example), please modify the script.\n'
        colour_LUT_loc = os.environ.get('FREESURFER_HOME')
        colour_LUT = os.path.join(colour_LUT_loc,'FreeSurferColorLUT.txt')
        connectome_LUT = os.path.join(connectome_LUT_loc,'fs_default.txt')
    
    nodes = os.path.dirname(sMRI_GIF) + '_nodes_' + algorithm + '.mif'
    if not os.path.isfile(nodes):
        try:
            labelconvert_check = subprocess.call(['labelconvert',sMRI_GIF,colour_LUT,connectome_LUT,nodes])
        except subprocess.CalledProcessError as err:
            print('convert_sMRI_labels(): labelconvert failed :-(')
                labelconvert_check = err.returncode
    
    return (labelconvert_check,nodes)

def preprocess_individual(sMRI,dMRI):
    """
    preprocess_individual(sMRI,dMRI)
    
    sMRI = raw structural MRI (T1)
    dMRI = raw diffusion MRI (DWI) in MRtrix Image Format (MIF), including bvals and bvecs in header
    
    Steps:
      dMRI
      1. denoise()     via dwidenoise
      2. preproc()     via dwipreproc
      3. biascorrect() via dwibiascorrect
    
    """
    algorithm = 'gif3'
    
    #*** 1. Denoise raw dMRI
    (dwidenoise_check,dMRI_denoised) = denoise(dMRI)
    #*** 2. Preprocess using MRtrix script: assumes no reverse phase encoding scan was performed
    (dwipreproc_check,dMRI_preproc) = preproc(dMRI_denoised)
    #*** 3. Bias-field correction using ANTS >>> requires ANTS to be installed
    checkForANTS = os.system('N4BiasFieldCorrection --version')
    if checkForANTS==0:
        (dwibiascorrect_check,dMRI_biascorrected) = biascorrect(dMRI_preproc)
    else:
        print("ANTS does not seem to be installed. Have a look here:\n\thttps://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/")
    
    #*** 4. Register sMRI to dMRI
    
    #*** 5. Run GIF3 on registered sMRI
    
    #*** 4. Generate 5TT image for ACT, assuming GIF parcellation exists
    (fiveTT_check,sMRI_5TT,sMRI_GIF) = fiveTT_GIF3(sMRI,algorithm)
    #*** 5. Use labelconvert to generate connectome-friendly numbering in the voxels
    (labelconvert_check,nodes) = convert_sMRI_labels(sMRI_GIF,algorithm)
    
    #*** 6. mtnormalise
    
    #*** 7. Generate 

def preprocess_dMRI(dMRI_raw,sMRI_raw,dMRI_mask,sMRI_outputFolder,dMRI_outputFolder):
    """
    preprocess_dMRI(path_to_raw_dMRI,path_to_raw_sMRI,path_to_dMRImask,outputFolder_sMRI,outputFolder_dMRI)
    performs some of the necessary preprocessing steps for group-wise analysis of DWI images,
    as prescribed by the good developers of mrtrix3, and customised for
    the ADNI/PPMI datasets:
    
    To enable robust quantitative comparisons of AFD across subjects (or AFD-derived
    quantities such as SIFT-filtered tractograms) there are three steps required:
    1. Bias field correction to eliminate low frequency intensity inhomogeneities across the image.
    2. Global intensity normalisation by normalising the median CSF or WM b=0 intensity across all
       subjects. This avoids the above noted issues with voxel-wise b=0 normalisation, such as CSF
       partial volume influencing restricted white matter DW signal (and therefore the AFD).
    3. Use the same single fibre response function in the spherical deconvolution for all subjects.
       This ensures differences in intra-axonal volume (and therefore DW signal) across subjects are
       detected as differences in the FOD amplitude (the AFD). See the AFD paper for more details.
    
    Usage:
    (dMRI_biascorrected,dMRI_mask,sMRI_reg,reg_aladin_output_transform) = preprocess_dMRI(path_to_raw_dMRI,path_to_raw_sMRI,path_to_dMRImask,outputFolder_sMRI,outputFolder_dMRI)
    """
    
    sMRI_ext = 'nii.gz'
    
    #*** 0. Preliminaries
    #* 0.1 Check if DICOM folder or not
    if os.path.isdir(dMRI_raw):
        dMRI_basename = dMRI_raw
        dMRI_mif = dMRI_basename + '.mif'
    else:
        dMRI_basename = os.path.splitext(os.path.basename(dMRI_raw).replace('.mif',''))[0]
        dMRI_mif = dMRI_raw
    if os.path.isdir(sMRI_raw):
        sMRI_basename = sMRI_raw
        sMRI_nii = '{0}.{1}'.format(sMRI_basename,sMRI_ext)
    else:
        sMRI_basename = os.path.splitext(os.path.basename(sMRI_raw).replace(sMRI_ext,''))[0]
        sMRI_nii = sMRI_raw
    #* 0.2 Filenames for later
    direc = os.path.split(os.path.split(os.path.abspath(dMRI_raw))[0])[0]
    dMRI_denoised = os.path.join(dMRI_outputFolder,'{0}_denoised.mif'.format(dMRI_basename))
    
    #*** 1. Convert raw sMRI (T1) to NII
    if not os.path.isfile(sMRI_nii):
        mrconvert_check = subprocess.call(['mrconvert',sMRI_raw,sMRI_nii])
    #*** 2. Denoise raw dMRI
    dwidenoise_check,dMRI_denoised = denoise(dMRI_mif,dMRI_denoised)
    #*** 3. Find phase encoding direction: PPMI MRI manual mentions L/R; ADNI is A/P
    isADNI = True
    phase_encoding_dir = getPhaseEncoding(dMRI_mif,isADNI)
    #*** 4. Preprocess using MRtrix script: assumes no reverse phase encoding scan was performed (error correction)
    (dwipreproc_check,dMRI_preproc) = preproc(dMRI_denoised,phase_encoding_dir)
    
    #*** I've removed this step as it occasionally made the dwibiascorrect script fail. dwibiascorrect estimates a mask itself using the same procedure, anyway. Weird
    # #*** 5. Estimate a brain mask from the dMRI, after preprocessing
    # if not os.path.isfile(dMRI_mask):
    #     dwi2mask_check = subprocess.call(['dwi2mask',dMRI_preproc,dMRI_mask])
    
    #*** 6. Bias-field correction using ANTS >>> requires ANTS to be installed
    checkForANTS = os.system('N4BiasFieldCorrection --version')
    if checkForANTS==0:
        (dwibiascorrect_check,dMRI_biascorrected) = biascorrect(dMRI_preproc,dMRI_mask)
    else:
        print("ANTS does not seem to be installed. Have a look here:\n\thttps://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/")
    #*** 7. Rigid-body registration of sMRI to dMRI - using NiftyReg (reg_aladin -rigOnly)
    (reg_aladin_check,sMRI_reg,reg_aladin_output_transform) = reg_sMRI_to_dMRI(sMRI_nii,dMRI_biascorrected)
    
    
    ########## TO DO: 
    # 1. remove T1-registration
    # 2. Add mtnormalise()
    # 3. CSD
    
    return (dMRI_biascorrected,dMRI_mask,sMRI_reg,reg_aladin_output_transform)


def preprocess_sMRI_GIF(sMRI_folder=None,wildcardRaw = '*.nii.gz'):
    """
    preprocess_sMRI_GIF(sMRI_folder,wildcardRaw = '*.nii.gz')
    
    Wrapper for GIF-processing of structural MRI for Anatomically Constrained Tractography using MRtrix3.
    
    GIF = Geodesic Information Flows for segmentation and parcellation.
    
    Developed on ADNI/PPMI data, so may require tweaking for your data.
    
    Usage:
      (filenameList_processed,filenameList_process,qsubGIFFile) = preprocess_sMRI_GIF(sMRI_folder,'*.nii.gz')
    
    """
    
    if sMRI_folder is None:
        sMRI_folder = os.getcwd()
    
    runDate = str(datetime.date.today()).replace('-','')
    
    gif_JOBNAME = 'gif3_' + sMRI_folder
    
    #*** Original images
    depth = 1 # use glob in the top level directory only
    rawMRI = traverseFolders(sMRI_folder,wildcardRaw,depth)
    rawMRI.sort()
    
    #*** Detect GIF-processed images
    wildcardProcessed = '*_NeuroMorph_Parcellation.nii.gz'
    depth = 2 # anything other than 1 uses os.walk down through the directory tree
    processedMRI = traverseFolders(workingDir,wildcardProcessed,depth) # find GIF output
    processedMRI.sort()
    #*** Match GIF parcellations to original images
    #* Requires that GIF puts results in a folder having the same basename as the original MRI (not guaranteed for future versions)
    processedMRI_resultsFolder = ['/'.join(pathToImage.split('/')[0:-1]) for pathToImage in processedMRI] 
    processedMRI_rawImageFile = ['/'.join(pathToImage.split('/')[0:-1])+'.nii.gz' for pathToImage in processedMRI] 
    processedSubjectsWildcard = ['_'.join(image_path.split('_')[0:-1]) for image_path in processedMRI_rawImageFile]
    if len(processedMRI_rawImageFile)>0:
        print('*** preprocess_sMRI_GIF(): GIF-processing complete for...')
        print('***   ...{0}'.format(processedMRI_rawImageFile))
    
    #*** Identify sMRI that are yet to be processed by GIF (write to text file for SGE qsub on a cluster)
    rawMRIYetToBeProcessed = [pathToImage for pathToImage in rawMRI if pathToImage not in processedMRI_rawImageFile]
    rawMRIYetToBeProcessed_folders = [pathToImage for pathToImage in rawMRI if pathToImage not in processedMRI_rawImageFile]
    if len(rawMRIYetToBeProcessed)>0:
        print('*** preprocess_sMRI_GIF(): GIF-processing NOT complete for...')
        print('  ...{0}'.format(rawMRIYetToBeProcessed))
    else:
        print('*** preprocess_sMRI_GIF(): GIF-processing complete for all images submitted')
    
    #*** Optionally reduce the number of jobs submitted (helpful in case of file storage limitations)
    N = len(rawMRIYetToBeProcessed)
    rawMRIYetToBeProcessed = rawMRIYetToBeProcessed[0:N]
    
    #*** Write list of image files for processing to text file
    filenameList_process = os.path.join(sMRI_folder,'{0}_filelist_{1}.txt'.format(gif_JOBNAME,runDate))
    outfile = open(filenameList_process, 'w')
    print >> outfile, "\n".join(str(i) for i in rawMRIYetToBeProcessed)
    outfile.close()
    
    #*** Generate qsub file
    print('*** preprocess_sMRI_GIF(): Generating qsub text file')
    qsubGIFFile = qsubGIF(workingDir,filenameList,JOBNAME)
    
    #*** Write list of processed image files to text file
    filenameList_processed = filenameList_process.replace('.txt','_processingComplete.txt')
    outfile = open(filenameList_processed, 'w')
    print >> outfile, "\n".join(str(i) for i in processedMRI_rawImageFile)
    outfile.close()
    
    #*** Print instructions for submitting the GIF processing to the cluster
    print("\n*** preprocess_sMRI_GIF():\n    Found {0} raw images,\n          {1} processed images (see {2}), and\n          {3} images were prepared for GIF processing (see {4})".format(len(rawMRI),len(processedMRI),filenameList_processed,len(rawMRIYetToBeProcessed),filenameList_process))
    print('    Submit your GIF jobs using:\n    qsub ' + qsubGIFFile) #+ ' -o ' + qsubGIFFile + '.out' + ' -e ' + qsubGIFFile + '.err'
    
    return (filenameList_processed,filenameList_process,qsubGIFFile)


def qsubGIF(workingDir, textFileListOfFilenames, GIFversion = 3, JOBNAME = 'bigGIF'):
    """This is used to automate qsub submission of an array job to the CS 
    cluster at UCL, for processing a list of structural MRI (T1) 
    provided in a text file."""
    
    if GIFversion==2:
        EXEC_LOC = '/cluster/project0/MS_LATA/GIF_old/bin'
        GIFDB = '/cluster/project0/MS_LATA/NewGIFDB/db.xml'
    elif GIFversion==3:
        EXEC_LOC = '/share/apps/cmic/GIF/gif_3ab700e/bin' #'/share/apps/cmic/GIF/bin'
        GIFDB = '/share/apps/cmic/GIF/db/db.xml'
    else:
        EXEC_LOC = 'ERROR: Unknown GIF version (use 2 or 3)'
        GIFDB = 'ERROR: Unknown GIF version (use 2 or 3)'
    
    with open(textFileListOfFilenames) as f:
        nImages = sum(1 for _ in f)
    print 'Found %i lines of text\n' % nImages
    qsubGIFText = """#$ -S /bin/bash
#$ -cwd # execute job from the current working directory
#$ -V # export environment variables
# $ -j y # merge stdout with stderr
# $ -R y # reservation y/n
#$ -pe smp 4 # parallel environment
#$ -l h_rt=23:59:00
#$ -l h_vmem=1.8G
#$ -l tmem=1.8G
#$ -t 1-%s
#$ -N %s
EXEC=seg_GIF
EXEC_LOC=%s
GIFDB=%s

DATA_PATH=%s
FILELIST=%s
FILE=$(awk "NR==$SGE_TASK_ID" $FILELIST)
# RESULTS_PATH=${FILE%%.*}
FILE_BASE=`basename ${FILE} .gz`
FILE_BASE=`basename ${FILE_BASE} .nii`
RESULTS_PATH=${FILE_BASE}

if [[ ! -d ${RESULTS_PATH} ]] ; then
mkdir -p ${RESULTS_PATH}
fi

DIR_TEMP=`mktemp -d ./seg_GIF.XXXXXXXXXXXXXXXXXXX`

if [ ! -e ${RESULTS_PATH}/${FILE_BASE}_NeuroMorph_Segmentation.nii.gz  ] || [ ! -e ${RESULTS_PATH}/${FILE_BASE}_NeuroMorph_Parcellation.nii.gz ]; then
    #*** V3: ${EXEC_LOC}/${EXEC} -in ${FILE} -out ${RESULTS_PATH} -db ${GIFDB} -cpp ${DIR_TEMP} -temper 0.05 -regNMI -v 1 -segPT 0.1 -lncc_ker -4 -omp 4 -regBE 0.001 -regJL 0.00005 > ${RESULTS_PATH}_log.$SGE_TASK_ID 2>&1
    #*** V2:
    ${EXEC_LOC}/${EXEC} -omp 4 -in ${FILE} -out ${RESULTS_PATH} -db ${GIFDB} -cpp ${DIR_TEMP} -temper 0.05 -regNMI -v 2 > ${RESULTS_PATH}_log.$SGE_TASK_ID 2>&1
    rm -rf ${DIR_TEMP}
fi
    """ % (nImages, JOBNAME, EXEC_LOC, GIFDB, workingDir, textFileListOfFilenames)
    qsubGIFFile = textFileListOfFilenames + '.sh'
    outfile = open(qsubGIFFile, 'w')
    print >> outfile, qsubGIFText
    outfile.close()
    return qsubGIFFile













#####################################################################################
cwd = '/Users/noxtoby/Documents/Research/NetMON/Analysis/Methods/ConnectomePipeline' #os.getcwd()
sMRI_reg_folder = os.path.join(cwd,'sMRI_reg')
dMRI_preprocessed_folder = os.path.join(cwd,'dMRI_preprocessed')
make_sure_path_exists(sMRI_reg_folder)
make_sure_path_exists(dMRI_preprocessed_folder)
#* Check for dependencies: ANTS, etc.
depCheckVal = checkForDependencies()
if depCheckVal!=0:
    print('ERROR: missing dependency {0}'.format(depCheckVal))
    print('see checkForDependencies() function for details')
    sys.exit(depCheckVal)

#* Input spreadsheet
inputCSV = 'inputSpreadsheet.csv'
df = pd.read_csv(os.path.join(cwd,inputCSV))

#* Match sMRI to dMRI using Subject ID
colsToKeep = ['Subject ID','Path']
df_sMRI = df[df['Modality']=='sMRI'].copy()
df_dMRI = df[df['Modality']=='dMRI'].copy()
df_sMRI = df_sMRI[colsToKeep].copy().rename(index=str, columns={"Path": "Path_sMRI"})
df_dMRI = df_dMRI[colsToKeep].copy().rename(index=str, columns={"Path": "Path_dMRI"})
df = pd.merge(df_sMRI,df_dMRI)

#*** TO DO: preprocess_individual()
#***        - ACT preprocessing on a per-individual basis
#*** 1. preprocess_dMRI: as below
#*** 2. register sMRI to dMRI_biascorrected
#*** 3. run GIF3 on sMRI_registered
#*** 4. ACT stuff: 5TT, nodes / WM_FODs (response functions & CSD)
#*** 5. mtnormalise
#*** 6. Tractography
#*** 7. Connectomes

###### Parallelise the pipeline for the cluster: submit a single job per individual: preprocess_dMRI() +  ######
#* Preprocess dMRI, using mostly MRtrix3: 
# - Convert NII to MIF, including b-values and b-vectors
# - preprocess_dMRI()
new_column_dMRI = 'Path_dMRI_biascorrected'
new_column_sMRI = 'Path_sMRI_registered'
for kf in range(0,df.shape[0]):
    dMRI = df['Path_dMRI'][kf]
    dMRI_mif = dMRI.replace('.nii.gz','.mif')
    #* Convert to MRtrix image format (MIF), including b-values and b-vectors
    if not os.path.isfile(dMRI_mif):
        bvec = dMRI.replace('.nii.gz','.bvec')
        bval = dMRI.replace('.nii.gz','.bval')
        mrconvert_check = subprocess.call(['mrconvert','-fslgrad',bvec,bval,dMRI,dMRI_mif])
    #*** I've removed the dwi2mask step because it occasionally makes the dwibiascorrect script fail
    dMRI_mask = dMRI_mif.replace('.mif','_mask.mif')
    # if not os.path.isfile(dMRI_mask):
    #     dwi2mask_check = subprocess.call(['dwi2mask',dMRI_mif,dMRI_mask])
    sMRI = df['Path_sMRI'][kf]
    print('************\n{0} of {1} - Preprocessing {2} and {3} using preprocess_dMRI()'.format(kf+1,df.shape[0],os.path.basename(dMRI_mif),os.path.basename(sMRI)))
    
    #************ preprocess dMRI ************
    (dMRI_biascorrected,dMRI_mask,sMRI_reg_aladin,reg_aladin_output_transform) = preprocess_dMRI(dMRI_mif,sMRI,dMRI_mask,sMRI_reg_folder,dMRI_preprocessed_folder)
    
    #*** Save new columns to dataframe
    df[new_column_dMRI][kf] = dMRI_biascorrected
    df[new_column_sMRI][kf] = sMRI_reg_aladin
    #--- Perhaps also write to input CSV to avoid the need to reprocess. Need to code it above: if [processedFileExists] then [skip]
    



######################################################################
######################################################################
### Code below is under development
######################################################################
######################################################################

#*** Loop through individuals for preprocessing
for kf in range(0,df.shape[0]):
    #*** GIF
    dMRI = df[new_column_dMRI][kf]
    sMRI = df[new_column_sMRI][kf]
    if os.path.isfile(dMRI) and os.path.isfile(sMRI):
        print("dMRI_biascorrected and sMRI_registered files exist. Submitting GIF processing job.")
        #qsub_GIF_individual(sMRI)
    else:
        print("Preprocessed dMRI and sMRI files not found. Skipping individual {0} of {1}: {2}".format(kf+1,df.shape[0],df['Subject ID'][kf]))
    #*** TO DO: dMRI intensity normalisation via mtnormalise
    dMRI_normalised = normalise_dMRI(dMRI)
    

#************ run GIF on sMRI_registered ************
#* Version 1: semiautomatic 
#*            Identify GIF-processed files and prepare for resubmission (qsub)
(filenameList_processed,filenameList_process,qsubGIFFile) = preprocess_sMRI_GIF(sMRI_reg_folder,'*.nii.gz')
#* Version 2: loop through dataframe (from inputSpreadsheet.csv)


###### dMRI: global intensity normalisation
### Option 1. mtnormalise
#* Note: requires a mask, so I probably have to run dwi2mask anyway - or edit dwiintensitynorm to call population_template without -mask_dir option
dMRI_normalised = normalise(dMRI = df[new_column_dMRI][kf])

### Option 2: dwiintensitynorm script (generate FA template using healthy controls)
#input_dir = 'dMRI_biascorrected'
#output_dir = 'dMRI_normalised'
#fa_template = os.path.join(os.getcwd(),'fa_template_group.mif')
#wm_mask = os.path.join(os.getcwd(),'wm_mask_group.mif')
#controls_ID = ['003_S_4081','003_S_4119','003_S_4288','003_S_4350','003_S_4441','003_S_4555','003_S_4644','003_S_4872','003_S_4900','007_S_4387','007_S_4488','007_S_4516','007_S_4620','007_S_4637','016_S_4097','016_S_4121','016_S_4688','016_S_4951','016_S_4952','021_S_4254','021_S_4421','029_S_4279','029_S_4290','029_S_4384','029_S_4585','029_S_4652','094_S_4234','094_S_4459','094_S_4503','094_S_4560','098_S_4003','098_S_4018','098_S_4050','098_S_4275','127_S_4148','127_S_4198','127_S_4604','127_S_4843','129_S_4369','129_S_4371','129_S_4396','129_S_4422']
#
##* Move masks into single folder
#mask_dir = 'dMRI_masks'
#
##* Move Nifti files from bias_corrected folder
#txt = 'mv {0}/*.nii.gz {1}/'.format(input_dir,input_dir + '_nii')
#os.system(txt)
#
##* dwiintensitynorm - takes ages!
#txt = 'dwiintensitynorm {0} {1} {2} {3} {4}'.format(input_dir,mask_dir,output_dir,fa_template,wm_mask)
#os.system(txt)










# #************ TCK2CONNECTOME ************#
# foreach /Volumes/maxtor201609/ADNI/201705_ADNI_SC/connectomes_gif3/*4M.tck : tck2connectome -assignment_all_voxels IN /Volumes/maxtor201609/ADNI/201705_ADNI_SC/connectomes_gif3/UNI_nodes_gif3.mif /Volumes/maxtor201609/ADNI/201705_ADNI_SC/connectomes_gif3/PRE_connectome_assignment_all_voxels.tsv
#
#
# #************ CONNETOME2TCK - for visualisation
# tck2connectome -assignment_all_voxels -out_assignments ${tckIn}_visualisation.tck ${tckIn} ${nodesMIF} ${connectomeOut_tsv}
# connectome2tck ${tckIn} ${assignments_txt} ${tckIn}_visualisation.tck

