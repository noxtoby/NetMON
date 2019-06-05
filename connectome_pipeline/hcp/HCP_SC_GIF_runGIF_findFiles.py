#!/usr/bin/env python
# encoding: utf-8
"""
HCP_SC_GIF_runGIF_findFiles.py

Created by Neil Oxtoby in October 2016.
Copyright (c) 2016 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
Parcellate structural MRI in preparation for generating 
Structural Connectomes from selected Human Connectome Project data.

Option 1 (latest version of GIF, subject to Jorge Cardoso's tweaks):
  cluster:/share/apps/cmic/GIF/runGIF_v3_folder.sh <INSERT_FOLDER_HERE>
Option 2 (frozen/stable version of GIF, courtesy of Arman Eshaghi):
  cluster:/cluster/project0/MS_LATA/GIF_old/bin/seg_GIF

    Usage:
        ./HCP_SC_GIF_runGIF_findFiles.py

'''
import fnmatch
import os
import shutil
from subprocess import call

def traverseFolders(topFolder='/Users/noxtoby/Documents/HCP/SelectedBySara',filenameWildcard='T1w*dc_restore.nii.gz'):
    """Traverses the file hierarchy within topFolder looking for files named T1w*.nii.gz
    and returns a list containing the paths, including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, filenameWildcard):
            matches.append(os.path.join(root, filename))
    return matches

def traverseFoldersDone(topFolder='/Users/noxtoby/Documents/HCP/SelectedBySara',filenameWildcard='T1w*NeuroMorph_Parcellation.nii.gz'):
    """Traverses the file hierarchy within topFolder looking for files named
    T1w*NeuroMorph_Parcellation.nii.gz
    and returns a list containing the paths, including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, filenameWildcard):
            matches.append(os.path.join(root, filename))
    return matches

def qsubGIF(workingDir, textFileListOfFilenames):
    """Generates a "qsub file" to be submitted to the UCL CS cluster for processing a 
    text-file-list of structural MR images using the GIF algorithm.
    """
    
    with open(textFileListOfFilenames) as f:
        nImages = sum(1 for _ in f)
    #print 'Found %i lines of text\n' % nImages
    qsubGIFText = """#$ -S /bin/bash
#$ -cwd # execute job from the current working directory
#$ -V # export environment variables
# $ -j y # merge stdout with stderr
# $ -R y # reservation y/n
#$ -pe smp 4 # parallel environment
#$ -l h_rt=23:59:00
#$ -l h_vmem=1.9G
#$ -l tmem=1.9G
#$ -l tscratch=10G
#$ -t 1-%i
#$ -N GIFHCP
EXEC=seg_GIF
EXEC_LOC=/cluster/project0/MS_LATA/GIF_old/bin
GIFDB=/cluster/project0/MS_LATA/NewGIFDB/db.xml
DATA_PATH=%s
FILELIST=%s
FILE=$(awk "NR==$SGE_TASK_ID" $FILELIST)
RESULTS_PATH=${FILE%%.*}
FILE_NAME=${FILE##*/}
FILE_BASE=${FILE_NAME%%.*}

if [[ ! -d ${RESULTS_PATH} ]] ; then
mkdir -p ${RESULTS_PATH}
fi

DIR_TEMP=`mktemp -d ./seg_GIF.XXXXXXXXXXXXXXXXXXX`

if [ ! -e ${RESULTS_PATH}/${FILE_BASE}_NeuroMorph_Segmentation.nii.gz  ] || [ ! -e ${RESULTS_PATH}/${FILE_BASE}_NeuroMorph_Parcellation.nii.gz ]; then
${EXEC_LOC}/${EXEC} -omp  4 -in ${FILE} -out ${RESULTS_PATH} -db ${GIFDB} -cpp ${DIR_TEMP} -temper 0.05 -regNMI -v 2 > ${RESULTS_PATH}_log.$SGE_TASK_ID 2>&1
rm -rf ${DIR_TEMP}
fi
    """ % (nImages, workingDir, textFileListOfFilenames)
    qsubGIFFile = textFileListOfFilenames + '.sh'
    outfile = open(qsubGIFFile, 'w')
    print >> outfile, qsubGIFText
    outfile.close()
    return qsubGIFFile

def unzipFiles(selectedSubjects):
    for s in selectedSubjects:
        if os.path.isfile('%s_3T_Structural_preproc.zip' % (s)):
            # unzip ${SUBJ}_3T_Structural_preproc ${SUBJ}/T1w/T1w_acpc_dc_restore.nii.gz ${SUBJ}/T1w/aparc+aseg.nii.gz
            t = "unzip -n %s_3T_Structural_preproc %s/T1w/T1w_acpc_dc_restore.nii.gz %s/T1w/aparc+aseg.nii.gz" % (s,s,s)
            os.system(t)
        if os.path.isfile('%s_3T_Diffusion_preproc.zip' % (s)):
            #unzip ${SUBJ}_3T_Diffusion_preproc ${SUBJ}/T1w/Diffusion/bvals ${SUBJ}/T1w/Diffusion/bvecs ${SUBJ}/T1w/Diffusion/data.nii.gz
            d = "unzip -n %s_3T_Diffusion_preproc %s/T1w/Diffusion/data.nii.gz %s/T1w/Diffusion/bvals %s/T1w/Diffusion/bvecs %s/T1w/Diffusion/nodif_brain_mask.nii.gz" % (s,s,s,s,s)
            os.system(d)



def main():
    """
    1. Unzips the files if necessary
    2. Decompresses gzipped Nifti files
    3. Saves as list in a text file
    4. Prints the rsync command for syncing the files to the cluster
    """
    
    workingDir = '/Users/noxtoby/Documents/HCP'
    topFolder = workingDir + '/SelectedBySara'
    #=== Unzip files
    SubjectsSelectedBySara = [139839, 141119] #[107018, 111009, 112314, 114217, 118023, 124624, 128026, 129331, 129937, 130619, 131419, 134223, 100206, 100610, 102513, 110613, 112112, 114641, 121416, 118225, 131823, 134728, 139839, 141119]
    unzipFiles(SubjectsSelectedBySara)
    
    #=== Decompress nifti files
    filenameWildcard = 'T1w*dc_restore.nii.gz'
    structuralMRI = traverseFolders(topFolder,filenameWildcard)
    for sMRI in structuralMRI:
        #== Option 1: convert datatype to UINT16 big endian (same datatype as PPMI images)
        #sMRInew = sMRI.replace('.nii.gz','_uint16be.nii')
        #comnd = 'mrconvert -datatype uint16be %s %s' % (sMRI, sMRInew)
        #filenameWildcard = 'T1w*dc_restore_uint16be.nii'
        #== Option 2: simply decompress
        sMRInew = sMRI.replace('.nii.gz','.nii')
        if not os.path.isfile(sMRInew):
            comnd = 'mrconvert %s %s' % (sMRI, sMRInew)
            os.system(comnd)
        filenameWildcard = 'T1w*dc_restore.nii'
    
    #=== Traverse folders and list diffusion files to be transferred to the cluster
    diffusionMRI = traverseFolders(topFolder,'data.nii.gz')
    diffusionMRI += traverseFolders(topFolder,'bvals')
    diffusionMRI += traverseFolders(topFolder,'bvecs')
    diffusionMRI += traverseFolders(topFolder,'nodif_brain_mask.nii.gz')
    diffusionMRI_relativePath = [kp.replace(topFolder + '/','') for kp in diffusionMRI]
    
    #=== Traverse folders and list T1w files to be transferred to the cluster
    filenameWildcard = 'T1w*dc_restore.nii'
    structuralMRI = traverseFolders(topFolder,filenameWildcard)
    structuralMRI_relativePath = [kp.replace(topFolder + '/','') for kp in structuralMRI]
    structuralMRIDone = traverseFoldersDone(topFolder)
    #= Now change filename from GIF result to original file
    structuralMRIDone = ['/'.join(sMRI.split('/')[0:-1])+'.nii' for sMRI in structuralMRIDone]
    structuralMRIDone_relativePath = [kp.replace(topFolder + '/','') for kp in structuralMRIDone]
    structuralMRIToDo_relativePath = [sMRI for sMRI in structuralMRI_relativePath if sMRI not in structuralMRIDone_relativePath]
    N = min(100,len(structuralMRIToDo_relativePath))
    if N!=0:
        #=== Choose only first N files
        print "%i Nifti files remaining to be processed." % (len(structuralMRIToDo_relativePath))
        print "Selected %i" % (N)
        structuralMRIToDo_relativePath = structuralMRIToDo_relativePath[0:N]
        #== Output list of files for GIF processing in a text file
        filenameList = '%s/HCP_T1w_files.txt' % (workingDir)
        outfile = open(filenameList, 'w')
        print >> outfile, "\n".join(str(i) for i in structuralMRIToDo_relativePath)
        outfile.close()
        print 'Transmit T1w Nifti files to the cluster using this command:'
        print "rsync -avh --progress -e 'ssh -p 22' --files-from=" + filenameList + " " + topFolder + "/ cluster:/SAN/medic/Net_Mod_MS/HCPData/"
    else:
        print 'No T1w images to be processed - all done.'
    #qsubGIFFile = qsubGIF(workingDir,filenameList)
    N = min(300,len(diffusionMRI_relativePath))
    if N!=0:
        #== Output list of diffusion files to be transferred to the cluster
        filenameListD = '%s/HCP_Diffusion_files.txt' % (workingDir)
        outfile = open(filenameListD, 'w')
        print >> outfile, "\n".join(str(i) for i in diffusionMRI_relativePath)
        outfile.close()
        print 'Transmit diffusion Nifti files to the cluster using this command:'
        print "rsync -avh --progress -e 'ssh -p 22' --files-from=" + filenameListD + " " + topFolder + "/ cluster:/SAN/medic/Net_Mod_MS/HCPData/"
    else:
        print 'No diffusion images to be transferred - all done.'



if __name__ == '__main__':
    main()


