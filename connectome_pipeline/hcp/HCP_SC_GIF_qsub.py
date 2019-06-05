#!/usr/bin/env python
# encoding: utf-8
"""
HCP_SC_GIF_qsub.py

Created by Neil Oxtoby in October 2016.
Copyright (c) 2016 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
Generates a qsub (Sun Grid Engine) file for processing Human Connectome Project (HCP) 
data to produce a structural connectome (SC) per individual, based on the 
Neuromorphometrics atlas parcellation produced using the Geodesic Information Flows 
(GIF) algorithm in NiftySeg.

    Usage:
        HCP_SC_GIF_qsub.py HCPDataFolder
    Outcome:
        writes HCP_SC_GIF_qsub.sh, which you need to submit using qsub:
        > qsub HCP_SC_GIF_qsub.sh -o HCP_SC_GIF.o -e HCP_SC_GIF.e

'''
import fnmatch
import subprocess
import argparse
import os

def traverseFolders(topFolder='/SAN/medic/Net_Mod_MS/HCPData',filenameWildcard='T1w*dc_restore.nii'):
    """Traverses the file hierarchy within topFolder looking for files matching
    the filenameWildcard argument and returns a list containing the paths,
    including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, filenameWildcard):
            matches.append(os.path.join(root, filename))
    return matches

def diff(first, second):
    """Returns items in the first list if not in the second list"""
    second = set(second)
    return [item for item in first if item not in second]

def intersec(first, second):
    """Returns items in the first list if also in the second list"""
    second = set(second)
    return [item for item in first if item in second]

def qsubGIF(workingDir, textFileListOfFilenames):
    """This can be used to automate qsub submission of an array job to the
    CS cluster at UCL, for processing a list of T1 images (structural MRI)
    provided in a text file."""

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
#$ -l h_vmem=1.9G
#$ -l tmem=1.9G
# $ -l tscratch=10G
#$ -t 1-%i
#$ -N GIF3
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

if [ ! -e ${RESULTS_PATH}/${FILE_BASE}_NeuroMorph_Segmentation.nii.gz  ] || [ ! -e ${RESULTS_PATH}/${FILE_BASE}_NeuroMorph_Parcellation.nii.gz ]\
; then
${EXEC_LOC}/${EXEC} -omp  4 -in ${FILE} -out ${RESULTS_PATH} -db ${GIFDB} -cpp ${DIR_TEMP} -temper 0.05 -regNMI -v 2 > ${RESULTS_PATH}_log.$SGE_\
TASK_ID 2>&1
rm -rf ${DIR_TEMP}
fi
    """ % (nImages, workingDir, textFileListOfFilenames)
    qsubGIFFile = textFileListOfFilenames + '.sh'
    outfile = open(qsubGIFFile, 'w')
    print >> outfile, qsubGIFText
    outfile.close()
    return qsubGIFFile

def qsubHCPSC(topFolder='/SAN/medic/Net_Mod_MS/HCPData',individualTopFolders=['1002006']):
    """
    (In progress) For submitting a Structural Connectome calculation to the 
    CS cluster at UCL. Currently focussed on HCP data (HCP_SC_GIF.py), but should be 
    adaptable.
    """
    qsubFiles = []
    for ki in range(0,len(individualTopFolders)):
        qsubFile = 'HCP_SC_GIF_%s.sh' % (individualTopFolders[ki])
        print 'Generating qsub script for %s: %s' % (individualTopFolders[ki],qsubFile)
        qsubFiles += qsubFile
        
        qsubFileText="""#!/bin/bash -l
# Batch script to run an array job under SGE.
# 0. MRtrix was built on this node, so we have to use this one
#$ -l hostname=burns*
# $ -l tscratch=100G
# 1. Force bash
#$ -S /bin/bash
# 2. Request hours of wallclock time (format hours:minutes:seconds)
#$ -l h_rt=119:59:00
# 3. Request RAM
#$ -l h_vmem=15G,tmem=7.5G
# 4. (not used here) Set up the job array, numbered 1-nJobs
#  -t 1-
# 5. Set the name of the job ($JOB_NAME in step 8).
#$  -N SC%i
# 7. Set the working directory.
#$ -wd %s/%s
# 8. Export my local variables to the compute node (may be unnecessary)
#$ -V
# 9. Run the application
../HCP_SC_GIF.py ../%s
""" % (individualTopFolders[ki],topFolder,individualTopFolders[ki],individualTopFolders[ki])
        
        outfile = open(qsubFile, 'w')
        print >> outfile, qsubFileText
        outfile.close()
    return qsubFiles

def main():
    """Scans the supplied HCP data folder checking for the required diffusion, structural, and parcellated structural images. 
    Then proceeds through the HCP SC steps."""
    #=== Parse inputs ===#
    parser = argparse.ArgumentParser(description='Process the script argument (should be a folder containing HCP data).')
    parser.add_argument('folder', metavar='F', type=str, help='top folder containing HCP data')
    args = parser.parse_args()
    topFolder = args.folder
    #print args.folder
    
    #=== Scan the folder for T1 images ===#
    filenameWildcardT1 = 'T1w*dc_restore.nii'
    structuralMRI = traverseFolders(topFolder,filenameWildcardT1)
    structuralMRI_relativePath = [kp.replace(topFolder + '/','') for kp in structuralMRI]
    #=== Scan the folder for parcellated T1 images ===#
    filenameWildcardParcellated = str.replace(filenameWildcardT1,'.nii','_NeuroMorph_Parcellation.nii.gz')
    structuralMRIParcellated = traverseFolders(topFolder,filenameWildcardParcellated)
    structuralMRIParcellated_relativePath = [kp.replace(topFolder + '/','') for kp in structuralMRIParcellated]
    #=== Scan the folder for diffusion images ===#
    filenameWildcardDiffusion = 'data.nii.gz'
    diffusionMRI = traverseFolders(topFolder,filenameWildcardDiffusion)
    diffusionMRI_relativePath = [kp.replace(topFolder + '/','') for kp in diffusionMRI]
    
    #=== Match structural to parcellation ===#
    structuralMRIMatcher = [os.path.split(sMRI)[0]+'.nii' for sMRI in structuralMRIParcellated]
    structuralMRIMatchedToParcellation = intersec(structuralMRI,structuralMRIMatcher)# [sMRI for sMRI in structuralMRI if sMRI in structuralMRIMatcher]
    structuralMRINotMatchedToParcellation = diff(structuralMRI,structuralMRIMatcher) #[sMRI for sMRI in structuralMRI if sMRI not in structuralMRIMatcher]
    # print '******** Structural MRI names to match with parcellations ********'
    # print structuralMRIMatcher
    # print ''
    # print '******** Structural MRI not matched to Parcellations ********'
    # print structuralMRINotMatchedToParcellation
    # print ''
    if len(structuralMRINotMatchedToParcellation)!=0:
        print 'Found some T1w images to be parcellated:'
        print structuralMRINotMatchedToParcellation
    else:
        print '=== All T1w images have been parcellated using GIF/Neuromorphometrics.'
    
    #=== Match structural to diffusion ===#
    diffusionMRIMatcher = ['/'.join(os.path.split(dMRI)[0].split('/')[0:-1])+'/T1w_acpc_dc_restore.nii' for dMRI in diffusionMRI]
    structuralMRIMatchedToDiffusion = intersec(structuralMRI,diffusionMRIMatcher) #[sMRI for sMRI in structuralMRI if sMRI in diffusionMRIMatcher]
    structuralMRINotMatchedToDiffusion = diff(structuralMRI,diffusionMRIMatcher) #[sMRI for sMRI in structuralMRI if sMRI not in diffusionMRIMatcher]
    # print '******** Structural MRI names to match with diffusion images ********'
    # print diffusionMRIMatcher
    # print ''
    # print '******** Structural MRI not matched to diffusion images ********'
    # print structuralMRINotMatchedToDiffusion
    # print ''
    
    #=== Find the individuals' top folder for running through the HCP SC tutorial pipeline (originally HCP_SC_GIF.py)
    structuralMRIMatchedToParcellationAndDiffusion = [sMRI for sMRI in structuralMRIMatchedToParcellation if sMRI in structuralMRIMatchedToDiffusion]
    HCPDataIndividualsTopFolder = [sMRI.replace(topFolder,'').split('/')[1] for sMRI in structuralMRIMatchedToParcellationAndDiffusion]
    print HCPDataIndividualsTopFolder
    #=== Identify files to be processed
    nJobsGIF = len(structuralMRINotMatchedToParcellation)
    nMissingDWI = len(structuralMRINotMatchedToDiffusion)
    if nMissingDWI!=0:
        print '=== However, I found missing diffusion images: %i missing to be precise' % nMissingDWI
    
    structuralMRIToDo_GIF = structuralMRINotMatchedToParcellation
    structuralMRIToDo_GIF_relativePath = [kp.replace(topFolder + '/','') for kp in structuralMRIToDo_GIF]
    
    #print structuralMRIMatchedToParcellationAndDiffusion
    
    if len(structuralMRIToDo_GIF)!=0:
        #=== Choose only first N files
        print "%i Nifti files remaining to be processed." % (len(structuralMRIToDo_GIF_relativePath))
        #print structuralMRIToDo_relativePath
        N = min(150,len(structuralMRIToDo_GIF_relativePath))
        print "Selected %i" % (N)
        structuralMRIToDo_GIF_relativePath = structuralMRIToDo_GIF_relativePath[0:N]
        structuralMRIToDo_GIF = structuralMRIToDo_GIF[0:N]
        #=== Write the list of T1 files to be processed using GIF
        filenameList = '%s/HCP_T1w_files_for_GIF.txt' % (topFolder)
        outfile = open(filenameList, 'w')
        print >> outfile, "\n".join(str(i) for i in structuralMRIToDo_GIF)
        outfile.close()
        qsubGIFFile = qsubGIF(topFolder,filenameList)
        print '=== 1. Submit your GIF jobs using qsub: '
        print 'qsub ' + qsubGIFFile + ' -o ' + qsubGIFFile + '.out' + ' -e ' + qsubGIFFile + '.err'
    else:
        print '' #'=== Hooray: no GIF jobs to submit!' 
        #=== qsub files for the cluster
        qsubFiles = qsubHCPSC(topFolder,HCPDataIndividualsTopFolder)
   

if __name__ == '__main__':
    main()
