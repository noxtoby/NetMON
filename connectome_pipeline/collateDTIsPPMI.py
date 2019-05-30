#!/usr/bin/env python
# encoding: utf-8
"""
collateDTIsPPMI.py

Modified from collateT1sPPMI.py by Neil Oxtoby in February 2017.
Copyright (c) 2017 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
1. Find a list of diffusion tensor images in a given folder
2. See if there exists a corresponding MRtrix results folder, containing a "_preproc" file for each
3. For those without, resubmit to the cluster as an array job

The list of images to be processed is written in a text file 
- one file name per line - to be used in a Sun Grid Engine array job that 
preprocesses the DTI in preparation for use in a structural connectome pipeline.
Preprocessing to be done using my script (to be called SCpreDTI.py).

    Usage:
	~/scripts/collateDTIsPPMI.py

'''
import fnmatch
import os
import shutil
import glob
import datetime

def traversePPMIFolders(topFolder='/SAN/medic/Net_Mod_MS/PPMI', wildcard='PPMI_*.mif', depth=1):
    """Traverses the file hierarchy within topFolder looking for files named PPMI_*.mif 
    (or other supplied wildcard) and returns a list containing the paths, including filename, 
    for each file."""
    if depth==1:
        print 'traversePPMIFolders: depth==1 => globbing'
        imageFiles = glob.glob('{0}/{1}'.format(topFolder,wildcard))
        return imageFiles
    elif depth==2:
        print 'traversePPMIFolders: depth==2 => globbing'
        imageFiles = glob.glob('{0}/PPMI*/{1}'.format(topFolder,wildcard))
        return imageFiles
    else:
        print 'traversePPMIFolders: depth!=1 => os.walk'
        matches = []
        for root, dirnames, filenames in os.walk(topFolder):
            for filename in fnmatch.filter(filenames, wildcard):
                matches.append(os.path.join(root, filename))
        return matches
    



def qsubDTI(workingDir, textFileListOfFilenames):
    """For automating qsub submission of an array job to the CMIC cluster at UCL, 
    for processing a list of diffusion tensor images provided in a text file."""
    
    with open(textFileListOfFilenames) as f:
        nImages = sum(1 for _ in f)
    print 'Found %i lines of text\n' % nImages
    qsubDTIText = """#$ -S /bin/bash
#$ -l hostname=burns*             # required: MRtrix build node
#$ -l h_rt=7:59:00               # Request wall time
#$ -l h_vmem=3.8G,tmem=3.8G       # Request memory
#$ -t 1-%s                        # Array job
#$ -N DTIpreproc                  # $JOB_NAME
#$ -wd /SAN/medic/Net_Mod_MS/PPMI # Working directory
# $ -l tscratch=100G
# $ -V # export environment variables
#$ -cwd # execute job from the current working directory
# $ -j y # merge stdout with stderr
# $ -R y # reservation y/n

EXEC=SCpreDTI_preproc.py
EXEC_LOC=~/scripts

DATA_PATH=%s
FILELIST=%s
FILE_PATH=$(awk "NR==$SGE_TASK_ID" $FILELIST)
FILE_BASE=`basename ${FILE_PATH} _denoised.mif`
RESULTS_PATH=${DATA_PATH}/${FILE_BASE}

FILE_NAME=`basename ${FILE_PATH}`
RESULTS_PATH_RELATIVE=${FILE_BASE}

if [[ ! -d ${RESULTS_PATH} ]] ; then
mkdir -p ${RESULTS_PATH}
fi

if [ ! -e ${RESULTS_PATH}/${FILE_BASE}_preproc.mif  ]; then
#    ${EXEC_LOC}/${EXEC} ${FILE_PATH} ${RESULTS_PATH}
    ${EXEC_LOC}/${EXEC} ${FILE_NAME} ${RESULTS_PATH_RELATIVE}
fi
    """ % (nImages, workingDir, textFileListOfFilenames)
    qsubDTIFile = textFileListOfFilenames + '.sh'
    outfile = open(qsubDTIFile, 'w')
    print >> outfile, qsubDTIText
    outfile.close()
    return qsubDTIFile


def main():
    """
    Put it all together and hope for the best. ;-)
    """
    workingDir = '/SAN/medic/Net_Mod_MS/PPMI'
    topFolder = workingDir
    runDate = str(datetime.date.today()).replace('-','')
    
    # Raw images
    wildcardRaw = 'PPMI_*denoised.mif'
    depth = 1 # use glob in the top level directory only
    rawDTI = traversePPMIFolders(topFolder,wildcardRaw,depth)
    
    # Processed images
    wildcardProcessed = 'PPMI_*preproc.mif'
    depth = 2 # anything other than 1 uses os.walk down through the directory tree
    processedDTI = traversePPMIFolders(topFolder,wildcardProcessed,depth) # find output from SCpreDTI.py
    processedDTI_resultsFolder = [os.path.split(pathToImage)[0] for pathToImage in processedDTI] # rename to original file for comparison
    processedDTI_rawImageFile = [os.path.split(pathToImage)[0]+'_denoised.mif' for pathToImage in processedDTI] # rename to original file for comparison
    
    # Move processed files
    print 'Moving processed files to PPMI_Done/'
    for k in range(0,len(processedDTI_rawImageFile)):
        tx = os.system('mv {0} /SAN/medic/Net_Mod_MS/PPMI_Done/'.format(processedDTI_rawImageFile[k]))
        tx = os.system('mv {0} /SAN/medic/Net_Mod_MS/PPMI_Done/'.format(processedDTI_resultsFolder[k]))
    
    # Yet to be processed
    rawDTIYetToBeProcessed = [pathToImage for pathToImage in rawDTI if pathToImage not in processedDTI_rawImageFile]
    rawDTIYetToBeProcessed_folders = [pathToImage for pathToImage in rawDTI if pathToImage not in processedDTI_rawImageFile]
    
    # Optionally reduce the number of jobs submitted (helpful in case of file storage limitations)
    N = len(rawDTIYetToBeProcessed)
    rawDTIYetToBeProcessed = rawDTIYetToBeProcessed[0:N]
    
    # Write list of image files for processing to text file
    filenameList = '{0}/PPMI_filesFor_SCpreDTI_preproc_{1}.txt'.format(workingDir,runDate)
    outfile = open(filenameList, 'w')
    print >> outfile, "\n".join(str(i) for i in rawDTIYetToBeProcessed)
    outfile.close()
    
    # Generate qsub file
    qsubDTIFile = qsubDTI(workingDir,filenameList)
    # Print instructions for submitting the SC pre-processing to the cluster
    print "\n * * * Found {0} raw images, {1} processed images, and prepared {2} DTIs for preprocessing with SCpreDTI_preproc.py * * * \n".format(len(rawDTI),len(processedDTI),len(rawDTIYetToBeProcessed))
    print 'Submit your job using:\n  qsub {0}'.format(qsubDTIFile)
    # + ' -o ' + qsubGIFFile + '.out' + ' -e ' + qsubGIFFile + '.err'
    
    # Write list of processed image files to text file
    filenameList = '{0}/PPMI_alreadyProcessedBy_SCpreDTI_preproc_{1}.txt'.format(workingDir,runDate)
    outfile = open(filenameList, 'w')
    print >> outfile, "\n".join(str(i) for i in processedDTI_rawImageFile)
    outfile.close()


if __name__ == '__main__':
    main()


