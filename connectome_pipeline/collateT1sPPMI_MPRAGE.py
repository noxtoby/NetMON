#!/usr/bin/env python
# encoding: utf-8
"""
collateT1sPPMI_MPRAGE.py

I may end up running the steps in this manually within a python command line.

Created by Neil Oxtoby in August 2016.
Copyright (c) 2016 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
Traverse imaging data folders and collate a list of the MRI filenames 
(requires nifti format).
Place this list into a text file - one file name per line - to be used 
in a Sun Grid Engine array job that processes the MRI using a 
Geodesic Information Flows (GIF) pipeline (seg_GIF) on the CMIC cluster.

Option 1 (latest version of GIF, subject to Jorge Cardoso's tweaks):
  cluster:/share/apps/cmic/GIF/runGIF_v3_folder.sh <INSERT_FOLDER_HERE>
Option 2 (frozen/stable version of GIF, courtesy of Arman Eshaghi):
  cluster:/cluster/project0/MS_LATA/GIF_old/bin/seg_GIF

	Preferred Usage:
		collateT1sPPMI_MPRAGE.py PPMIDownloadFolder 
	Actual Usage:
		collateT1sPPMI_MPRAGE.py

'''
import fnmatch
import os
import shutil

def traversePPMIFolders(topFolder='/home/noxtoby/pond/PPMI/MPRAGE'):
    """Traverses the file hierarchy within topFolder looking for files named PPMI*MPRAGE*.nii 
    and returns a list containing the paths, including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, 'PPMI*MPRAGE*.nii'):
            matches.append(os.path.join(root, filename))
    return matches

def traversePPMIFoldersDone(topFolder='/home/noxtoby/pond/PPMI/MPRAGE'):
    """Traverses the file hierarchy within topFolder looking for files named
    PPMI*MPRAGE*NeuroMorph_Parcellation.nii.gz
    and returns a list containing the paths, including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, 'PPMI*MPRAGE*NeuroMorph_Parcellation.nii.gz'):
            matches.append(os.path.join(root, filename))
    return matches

def qsubGIF(workingDir, textFileListOfFilenames):
    """This will eventually be used to automate qsub submission of an array job to the 
    CMIC cluster at UCL, for processing a list of T1 images provided in a text file."""
    
    with open(textFileListOfFilenames) as f:
        nImages = sum(1 for _ in f)
    print 'Found %i lines of text\n' % nImages
    qsubGIFText = """#$ -S /bin/bash
#$ -cwd # execute job from the current working directory
#$ -V # export environment variables
# $ -j y # merge stdout with stderr
# $ -R y # reservation y/n
#$ -pe smp 4 # parallel environment
#$ -l h_rt=23:59:0
#$ -l h_vmem=1.9G
#$ -l tmem=1.9G
# $ -l tscratch=10G
#$ -t 1-%i
#$ -N GIFPPMI
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


def main():
    """docstring for main
    
    """
    workingDir = '/home/noxtoby/pond/PPMI'
    topFolder = workingDir + '/MPRAGE'
    structuralMRI = traversePPMIFolders(topFolder)
    structuralMRI_relativePath = [kp.replace(topFolder + '/','') for kp in structuralMRI]
    structuralMRIDone = traversePPMIFoldersDone(topFolder)
    #print structuralMRIDone[0]
    
    #= Now change filename from GIF result to original file
    structuralMRIDone = ['/'.join(sMRI.split('/')[0:-1])+'.nii' for sMRI in structuralMRIDone]
    #print structuralMRIDone[0]
    
    structuralMRIDone_relativePath = [kp.replace(topFolder + '/','') for kp in structuralMRIDone]
    
    structuralMRIToDo = [sMRI for sMRI in structuralMRI if sMRI not in structuralMRIDone]
    structuralMRIToDo_relativePath = [sMRI for sMRI in structuralMRI_relativePath if sMRI not in structuralMRIDone_relativePath]
    #=== Choose only first 50 files
    #N = min(50,len(structuralMRIToDo))
    N = len(structuralMRIToDo)
    structuralMRIToDo = structuralMRIToDo[0:N]
    structuralMRIToDo_relativePath = structuralMRIToDo_relativePath[0:N]
    
    #print "Done[0] = " + structuralMRIDone_relativePath[0]
    #print "ToDo[0] = " + structuralMRIToDo_relativePath[0]
    
    filenameList = '%s/PPMI_MPRAGE_files.txt' % (workingDir)
    outfile = open(filenameList, 'w')
    print >> outfile, "\n".join(str(i) for i in structuralMRIToDo)
    outfile.close()
    
    qsubGIFFile = qsubGIF(workingDir,filenameList)
    
    print 'Submit your job using:\nqsub ' + qsubGIFFile + ' -o ' + qsubGIFFile + '.out' + ' -e ' + qsubGIFFile + '.err'

if __name__ == '__main__':
    main()


