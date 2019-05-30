#!/usr/bin/env python
# encoding: utf-8
"""
collateT1sPPMI_MPRAGE_GIFresults.py

Created by Neil Oxtoby in August 2016.
Copyright (c) 2016 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
Traverse imaging data folders and collate a list of the MRI filenames 
(assumes nifti format).
Place this list into a text file - one file name per line - to be used 
in a Sun Grid Engine array job that processes the MRI using a 
Geodesic Information Flows (GIF) pipeline (seg_GIF) on the CMIC cluster.

Option 1 (latest version of GIF, subject to Jorge Cardoso's tweaks):
  cluster:/share/apps/cmic/GIF/runGIF_v3_folder.sh <INSERT_FOLDER_HERE>
Option 2 (frozen/stable version of GIF, courtesy of Arman Eshaghi):
  cluster:/cluster/project0/MS_LATA/GIF_old/bin/seg_GIF

	Ideal usage (actual usage hard codes the argument):
		collateT1sPPMI_MPRAGE_GIFresults.py PPMIDownloadFolder 

'''
import fnmatch
import os
import shutil

def traversePPMIFolders(topFolder='/home/noxtoby/pond/PPMI/MPRAGE'):
    """Traverses the file hierarchy within topFolder looking for files named 
    PPMI*MPRAGE*NeuroMorph_Parcellation.nii.gz 
    and returns a list containing the paths, including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, 'PPMI*MPRAGE*NeuroMorph_Parcellation.nii.gz'):
            matches.append(os.path.join(root, filename))
    return matches

def qsubGIF(workingDir, textFileListOfFilenames):
    """This will eventually be used to automate qsub submission of an array job to the CMIC cluster at UCL, for 
    processing a list of T1 images provided in a text file."""
    
    with open(textFileListOfFilenames) as f:
        nImages = sum(1 for _ in f)
    print 'Found %i lines of text\n' % nImages
    qsubGIFText = """#$ -S /bin/bash
# $ -cwd # execute job from the current working directory
#$ -V # export environment variables
# $ -j y # merge stdout with stderr
# $ -R y # reservation y/n
#$ -pe smp 4 # parallel environment
#$ -l h_rt=7:59:0
#$ -l h_vmem=1.8G
#$ -l tmem=1.8G
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
    
    filenameList = '%s/PPMI_MPRAGE_GIFresults_files.txt' % (workingDir)
    outfile = open(filenameList, 'w')
    # imageFilenames = [sMRI[87:150] for sMRI in structuralMRI]
    imageFilenames = [sMRI.split('/')[-1].split('_NeuroMorph_Parcellation.nii.gz')[0] for sMRI in structuralMRI]
    print >> outfile, "\n".join(str(i) for i in imageFilenames)
    outfile.close()
    
    folderList = '%s/PPMI_MPRAGE_GIFresults_foldersToDelete.txt' % (workingDir)
    outfile = open(folderList, 'w')
    ##foldersToDeleteAfterSyncing = [sMRI[0:78] for sMRI in structuralMRI]
    #foldersToDeleteAfterSyncing = ['/'.join(sMRI.split('/')[0:-3]) for sMRI in structuralMRI]
    foldersToDeleteAfterSyncing = ['/'.join(sMRI.split('/')[6:-3]) for sMRI in structuralMRI]
    print >> outfile, "\n".join(str(i) for i in foldersToDeleteAfterSyncing)
    outfile.close()
    
    #qsubGIFFile = qsubGIF(workingDir,filenameList)
    #
    #print 'Submit your job using qsub ' + qsubGIFFile + ' -o ' + qsubGIFFile + '.out' + ' -e ' + qsubGIFFile + '.err' 
    print "=== 1. Download the GIF results using this command:"
    desktopMachine = "noxtoby2.cs.ucl.ac.uk"
    desktopFolder = "~/Documents/PPMI/Images/MPRAGE/"
    print "         rsync -avh --progress -e 'ssh -p 22' --include-from=" + folderList + " " + topFolder + "/ "+desktopMachine+":"+desktopFolder
    print "=== 2. Then delete the folders (use 'rm -r' with caution!!!):"
    print "         cd " + topFolder + "; xargs rm -r < " + folderList

if __name__ == '__main__':
    main()


