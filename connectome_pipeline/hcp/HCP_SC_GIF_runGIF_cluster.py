#!/usr/bin/env python
# encoding: utf-8
"""
HCP_SC_GIF_runGIF_cluster.py

Created by Neil Oxtoby in October 2016.
Copyright (c) 2016 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
Generate Structural Connectomes from selected Human Connectome Project data.

1. Unzip selected data from each subject
2. Generate a text file containing the list of T1w MR images to be processed 
using the Geodesic Information Flows (GIF) algorithm (seg_GIF) on the UCL CS 
cluster.

Option 1 (latest version of GIF, subject to Jorge Cardoso's tweaks):
  cluster:/share/apps/cmic/GIF/runGIF_v3_folder.sh <INSERT_FOLDER_HERE>
Option 2 (frozen/stable version of GIF, courtesy of Arman Eshaghi):
  cluster:/cluster/project0/MS_LATA/GIF_old/bin/seg_GIF

=== Next steps would be... ===
X. Perform the local processing in my HCP_SC_GIF.sh "pipeline", i.e., 
everything up to the track generation and SIFTing (too memory-intensive for 
my desktop).
Y. Perform the rest on the cluster: tckgen, tcksift, tck2connectome

	Usage:
      	   ./HCP_SC_GIF_runGIF_cluster.py

'''
import fnmatch
import os
import shutil
from subprocess import call

def traverseFolders(topFolder='/SAN/medic/Net_Mod_MS/HCPData',filenameWildcard='T1w*dc_restore.nii'): 
    """Traverses the file hierarchy within topFolder looking for files matching 
    the filenameWildcard argument and returns a list containing the paths, 
    including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, filenameWildcard):
            matches.append(os.path.join(root, filename))
    return matches

def traverseFoldersDone(topFolder='/SAN/medic/Net_Mod_MS/HCPData',filenameWildcard='T1w*dc_restore_NeuroMorph_Parcellation.nii.gz'):
    """Traverses the file hierarchy within topFolder looking for files named
    T1w*NeuroMorph_Parcellation.nii.gz (or other supplied filenameWildcard)
    and returns a list containing the paths, including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, filenameWildcard):
            matches.append(os.path.join(root, filename))
    return matches

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
#$ -N GIFHCP2
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
    """
    1. Identify yet-to-be-parcellated T1w images.
    2. Save as list.
    3. Generate qsub file and print instructions for SGE array job submission.
    """
    
    workingDir = '/SAN/medic/Net_Mod_MS/HCPData' #'/home/noxtoby/pond/HCP'
    topFolder = workingDir #+ '/SelectedBySara'
    
    #=== Traverse folders and list T1w files
    filenameWildcardT1 = 'T1w*dc_restore.nii'
    structuralMRI = traverseFolders(topFolder,filenameWildcardT1)
    structuralMRI_relativePath = [kp.replace(topFolder + '/','') for kp in structuralMRI]
    #=== Traverse folders and list T1w files
    filenameWildcardT1Done = str.replace(filenameWildcardT1,'.nii','_NeuroMorph_Parcellation.nii.gz')
    structuralMRIDone = traverseFoldersDone(topFolder,filenameWildcardT1Done)
    
    #=== Iterate through the sMRI and check if there is a corresponding GIF output file
    structuralMRIDoneMatcher = [os.path.split(sMRI)[0]+'.nii' for sMRI in structuralMRIDone]
    structuralMRIToDo = [sMRI for sMRI in structuralMRI if sMRI not in structuralMRIDoneMatcher]
    
    if len(structuralMRIToDo)!=0:
        #=== Choose only first N files
        print "%i Nifti files remaining to be processed." % (len(structuralMRIToDo))
        print structuralMRIToDo
        print structuralMRIDone
        N = min(150,len(structuralMRIToDo))
        print "Selected %i" % (N)
        structuralMRIToDo = structuralMRIToDo[0:N]
        structuralMRIToDo = structuralMRIToDo[0:N]
        
        #=== Write the list of T1 files to be processed using GIF
        filenameList = '%s/HCP_T1w_files_for_GIF.txt' % (workingDir)
        outfile = open(filenameList, 'w')
        print >> outfile, "\n".join(str(i) for i in structuralMRIToDo)
        outfile.close()
        qsubGIFFile = qsubGIF(workingDir,filenameList)
        print '=== 1. Submit your GIF jobs using qsub: '
        print 'qsub ' + qsubGIFFile + ' -o ' + qsubGIFFile + '.out' + ' -e ' + qsubGIFFile + '.err'
    else:
        print '=== You are done: no GIF jobs to submit!'
    
    
    #=== Do the same thing as above, but for tractography (a la the MRTrix HCP Structural Connectome tutorial)
    #== Parcellated T1s
    #== Preprocessed diffusion images
    #== Check for completed tractography jobs via the presence of SIFT.tck files
    #filenameWildcardTCK = '*SIFT.tck'
    #filenameList = '%s/HCP_SIFT_files.txt' % (workingDir)
    #outfile = open(filenameList, 'w')
    #print >> outfile, "\n".join(str(i) for i in structuralMRIToDo)
    #outfile.close()
    #
    #qsubHCPSCFile = qsubHCPSC(workingDir,filenameList,filenameWildcardTCK)
    #print '=== 2. (in progress) Submit your HCP Structural Connectome jobs using qsub: '
    #print 'qsub ' + qsubHCPSCFile + ' -o ' + qsubHCPSCFile + '.out' + ' -e ' + qsubHCPSCFile + '.err'
    


if __name__ == '__main__':
    main()

