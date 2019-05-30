#!/usr/bin/env python
# encoding: utf-8
"""
collateT1sADNI.py
Modified from collateT1sHCP.py
  which was modified from collateT1sPPMI.py

Created by Neil Oxtoby in August 2016.
Updated February 2017.
Copyright (c) 2016-2017 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
1. Find a list of structural images in a given folder
2. See if there exists a corresponding GIF results folder, containing a NeuroMorph_Parcellation file for each
3. For those without, resubmit to the cluster as an array job (/share/apps/cmic/GIF/runGIF_v3_folder.sh submits individual jobs, which the cluster managers dislike)

The list of images to be processed is written in a text file 
- one file name per line - to be used in a Sun Grid Engine array job that 
processes the MRI using a Geodesic Information Flows (GIF) pipeline 
(seg_GIF) on the CMIC cluster.

Option 1 (latest version of GIF, subject to Jorge Cardoso's tweaks):
  cluster:/share/apps/cmic/GIF/bin/seg_GIF.sh
Option 2 (frozen/stable version of GIF, courtesy of Arman Eshaghi):
  cluster:/cluster/project0/MS_LATA/GIF_old/bin/seg_GIF

	Usage:
		/path/to/collateT1sADNI.py /path/to/ADNIFolder

'''
import fnmatch
import os
import shutil
import glob
import datetime

def traverseFolders(topFolder='/SAN/medic/Net_Mod_MS/HCPData/ADNI2GIF3_Images', wildcard='*_T1.nii.gz', depth=1):
    """Traverses the file hierarchy within topFolder looking for files named *_T1.nii.gz 
    and returns a list containing the paths, including filename, for each file."""
    if depth==1:
        print 'traverseFolders: depth==1 => globbing'
        imageFiles = glob.glob('{0}/{1}'.format(topFolder,wildcard))
        return imageFiles
    elif depth==2:
        print 'traverseFolders: depth==2 => globbing'
        imageFiles = glob.glob('{0}/*T1*/{1}'.format(topFolder,wildcard))
        return imageFiles
    else:
        print 'traverseFolders: depth!=1 => os.walk'
        matches = []
        for root, dirnames, filenames in os.walk(topFolder):
            for filename in fnmatch.filter(filenames, wildcard):
                matches.append(os.path.join(root, filename))
        return matches
    



def qsubGIF(workingDir, textFileListOfFilenames, jobName):
    """This will eventually be used to automate qsub submission 
    of an array job to the CMIC cluster at UCL, for processing a 
    list of T1 images provided in a text file."""
    
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
#$ -t 1-%s
#$ -N %s
EXEC=seg_GIF
# EXEC_LOC=/cluster/project0/MS_LATA/GIF_old/bin
# GIFDB=/cluster/project0/MS_LATA/NewGIFDB/db.xml
EXEC_LOC=/share/apps/cmic/GIF/gif_3ab700e/bin
GIFDB=/share/apps/cmic/GIF/db/db.xml

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
    ${EXEC_LOC}/${EXEC} -in ${FILE} -out ${RESULTS_PATH} -db ${GIFDB} -cpp ${DIR_TEMP} -temper 0.05 -regNMI -v 1 -lncc_ker -4 -omp 4 -regBE 0.001 -regJL 0.00005 > ${RESULTS_PATH}_log.$SGE_TASK_ID 2>&1
    ### OLD: ${EXEC_LOC}/${EXEC} -omp 4 -in ${FILE} -out ${RESULTS_PATH} -db ${GIFDB} -cpp ${DIR_TEMP} -temper 0.05 -regNMI -v 2 > ${RESULTS_PATH}_log.$SGE_TASK_ID 2>&1
    rm -rf ${DIR_TEMP}
fi
    """ % (nImages, jobName, workingDir, textFileListOfFilenames)
    qsubGIFFile = textFileListOfFilenames + '.sh'
    outfile = open(qsubGIFFile, 'w')
    print >> outfile, qsubGIFText
    outfile.close()
    return qsubGIFFile


def main():
    """docstring for main
    
    """
    runDate = str(datetime.date.today()).replace('-','')
    
    workingDir = '/SAN/medic/Net_Mod_MS/simADNI'
    topFolder = os.path.join(workingDir,'ADNI/035_S_4414/MT1__GradWarp__N3m/2012-03-19_15_27_42.0/S144404')
    
    #tfpAD = 'AD'
    #topFolderAD = topFolder + '/' + tfpAD
    #
    #tfpCN = 'Control'
    #topFolderCN = topFolder + '/' + tfpCN
    #
    #tfpLMCI = 'LMCI'
    #topFolderLMCI = topFolder + '/' + tfpLMCI
    #
    #tfpEMCI = 'EMCI'
    #topFolderEMCI = topFolder + '/' + tfpEMCI
    #
    
    tf = topFolder #topFolderEMCI
    tfp = '' #tfpEMCI
    workingDir = tf
    JOBNAME = 'simADNIGIF3' + tfp
    
    #*** Raw images
    wildcardRaw = '*sim_20180212_*.nii.gz'
    depth = 1 # use glob in the top level directory only
    rawMRI = traverseFolders(tf,wildcardRaw,depth)
    rawMRI.sort()
    
    #*** Processed images
    wildcardProcessed = '*T1_NeuroMorph_Parcellation.nii.gz'
    depth = 2 # anything other than 1 uses os.walk down through the directory tree
    processedMRI = traverseFolders(workingDir,wildcardProcessed,depth) # find GIF output
    processedMRI.sort()
    processedMRI_resultsFolder = ['/'.join(pathToImage.split('/')[0:-1]) for pathToImage in processedMRI] # rename to original file for comparison
    processedMRI_rawImageFile = ['/'.join(pathToImage.split('/')[0:-1])+'.nii.gz' for pathToImage in processedMRI] # rename to original file for comparison
    
    processedSubjectsWildcard = ['_'.join(t1.split('_')[0:-1]) for t1 in processedMRI_rawImageFile]
    
    #*** Move processed files
    #for k in range(0,len(processedMRI_rawImageFile)):
    #    mv_cmd = 'mv {0}* {1}/Done/{2}/'.format(processedSubjectsWildcard[k],topFolder,tfp)
    #    tx = os.system(mv_cmd)
    
    # #*** Delete old folders
    # for k in range(0,len(rawMRI)):
    #     cmd_ = 'rm -r {0}'.format(rawMRI[k].replace('.nii.gz',''))
    #     #print cmd_
    #     tx = os.system(cmd_)
    
    #*** Yet to be processed
    rawMRIYetToBeProcessed = [pathToImage for pathToImage in rawMRI if pathToImage not in processedMRI_rawImageFile]
    rawMRIYetToBeProcessed_folders = [pathToImage for pathToImage in rawMRI if pathToImage not in processedMRI_rawImageFile]
    
    #*** Optionally reduce the number of jobs submitted (helpful in case of file storage limitations)
    N = len(rawMRIYetToBeProcessed)
    rawMRIYetToBeProcessed = rawMRIYetToBeProcessed[0:N]
    
    #*** Write list of image files for processing to text file
    filenameList = '{0}/{1}_filelist_{2}.txt'.format(workingDir,JOBNAME,runDate)
    outfile = open(filenameList, 'w')
    print >> outfile, "\n".join(str(i) for i in rawMRIYetToBeProcessed)
    outfile.close()
    
    #*** Generate qsub file
    qsubGIFFile = qsubGIF(workingDir,filenameList,JOBNAME)
    #*** Print instructions for submitting the GIF processing to the cluster
    print "\n * * * ADNI-{0}: found {1} raw images, {2} processed images, and prepared {3} images for GIF processing * * * \n".format(tfp,len(rawMRI),len(processedMRI),len(rawMRIYetToBeProcessed))
    print '         Submit your job using:\n   qsub ' + qsubGIFFile  #+ ' -o ' + qsubGIFFile + '.out' + ' -e ' + qsubGIFFile + '.err'
    
    #*** Write list of processed image files to text file
    filenameList = filenameList.replace('.txt','_processingComplete.txt')
    outfile = open(filenameList, 'w')
    print >> outfile, "\n".join(str(i) for i in processedMRI_rawImageFile)
    outfile.close()


if __name__ == '__main__':
    main()


