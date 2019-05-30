import fnmatch
import subprocess
import argparse
import os, sys
import errno
from shutil import copyfile
import glob
import csv
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup

def traverseFolders(topFolder='.', wildcard='*.dcm', depth=1):
    """Traverses the file hierarchy within topFolder looking for files named *_5TT_gif3.mif 
    (or other supplied wildcard) and returns a list containing the paths, including filename, 
    for each file."""
    if depth==1:
        print('traverseFolders: depth==1 => globbing')
        imageFiles = glob.glob('{0}/{1}'.format(topFolder,wildcard))
        return imageFiles
    elif depth==2:
        print('traverseFolders: depth==2 => globbing')
        imageFiles = glob.glob('{0}/*/{1}'.format(topFolder,wildcard))
        return imageFiles
    else:
        print('traverseFolders: depth!=1 => os.walk')
        matches = []
        for root, dirnames, filenames in os.walk(topFolder):
            for filename in fnmatch.filter(filenames, wildcard):
                matches.append(os.path.join(root, filename))
        return matches

def writeBfiles(dMRI,bvecs,bvals):
    bvecsFile = dMRI.replace('.nii.gz','.bvec')
    bvalsFile = bvecsFile.replace('.bvec','.bval')
    with open(bvalsFile,'w') as f:
        f.write(' '.join(bvals))
        f.close()
    with open(bvecsFile,'w') as f:
        #* x/y/z-gradients
        xg = []
        yg = []
        zg = []
        for kb in range(0,len(bvecs)):
            xg.append(bvecs[kb][0])
            yg.append(bvecs[kb][1])
            zg.append(bvecs[kb][2])
        f.write(' '.join(xg))
        f.write('\n')
        f.write(' '.join(yg))
        f.write('\n')
        f.write(' '.join(zg))
        f.close()

#* Name of CSV file containing image file paths and IDs
csvFile_Input = 'inputSpreadsheet.csv'
createCSV = True

#* File locations
topFolder = '/Volumes/maxtor201609/ADNI/The-168-PostHoc'
sMRIFolder = os.path.normpath(os.path.join(topFolder,'sMRI'))
dMRIFolder = os.path.normpath(os.path.join(topFolder,'dMRI'))
dataFolder = os.path.join(topFolder,'ADNI')

#* Search for T1 files - this can be done manually if you already have nifti files
filenameWildcard = '*SPGR*.nii'
depth = 4
sMRI_list = traverseFolders(dataFolder,filenameWildcard,depth)
#* Structural MRI: custom renaming of ADNI files - your mileage may vary
sMRI_renamed = [os.path.basename(n) for n in sMRI_list]
sMRI_renamed = [bn.split('br_raw')[0] + 'S' + bn.split('_S')[-1] + '.gz' for bn in sMRI_renamed]
sMRI_renamed_path = [os.path.join(sMRIFolder,sMRI_renamed[k]) for k in range(0,len(sMRI_list))]
for k in range(0,len(sMRI_renamed_path)):
    if not os.path.exists(sMRI_renamed_path[0]):
        oot = subprocess.call(['mrconvert',sMRI_list[k],sMRI_renamed_path[k]])
    else:
        print('File exists, moving along: {0}'.format(sMRI_renamed_path[k]))

#* Search for DWI files - this can be done manually if you already have concatenated nifti files
filenameWildcard = '*DTI*.nii'
depth = 4
dMRI_list = traverseFolders(dataFolder,filenameWildcard,depth)
#* DWI: concatenate diffusion image series into a single file
dMRI_list_containingFolder = list(set([os.path.split(dMRI_list[k])[0] for k in range(0,len(dMRI_list))]))
#* DWI: custom renaming of ADNI files - your mileage may vary
dMRI_renamed = ['_'.join(['ADNI',fold.split('/')[6],'MR',fold.split('/')[7],os.path.basename(glob.glob(os.path.join(fold,'S*.xml'))[0]).split('.')[0]]) + '.nii.gz' for fold in dMRI_list_containingFolder]
dMRI_renamed_path = [os.path.join(dMRIFolder, dMRI_renamed[k]) for k in range(0,len(dMRI_renamed))]
for k in range(0,len(dMRI_renamed_path)):
    shell_cmd = ' '.join(['mrcat' , os.path.join(dMRI_list_containingFolder[k],'*.nii') , dMRI_renamed_path[k]])
    if not os.path.exists(dMRI_renamed_path[k]):
        oot = subprocess.call(shell_cmd, shell=True)
    else:
        print('File exists, moving along: {0}'.format(dMRI_renamed_path[k]))

#* Extract bvecs and bvals
filenameWildcard = '*.xml'
dMRI_xml_list = [] #traverseFolders(dataFolder,filenameWildcard,depth)
for k in range(0,len(dMRI_list_containingFolder)):
    dMRI_xml_list.append(glob.glob(os.path.join(dMRI_list_containingFolder[k],filenameWildcard)))
    bvecs = []
    bvals = []
    print('{0} of {1}: Parsing XML files in {2}'.format(k+1,len(dMRI_list_containingFolder),dMRI_list_containingFolder[k]))
    for j in range(0,len(dMRI_xml_list[k])):
        soup = BeautifulSoup(open(dMRI_xml_list[k][j]),'html.parser')
        if soup.metadata == None:
            bv = soup.diffusion.attrs['bvalue']
            bvals.append(bv)
            x = soup.diffusion.attrs['xgradient']
            y = soup.diffusion.attrs['ygradient']
            z = soup.diffusion.attrs['zgradient']
            bvecs.append([x,y,z])
    #* Write to files
    writeBfiles(dMRI_renamed_path[k],bvecs,bvals)

#* Save filepaths to CSV: match using Subject ID, then write to file
sMRI_SubjectID = [os.path.split(d)[-1][5:15] for d in sMRI_renamed_path]
sMRI_ImageID = [os.path.split(d)[-1].split('_I')[-1].split('.nii.gz')[0] for d in sMRI_renamed_path]
dMRI_SubjectID = [os.path.split(d)[-1][5:15] for d in dMRI_renamed_path]
dMRI_ImageID = [os.path.split(d)[-1].split('I')[-1].split('.nii.gz')[0] for d in dMRI_renamed_path]
if createCSV:
    print('Creating CSV file: {0}'.format(csvFile_Input))
    with open(os.path.join(topFolder,csvFile_Input), 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerow(['Subject ID','Image ID','Modality','Path'])
        for k in range(0,len(dMRI_renamed_path)):
            dwi = dMRI_renamed_path[k]
            dwi_iid = dMRI_ImageID[k]
            i = dMRI_SubjectID[k]
            for j in range(0,len(sMRI_SubjectID)):
                if sMRI_SubjectID[j]==i:
                    t1 = sMRI_renamed_path[j]
                    t1_iid = sMRI_ImageID[j]
                    break
                else:
                    t1 = None
                    t1_iid = None
            # Write DWI
            csvwriter.writerow([i,dwi_iid,'dMRI',dwi])
            # Write T1
            csvwriter.writerow([i,t1_iid,'sMRI',t1])
else:
    print('Not creating CSV file {0}'.format(csvFile_Input))

print('All done')
