#!/usr/bin/env python
# encoding: utf-8
"""
HCP_SC_GIF3.py

Human Connectome Project Structural Connectome using parcellation by GIF version 3.

Created by Neil Oxtoby in February 2017.
Copyright (c) 2016 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
Script for automatic generation of structural connectomes 
using a GIF/Neuromorphometrics parcellation of selected, preprocessed structural and 
diffusion data from the Human Connectome PRoject (HCP).

GIF3 has additional white matter regions of interest to GIF2, and so required 
modifying the relevant connectome lookup tables in MRtrix.

    Usage:
        HCP_SC_GIF.py HCPIndividualDataParentFolder
    Outcome:
        Produces a bunch of intermediate mrtrix files (100GB+) and a CSV connectome

'''
import fnmatch
import subprocess
import argparse
import os
import glob

# # def traverseFolders(topFolder='/SAN/medic/Net_Mod_MS/HCPData',filenameWildcard='T1w*dc_restore.nii'):
# def traverseFolders(topFolder='/Volumes/datapond/datasets/HCP/SelectedBySara/GIF3connectomes',filenameWildcard='*T1w*dc_restore.nii.gz'):
#     """Traverses the file hierarchy within topFolder looking for files matching
#     the filenameWildcard argument and returns a list containing the paths,
#     including filename, for each file."""
#     matches = []
#     for root, dirnames, filenames in os.walk(topFolder):
#         for filename in fnmatch.filter(filenames, filenameWildcard):
#             matches.append(os.path.join(root, filename))
#     return matches

def traverseFolders(topFolder='/Volumes/datapond/datasets/HCP/SelectedBySara/GIF3connectomes', wildcard='*T1w*.nii.gz', depth=1):
    """Traverses the file hierarchy within topFolder looking for files named *T1w*.nii.gz
    and returns a list containing the paths, including filename, for each file."""
    if depth==1:
        print 'traverseFolders: depth==1 => globbing'
        imageFiles = glob.glob('{0}/{1}'.format(topFolder,wildcard))
        return imageFiles
    elif depth==2:
        print 'traverseFolders: depth==2 => globbing'
        imageFiles = glob.glob('{0}/*T1w*/{1}'.format(topFolder,wildcard))
        return imageFiles
    else:
        print 'traverseFolders: depth!=1 => os.walk'
        matches = []
        for root, dirnames, filenames in os.walk(topFolder):
            for filename in fnmatch.filter(filenames, wildcard):
                matches.append(os.path.join(root, filename))
        return matches


def HCP_SC(individualFilePrefix, T1w, T1w_parcellation, dMRI, DWIbvals, DWIbvecs, dMRI_mask, algorithm='gif3'):
    """docstring for HCP_SC"""
    runPrelims = True
    runTractography = True
    mrtrixLoc = '/Users/noxtoby/Documents/mrtrix3'
    connectomeLUTloc = '%s/src/connectome/tables/' % mrtrixLoc
    mrtrix_version = subprocess.check_output(['mrconvert','--version'])
    if '0.3.14' in mrtrix_version:
        mrtrix_version = '0.3.14'
    elif '0.3.15' in mrtrix_version:
        mrtrix_version = '0.3.15'
    else:
        mrtrix_version = '0.3.15'
    #print 'mrtrix version %s' % mrtrix_version
    #print algorithm

    T1w_parc = T1w_parcellation
    if algorithm == 'gif3':
        colourLUTloc = os.environ.get('GIFDB_HOME') # This environment variable is defined in my .bash_profile
        colourLUT = '%s/GIFColourLUTv3.txt' % colourLUTloc
        connectomeLUT = '%s/%s_default.txt' % (connectomeLUTloc,algorithm)
    elif algorithm == 'gif2':
        colourLUTloc = os.environ.get('GIFDB_HOME') # This environment variable is defined in my .bash_profile
        colourLUT = '%s/GIFColourLUTv2.txt' % colourLUTloc
        connectomeLUT = '%s/%s_default.txt' % (connectomeLUTloc,algorithm)
    elif algorithm == 'freesurfer':
        colourLUTloc = os.environ.get('FREESURFER_HOME')
        colourLUT = '%s/FreeSurferColorLUT.txt' % colourLUTloc
        connectomeLUT = '%s/fs_default.txt' % connectomeLUTloc
    else:
        print 'ERROR: this script was written specifically for GIF versions 3 or 2: the GIF+Neuromorphometrics parcellation.\n       To use a FreeSurfer parcellation (for example), please modify the script (Or use algorithm = "freesurfer" at your own risk).\n'
        print ''
        #T1w_parc = T1_FREESURFER

    fiveTT = '{0}_5TT_{1}.mif'.format(individualFilePrefix,algorithm)
    fiveTTvis = fiveTT.replace('.mif','_vis.mif')
    nodes = '{0}_nodes_{1}.mif'.format(individualFilePrefix,algorithm)
    nodes_fixSGM = nodes.replace('.mif','_fixSGM.mif')
    #** Choose to keep GIF's subcortical GM parcellation
    nodesmif = nodes
    #** Or not...
    # nodesmif = nodes_fixSGM
    nodesmifuint32 = nodesmif.replace('.mif','_uint32.mif')

    dwi = dMRI #'{0}_DWI.mif'.format(individualFilePrefix)
    meanb0 = '{0}_meanb0.mif'.format(individualFilePrefix)
    rfvoxels = '{0}_RF_voxels_{1}.mif'.format(individualFilePrefix,algorithm)
    dwi2response_algorithm = 'msmt_5tt'
    rfwm = '{0}_RF_WM_{1}.txt'.format(individualFilePrefix,algorithm)
    rfgm = '{0}_RF_GM_{1}.txt'.format(individualFilePrefix,algorithm)
    rfcsf = '{0}_RF_CSF_{1}.txt'.format(individualFilePrefix,algorithm)
    wmfods = '{0}_WM_FODs_{1}.mif'.format(individualFilePrefix,algorithm)
    gm = '{0}_GM_{1}.mif'.format(individualFilePrefix,algorithm)
    csf = '{0}_CSF_{1}.mif'.format(individualFilePrefix,algorithm)
    tissueRGB = '{0}_tissueRGB_{1}.mif'.format(individualFilePrefix,algorithm)

    downsampleWMFODs = True
    maxLengthTracks = '250'
    numTracks = '32M'
    numTracksSIFT = '8M'
    tck = '{0}_tckgen{1}_{2}_{3}.tck'.format(individualFilePrefix,numTracks,algorithm,dwi2response_algorithm) 
    sift = tck.replace('.tck','_tcksift{0}.tck'.format(numTracksSIFT))
    connectomeOut = sift.replace('.tck','_connectome.tsv')

    print '\n\nProcessing structural images: %s, %s\n' % (T1w,T1w_parc)
    #*****************************************%
    #****** Structural Image Processing ******%
    #** 1. Generate a five-tissue-type segmented image appropriate for Anatomically-Constrained Tractography.
    if runPrelims:
        if not os.path.isfile(fiveTT):
            try:
                fiveTT_check = subprocess.call(['5ttgen',algorithm,T1w_parc,fiveTT])
                #print 'Skipping 5ttgen'
            except subprocess.CalledProcessError as err:
                print '5ttgen failed :-('
                fiveTT_check = err.returncode
        else:
            print "{0} already exists!".format(fiveTT)
    else:
        print 'Skipping 5ttgen'

    #** 2. Collapse the multi-tissue image into a 3D greyscale image for visualisation.
    # If the tissue segmentation image contains clearly erroneous tissue labels, you can delineate them 
    # manually using mrview‘s ROI editor tool, then apply your corrections to the tissue data using the 5ttedit command.
    if runPrelims:
        if not os.path.isfile(fiveTTvis):
            try:
                fiveTTvis_check = subprocess.call(['5tt2vis',fiveTT,fiveTTvis])
            except subprocess.CalledProcessError as err:
                print '5tt2vis failed :-('
                fiveTTvis_check = err.returncode
        else:
            print "{0} already exists!".format(fiveTTvis)
    else:
        print 'Skipping 5tt2vis'

    #** 3. Modify the integer values in the parcellated image, such that the numbers in the image no longer 
    #   correspond to entries in the colour lookup table, but rows and columns of the connectome.
    # - The configuration file (fs_default.txt/gif_default.txt) is also a handy text file that provides a 
    #   structure name for every row / column of the connectome matrix.
    # - NOTE: mrtrix v0.3.15 uses labelconvert, mrtrix v0.3.14 uses labelconfig
    if runPrelims:
        if not os.path.isfile(nodes):
            if mrtrix_version == '0.3.15':
                try:
                    labelconvert_check = subprocess.call(['labelconvert',T1w_parc,colourLUT,connectomeLUT,nodes])
                except subprocess.CalledProcessError as err:
                    print 'labelconvert failed :-('
                    labelconvert_check = err.returncode
            elif mrtrix_version == '0.3.14':
                try:
                    labelconfig_check = subprocess.call(['labelconfig',T1w_parc,connectomeLUT,nodes,'-lut_'+algorithm,colourLUT])
                    #print 'Skipping nodes_gif.mif generation'
                except subprocess.CalledProcessError as err:
                    print 'labelconfig failed :-('
                    labelconfig_check = err.returncode
        else:
            print "{0} already exists!".format(nodes)
    else:
        print 'Skipping %s generation' % nodes

    #** (optional) 4. Replace (GIF/FreeSurfer) estimates of sub-cortical grey matter structures with estimates from FSL’s FIRST tool:
    if runPrelims:
        if not os.path.isfile(nodes_fixSGM):
            try:
                print 'Skipping labelsgmfix'
                #labelsgmfix_check = subprocess.call(['labelsgmfix',nodes,T1w_parc,connectomeLUT,nodes_fixSGM,'-premasked'])
            except subprocess.CalledProcessError as err:
                print 'labelsgmfix failed :-('
                labelsgmfix_check = err.returncode
        else:
            print "{0} already exists!".format(nodes_fixSGM)
    else:
        print 'Skipping labelsgmfix'

    print '\n\nProcessing diffusion image: %s\n' % (dMRI)
    #****************************************%
    #****** Diffusion Image Processing ******%
    #** 1. Convert the diffusion images into a non-compressed format (not strictly necessary, but will make subsequent processing faster), 
    #  embed the diffusion gradient encoding information within the image header, 
    #  re-arrange the data strides to make volume data contiguous in memory for each voxel, 
    #  and convert to floating-point representation (makes data access faster in subsequent commands).
    # if runPrelims:
    #     if not os.path.isfile(dwi):
    #         try:
    #             mrconvert_dwi_check = subprocess.call(['mrconvert','-fslgrad',DWIbvecs,DWIbvals,'-datatype','float32','-stride','0,0,0,1',dMRI,dwi])
    #         except subprocess.CalledProcessError as err:
    #             print 'mrconvert failed on DWI :-('
    #             mrconvert_dwi_check = err.returncode
    #     else:
    #         print "{0} already exists!".format(dwi)
    # else:
    #     print 'Skipping DWI.mif creation'

    #** 2. Generate a mean b=0 image (useful for visualisation).
    if runPrelims:
        if not os.path.isfile(meanb0):
            try: #dwiextract -bzero ../../DWI.mif - | mrmath - mean ../../meanb0.mif -axis 3
                dwiextract_check = subprocess.Popen(('dwiextract', '-bzero',dwi,'-'), stdout=subprocess.PIPE)
                output = subprocess.check_output(('mrmath','-','mean',meanb0,'-axis','3'), stdin=dwiextract_check.stdout)
                dwiextract_check.wait()
            except subprocess.CalledProcessError as err:
                print 'dwiextract failed to extract a mean b=0 image :-('
                dwiextract_check = err.returncode
        else:
            print "{0} already exists!".format(meanb0)
    else:
        print 'Skipping creation of %s' % meanb0
    #** 3. Estimate the multi-shell, multi-tissue response function.
    #  (((((( Need to update this part with single-shell version - start with Ashkan's code ))))))
    #  (((((( Take care! Order matters: WM, GM, CSF ))))))
    if runPrelims:
        if not os.path.isfile(rfvoxels):
            if dwi2response_algorithm == 'msmt_5tt':
                print 'Using the parametric MSMT response function algorithm: dwi2response msmt_5tt.'
                try:
                    dwi2response_check = subprocess.call(['dwi2response',dwi2response_algorithm,dwi,fiveTT,rfwm,rfgm,rfcsf,'-voxels',rfvoxels])
                except subprocess.CalledProcessError as err:
                    print 'dwi2response failed to generate multi tissue response functions :-('
                    dwi2response_check = err.returncode
            else:
                #** Nonparametric version that doesn't require T1, nor 5TT segmentation
                print 'Using the dhollander nonparametric response function algorithm: dwi2response msmt_5tt => dwi2fod msmt_csd.'
                print ' *** This is incomplete, i.e., don''t rely on it to work. ***'
                dwi2response_algorithm = 'dhollander'
                rfwm = '%s_RF_WM_%s.txt' % (individualFilePrefix,dwi2response_algorithm)
                rfgm = '%s_RF_GM_%s.txt' % (individualFilePrefix,dwi2response_algorithm)
                rfcsf = '%s_RF_CSF_%s.txt' % (individualFilePrefix,dwi2response_algorithm)
                wmfods = '%s_WM_FODs_%s.mif' % (individualFilePrefix,dwi2response_algorithm)
                gm = '%s_GM_%s.mif' % (individualFilePrefix,dwi2response_algorithm)
                csf = '%s_CSF_%s.mif' % (individualFilePrefix,dwi2response_algorithm)
                tissueRGB = '%s_tissueRGB_%s_%s.mif' % (individualFilePrefix,algorithm,dwi2response_algorithm)
                try:
                    dwi2response_check = subprocess.call(['dwi2response',dwi2response_algorithm,dwi,rfwm,rfgm,rfcsf])
                except subprocess.CalledProcessError as err:
                    print 'dwi2response failed to generate multi tissue response functions :-('
                    dwi2response_check = err.returncode
                #** NOTE: How to check the appropriateness of response function voxel selections: 
                #mrview meanb0.mif -overlay.load RF_voxels_${algorithm}.mif -overlay.opacity 0.5 
        else:
            print "{0} already exists! Skipping dwi2response.".format(rfvoxels)
    
        #** 4. Perform Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution:
        if not os.path.isfile(wmfods):
            try:
                dwi2fod_check = subprocess.call(['dwi2fod','msmt_csd',dwi,rfwm,wmfods,rfgm,gm,rfcsf,csf,'-mask',dMRI_mask])
                #print 'Skipping dwi2fod'
            except subprocess.CalledProcessError as err:
                print 'dwi2fod failed to generate FODs :-('
                dwi2fod_check = err.returncode
        else:
            print "{0} already exists! Skipping dwi2fod.".format(wmfods)
    
        #=== Generate 4D RGB image of tissue types: CSF, WM, GM
        if not os.path.isfile(tissueRGB):
            try:
                mrconvert_4D_check = subprocess.Popen(('mrconvert',wmfods,'-','-coord','3','0'), stdout=subprocess.PIPE)
                output = subprocess.check_output(('mrcat',csf,gm,'-',tissueRGB,'-axis','3'), stdin=mrconvert_4D_check.stdout)
                mrconvert_4D_check.wait()
                # print 'Skipping %s creation' % tissueRGB
            except subprocess.CalledProcessError as err:
                print 'mrconvert | mrcat failed to concatenate tissue-types into a 4D, three-volumes image :-('
                mrconvert_4D_check = err.returncode
        else:
            print "{0} already exists!".format(tissueRGB)
        #= This generates a 4D image with 3 volumes, corresponding to the tissue densities of CSF, GM and WM, 
        #= which will then be displayed in mrview as an RGB image with CSF as red, GM as green and WM as blue 
        #= (as was presented in the MSMT CSD manuscript).
        #mrview tissueRGB_${algorithm}.mif -odf.load_sh WM_FODs.mif # visual check
    else:
        print 'Skipping prelims for diffusion image'

    #***********************************%
    #****** Connectome Generation ******%
    #** 1. Generate the initial tractogram: - can take AGES (1% after a few hours)
    if runTractography:
        print '\n\nPerforming tractography\n'
        if downsampleWMFODs:
            print 'Down-sampling WM FODs to save memory when SIFTing'
            wm = wmfods.replace('.mif','_downsampled.mif')
            try:
                if not(os.path.isfile(wm)):
                    mrresize_check = subprocess.call(['mrresize',wmfods,wm,'-scale','0.5','-interp','sinc'])
                wmfods = wm
            except subprocess.CalledProcessError as err:
                print 'mrresize failed to downsample WM FODs :-('
                mrresize_check = err.returncode
        #print 'NEXT: run tckgen on cluster for this participant'
    
        try:
            if not(os.path.isfile(tck)):
                tckgen_check = subprocess.call(['tckgen','-act',fiveTT,'-backtrack','-crop_at_gmwmi','-seed_dynamic',wmfods,'-maxlength',maxLengthTracks,'-number',numTracks,'-cutoff','0.06',wmfods,tck])
            else:
                print "Skipping tckgen: %s file already exists" % (tck)
        except subprocess.CalledProcessError as err:
            print 'tckgen failed to complete tractography :-('
            tckgen_check = err.returncode
        #=== Explicitly setting the maximum length is highly recommended for HCP data, as the default heuristic - 100 times the voxel size - would result in a maximum length of 125mm, which would preclude the reconstruction of some longer pathways.
        #** 2. Apply the Spherical-deconvolution Informed Filtering of Tractograms (SIFT) method, which reduces the overall streamline count, but provides more biologically meaningful estimates of structural connection density:
        try:
            if not os.path.isfile(sift):
                tcksift_check = subprocess.call(['tcksift','-act',fiveTT,'-term_number',numTracksSIFT,tck,wmfods,sift])
            else:
                print "Skipping tcksift because output already exists: {0}".format(sift)
        except subprocess.CalledProcessError as err:
            print 'tcksift failed to SIFT tracks :-('
            tcksift_check = err.returncode
        #=== If your system does not have adequate RAM to perform this process, the first recommendation is to reduce the spatial resolution of the FOD image and provide this alternative FOD image to SIFT (this should have little influence on the outcome of the algorithm, but will greatly reduce memory consumption). This has been handled above by "if downsampleWMFODs..."
        #=== If this still does not adequately reduce RAM usage, you will need to reduce the number of input streamlines to a level where your processing hardware can successfully execute the tcksift command, e.g.: tckedit 100M.tck 50M.tck -number 50M
    
        #** 3. Map streamlines to the parcellated image to produce a connectome:
        try:
            if not os.path.isfile(nodesmifuint32):
                mrconvertuint32_check = subprocess.call(['mrconvert','-datatype','uint32',nodesmif,nodesmifuint32])
            else:
                print "Skipping: {0} exists already".format(nodesmifuint32)
        except subprocess.CalledProcessError as err:
            print 'mrconvert failed to convert %s to UINT32 :-(' % nodesmif
            mrconvertuint32_check = err.returncode
        try:
            if not os.path.isfile(connectomeOut):
                tck2connectome_check = subprocess.call(['tck2connectome','-assignment_all_voxels',sift,nodesmifuint32,connectomeOut])
            else:
                print "Connectome file already exists: {0}".format(connectomeOut)
        except subprocess.CalledProcessError as err:
            print 'tck2connectome failed to create the connectome :-('
            tck2connectome_check = err.returncode
        #=== View tracks as streamlines
        #mrview ${T1w_parc} -tractography.load ${SIFTorTCK} -mode 3
        #=== View tractography and connectome
        #mrview nodes_${algorithm}_uint32.mif -connectome.init nodes_${algorithm}_uint32.mif -connectome.load connectome_${algorithm}.csv -tractography.load ${SIFTorTCK}
        return 'Your connectome has been created! See file {0}'.format(connectomeOut)
    else:
        print 'Skipping tractography'



def main():
    """Scans the supplied HCP data folder checking for the required diffusion, structural, and parcellated structural images. 
    Then proceeds through the HCP SC steps."""
    #=== Parse inputs ===#
    parser = argparse.ArgumentParser(description='Process the script argument (should be a folder containing HCP data).')
    parser.add_argument('folder', metavar='F', type=str, help='top folder containing HCP data')
    args = parser.parse_args()
    #print args.folder
    topFolder = args.folder
    
    algorithm = 'gif2'
    
    #=== Scan the folder for T1 images ===#
    filenameWildcardT1 = '*T1w*dc_restore.nii.gz'
    t1s = traverseFolders(topFolder+'/t1',filenameWildcardT1)
    t1s_relativePath = [kp.replace(topFolder + '/','') for kp in t1s]
    
    if algorithm=='gif3' or algorithm=='gif2':
        t1s_parc = ['{0}/{1}'.format(t1.replace('.nii.gz','').replace('/t1/','/t1{0}/'.format(algorithm)),os.path.basename(t1).replace('.nii.gz','_NeuroMorph_Parcellation.nii.gz')) for t1 in t1s]
    elif algorithm=='freesurfer':
        t1s_parc = [t1.replace('_T1w_acpc_dc_restore.nii.gz','_aparc+aseg.nii.gz').replace('/t1/','/t1{0}/'.format(algorithm)) for t1 in t1s]
    #* Existing parcellation files only
    t1s_parc_ = [t1 for t1 in t1s_parc if os.path.isfile(t1)]
    
    subjects = [os.path.basename(t1).replace('_T1w_acpc_dc_restore.nii.gz','') for t1 in t1s_relativePath ]
    
    filenameWildcardDWI = "*DWI_normalised.mif"
    dwis = traverseFolders(topFolder+'/dwinormalised',filenameWildcardDWI)
    dwimasks = traverseFolders(topFolder+'/dwimasks',"*mask*")
    #dwibvals = traverseFolders(topFolder+'/dwibvals',filenameWildcardDWI)
    
    #=== Call HCP_SC()
    print '* * * Calling HCP_SC() * * *'
    for ki in range(0,len(subjects)): 
        subj = subjects[ki]
        indivFilePrefix = topFolder + "/connectomes_" + algorithm + "/" + subj
        print "============ Subject {0} ============".format(subj)
        t1 = t1s[ki]
        t1parc = t1s_parc[ki]
        dMRI = dwis[ki]
        bval = None #dwibvals[ki]
        bvec = None #dwibvecs[ki]
        dMRI_mask = dwimasks[ki]
        if not os.path.isfile(t1parc):
            print "  Parcellation not found.".format(subj)
            continue
        if not os.path.isfile(dMRI):
            print "  DWI file not found.".format(subj)
            continue
        print "  Calling HCP_SC()"
        try:
            HCP_SC_GIF_out = HCP_SC(indivFilePrefix , t1, t1parc, dMRI, bval, bvec, dMRI_mask, algorithm)
        except subprocess.CalledProcessError as err:
            print '  => HCP_SC() failed :-('
            HCP_SC_GIF_out = err.returncode
        print "  => HCP_SC returned: {0}".format(HCP_SC_GIF_out)
    

if __name__ == '__main__':
    main()
