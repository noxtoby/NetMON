#!/usr/bin/env python
# encoding: utf-8
"""
HCP_SC_GIF.py

Created by Neil Oxtoby in November 2016.
Copyright (c) 2016 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
Script for automating the MRtrix tutorial on HCP structural connectomes 
using a GIF/Neuromorphometrics parcellation of select HCP data.

That is, this should produce a structural connectome (SC) per individual, based on 
their pre-processed T1 structural and diffusion images from the HCP database. 

Should be submitted as a job on the UCL CS cluster after first running HCP_SC_GIF_qsub.py, 
which generates a useful qsub file after first searching for the necessary 
Neuromorphometrics atlas parcellation produced by processing the T1w image using the 
Geodesic Information Flows (GIF) algorithm from NiftySeg.

    Usage:
        HCP_SC_GIF.py HCPIndividualDataParentFolder
    Outcome:
        Produces a bunch of intermediate mrtrix files (100GB+) and a CSV connectome

'''
import fnmatch
import subprocess
import argparse
import os

# def traverseFolders(topFolder='/SAN/medic/Net_Mod_MS/HCPData',filenameWildcard='T1w*dc_restore.nii'):
def traverseFolders(topFolder='/Volumes/oxPOND/data/HCP/SelectedBySara',filenameWildcard='T1w*dc_restore.nii'):
    """Traverses the file hierarchy within topFolder looking for files matching
    the filenameWildcard argument and returns a list containing the paths,
    including filename, for each file."""
    matches = []
    for root, dirnames, filenames in os.walk(topFolder):
        for filename in fnmatch.filter(filenames, filenameWildcard):
            matches.append(os.path.join(root, filename))
    return matches

def HCP_SC(topFolder, T1w, T1w_GIF, dMRI, DWIbvals, DWIbvecs, dMRI_mask, algorithm='gif'):
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
    
    if algorithm == 'gif':
        T1w_parc = T1w_GIF
        colourLUTloc = os.environ.get('GIFDB_HOME') # This environment variable is defined in my .bash_profile
        colourLUT = '%s/GIFColourLUT.txt' % colourLUTloc
        connectomeLUT = '%s/%s_default.txt' % (connectomeLUTloc,algorithm)
    else:
        print 'ERROR: this script was written specifically for the GIF+Neuromorphometrics parcellation.\n       To use a FreeSurfer parcellation (for example), please modify the script.\n'
        #T1w_parc = T1_FREESURFER
        colourLUTloc = os.environ.get('FREESURFER_HOME')
        colourLUT = '%s/FreeSurferColorLUT.txt' % colourLUTloc
        connectomeLUT = '%s/fs_default.txt' % connectomeLUTloc
    
    print '\n\nProcessing structural images: %s, %s\n' % (T1w,T1w_parc)
    #*****************************************%
    #****** Structural Image Processing ******%
    #** 1. Generate a five-tissue-type segmented image appropriate for Anatomically-Constrained Tractography.
    fiveTT = '%s/5TT_%s.mif' % (topFolder,algorithm)
    if runPrelims:
        try:
            fiveTT_check = subprocess.call(['5ttgen',algorithm,T1w_parc,fiveTT])
            #print 'Skipping 5ttgen'
        except subprocess.CalledProcessError as err:
            print '5ttgen failed :-('
            fiveTT_check = err.returncode
    else:
        print 'Skipping 5ttgen'
    
    #** 2. Collapse the multi-tissue image into a 3D greyscale image for visualisation.
    # If the tissue segmentation image contains clearly erroneous tissue labels, you can delineate them 
    # manually using mrview‘s ROI editor tool, then apply your corrections to the tissue data using the 5ttedit command.
    fiveTTvis = '%s/5TT_%s_vis.mif' % (topFolder,algorithm)
    if runPrelims:
        try:
            fiveTTvis_check = subprocess.call(['5tt2vis',fiveTT,fiveTTvis])
        except subprocess.CalledProcessError as err:
            print '5tt2vis failed :-('
            fiveTTvis_check = err.returncode
    else:
        print 'Skipping 5tt2vis'
    
    #** 3. Modify the integer values in the parcellated image, such that the numbers in the image no longer 
    #   correspond to entries in the colour lookup table, but rows and columns of the connectome.
    # - The configuration file (fs_default.txt/gif_default.txt) is also a handy text file that provides a 
    #   structure name for every row / column of the connectome matrix.
    # - NOTE: mrtrix v0.3.15 uses labelconvert, mrtrix v0.3.14 uses labelconfig
    nodes = '%s/nodes_%s.mif' % (topFolder,algorithm)
    if runPrelims:
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
        print 'Skipping %s generation' % nodes
    
    #** (optional) 4. Replace (GIF/FreeSurfer) estimates of sub-cortical grey matter structures with estimates from FSL’s FIRST tool:
    nodes_fixSGM = '%s/nodes_%s_fixSGM.mif' % (topFolder,algorithm)
    if runPrelims:
        try:
            print 'Skipping labelsgmfix'
            #labelsgmfix_check = subprocess.call(['labelsgmfix',nodes,T1w_parc,connectomeLUT,nodes_fixSGM,'-premasked'])
        except subprocess.CalledProcessError as err:
            print 'labelsgmfix failed :-('
            labelsgmfix_check = err.returncode
    else:
        print 'Skipping labelsgmfix'
    #** Choose to keep GIF's subcortical GM parcellation
    nodesmif = nodes
    #** Or not...
    # nodesmif = nodes_fixSGM
    
    print '\n\nProcessing diffusion image: %s\n' % (dMRI)
    #****************************************%
    #****** Diffusion Image Processing ******%
    #** 1. Convert the diffusion images into a non-compressed format (not strictly necessary, but will make subsequent processing faster), 
    #  embed the diffusion gradient encoding information within the image header, 
    #  re-arrange the data strides to make volume data contiguous in memory for each voxel, 
    #  and convert to floating-point representation (makes data access faster in subsequent commands).
    dwi = '%s/DWI.mif' % (topFolder)
    if runPrelims:
        try:
            mrconvert_dwi_check = subprocess.call(['mrconvert','-fslgrad',DWIbvecs,DWIbvals,'-datatype','float32','-stride','0,0,0,1',dMRI,dwi])
        except subprocess.CalledProcessError as err:
            print 'mrconvert failed on DWI :-('
            mrconvert_dwi_check = err.returncode
    else:
        print 'Skipping DWI.mif creation'
    #** 2. Generate a mean b=0 image (useful for visualisation).
    meanb0 = '%s/meanb0.mif' % (topFolder)
    if runPrelims:
        try: #dwiextract -bzero ../../DWI.mif - | mrmath - mean ../../meanb0.mif -axis 3
            dwiextract_check = subprocess.Popen(('dwiextract', '-bzero',dwi,'-'), stdout=subprocess.PIPE)
            output = subprocess.check_output(('mrmath','-','mean',meanb0,'-axis','3'), stdin=dwiextract_check.stdout)
            dwiextract_check.wait()
        except subprocess.CalledProcessError as err:
            print 'dwiextract failed to extract a mean b=0 image :-('
            dwiextract_check = err.returncode
    else:
        print 'Skipping creation of %s' % meanb0
    #** 3. Estimate the multi-shell, multi-tissue response function.
    #  (((((( Need to update this part with single-shell version - start with Ashkan's code ))))))
    #  (((((( Take care! Order matters: WM, GM, CSF ))))))
    rfvoxels = '%s/RF_voxels_%s.mif' % (topFolder,algorithm)
    dwi2response_algorithm = 'msmt_5tt'
    rfwm = '%s/RF_WM.txt' % topFolder
    rfgm = '%s/RF_GM.txt' % topFolder
    rfcsf = '%s/RF_CSF.txt' % topFolder
    wmfods = '%s/WM_FODs.mif' % topFolder
    gm = '%s/GM.mif' % topFolder
    csf = '%s/CSF.mif' % topFolder
    tissueRGB = '%s/tissueRGB_%s.mif' % (topFolder,algorithm)
    if runPrelims:
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
            rfwm = '%s/RF_WM_%s.txt' % (topFolder,dwi2response_algorithm)
            rfgm = '%s/RF_GM_%s.txt' % (topFolder,dwi2response_algorithm)
            rfcsf = '%s/RF_CSF_%s.txt' % (topFolder,dwi2response_algorithm)
            wmfods = '%s/WM_FODs_%s.mif' % (topFolder,dwi2response_algorithm)
            gm = '%s/GM_%s.mif' % (topFolder,dwi2response_algorithm)
            csf = '%s/CSF_%s.mif' % (topFolder,dwi2response_algorithm)
            tissueRGB = '%s/tissueRGB_%s_%s.mif' % (topFolder,algorithm,dwi2response_algorithm)
            try:
                dwi2response_check = subprocess.call(['dwi2response',dwi2response_algorithm,dwi,rfwm,rfgm,rfcsf])
            except subprocess.CalledProcessError as err:
                print 'dwi2response failed to generate multi tissue response functions :-('
                dwi2response_check = err.returncode
            #** NOTE: How to check the appropriateness of response function voxel selections: 
            #mrview meanb0.mif -overlay.load RF_voxels_${algorithm}.mif -overlay.opacity 0.5 
        #** 4. Perform Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution:
        try:
            dwi2fod_check = subprocess.call(['dwi2fod','msmt_csd',dwi,rfwm,wmfods,rfgm,gm,rfcsf,csf,'-mask',dMRI_mask])
            #print 'Skipping dwi2fod'
        except subprocess.CalledProcessError as err:
            print 'dwi2fod failed to generate FODs :-('
            dwi2fod_check = err.returncode
        #=== Generate 4D RGB image of tissue types: CSF, WM, GM
        try:
            mrconvert_4D_check = subprocess.Popen(('mrconvert',wmfods,'-','-coord','3','0'), stdout=subprocess.PIPE)
            output = subprocess.check_output(('mrcat',csf,gm,'-',tissueRGB,'-axis','3'), stdin=mrconvert_4D_check.stdout)
            mrconvert_4D_check.wait()
            print 'Skipping %s creation' % tissueRGB
        except subprocess.CalledProcessError as err:
            print 'mrconvert | mrcat failed to concatenate tissue-types into a 4D, three-volumes image :-('
            mrconvert_4D_check = err.returncode
        #= This generates a 4D image with 3 volumes, corresponding to the tissue densities of CSF, GM and WM, 
        #= which will then be displayed in mrview as an RGB image with CSF as red, GM as green and WM as blue 
        #= (as was presented in the MSMT CSD manuscript).
        #mrview tissueRGB_${algorithm}.mif -odf.load_sh WM_FODs.mif # visual check
    else:
        print 'Skipping dwi2response'
    
    if runTractography:
        print '\n\nPerforming tractography\n'
        connectomeOut = '%s/connectome_%s_%s.tsv' % (topFolder,algorithm,dwi2response_algorithm)
        downsampleWMFODs = True
        #***********************************%
        #****** Connectome Generation ******%
        #** 1. Generate the initial tractogram: - TAKES AGES (1% after a few hours)
        maxLengthTracks = '250'
        numTracks = '30M'
        numTracksSIFT = '10M'
        if dwi2response_algorithm == 'msmt_5tt':
            tck = '%s/%s_%s.tck' % (topFolder,numTracks,algorithm) 
            sift = '%s/%s_%s_SIFT.tck' % (topFolder,numTracksSIFT,algorithm)
        else :
            tck = '%s/%s_%s_%s.tck' % (topFolder,numTracks,algorithm,dwi2response_algorithm) 
            sift = '%s/%s_%s_%s_SIFT.tck' % (topFolder,numTracksSIFT,algorithm,dwi2response_algorithm)
        SIFTorTCK = sift # For the connectome
        if downsampleWMFODs:
            print 'Down-sampling WM FODs to save memory'
            wm = '%s/WM_FODs_downsampled.mif' % topFolder
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
            tcksift_check = subprocess.call(['tcksift','-act',fiveTT,'-term_number',numTracksSIFT,tck,wmfods,sift])
        except subprocess.CalledProcessError as err:
            print 'tcksift failed to SIFT tracks :-('
            tcksift_check = err.returncode
        #=== If your system does not have adequate RAM to perform this process, the first recommendation is to reduce the spatial resolution of the FOD image and provide this alternative FOD image to SIFT (this should have little influence on the outcome of the algorithm, but will greatly reduce memory consumption). This has been handled above by "if downsampleWMFODs..."
        #=== If this still does not adequately reduce RAM usage, you will need to reduce the number of input streamlines to a level where your processing hardware can successfully execute the tcksift command, e.g.: tckedit 100M.tck 50M.tck -number 50M
        
        #** 3. Map streamlines to the parcellated image to produce a connectome:
        nodesmifuint32 = str.replace(nodesmif,'.mif','_uint32.mif')
        try:
            mrconvertuint32_check = subprocess.call(['mrconvert','-datatype','uint32',nodesmif,nodesmifuint32])
        except subprocess.CalledProcessError as err:
            print 'mrconvert failed to convert %s to UINT32 :-(' % nodesmif
            mrconvertuint32_check = err.returncode
        try:
            tck2connectome_check = subprocess.call(['tck2connectome','-assignment_all_voxels',SIFTorTCK,nodesmifuint32,connectomeOut])
        except subprocess.CalledProcessError as err:
            print 'tck2connectome failed to create the connectome :-('
            tck2connectome_check = err.returncode
        #=== View tracks as streamlines
        #mrview ${T1w_parc} -tractography.load ${SIFTorTCK} -mode 3
        #=== View tractography and connectome
        #mrview nodes_${algorithm}_uint32.mif -connectome.init nodes_${algorithm}_uint32.mif -connectome.load connectome_${algorithm}.csv -tractography.load ${SIFTorTCK}
        return 'Your connectome has been created! See file %s' %(connectomeOut)
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
    
    #=== Scan the folder for T1 image ===#
    filenameWildcardT1 = 'T1w*dc_restore.nii'
    T1w = traverseFolders(topFolder,filenameWildcardT1)
    T1w_relativePath = [kp.replace(topFolder + '/','') for kp in T1w]
    #=== Scan the folder for parcellated T1 image ===#
    filenameWildcardParcellated = str.replace(filenameWildcardT1,'.nii','_NeuroMorph_Parcellation.nii.gz')
    T1w_GIF = traverseFolders(topFolder,filenameWildcardParcellated)
    T1w_GIF_relativePath = [kp.replace(topFolder + '/','') for kp in T1w_GIF]
    #=== Scan the folder for diffusion image ===#
    filenameWildcardDiffusion = 'data.nii.gz'
    dMRI = traverseFolders(topFolder,filenameWildcardDiffusion)
    dMRI_relativePath = [kp.replace(topFolder + '/','') for kp in dMRI]
    filenameWildcardDiffusion = 'bvals'
    DWIbvals = traverseFolders(topFolder,filenameWildcardDiffusion)
    filenameWildcardDiffusion = 'bvecs'
    DWIbvecs = traverseFolders(topFolder,filenameWildcardDiffusion)
    filenameWildcardDiffusion = 'nodif_brain_mask.nii.gz'
    dMRI_mask = traverseFolders(topFolder,filenameWildcardDiffusion)
    
    #=== Call HCP_SC()
    print '* * * Calling HCP_SC() * * *'
    try:
        algorithm = 'gif'
        #=== Extract the list entries as strings
        i = 0
        T1w = T1w.pop(i)
        print T1w_GIF
        T1w_GIF = T1w_GIF.pop(i)
        dMRI = dMRI.pop(i)
        DWIbvals = DWIbvals.pop(i)
        DWIbvecs = DWIbvecs.pop(i)
        dMRI_mask = dMRI_mask.pop(i)
        HCP_SC_GIF_out = HCP_SC(topFolder, T1w, T1w_GIF, dMRI, DWIbvals, DWIbvecs, dMRI_mask, algorithm)
        print HCP_SC_GIF_out
    except subprocess.CalledProcessError as err:
        print 'HCP_SC_GIF.py: HCP_SC() failed :-('
        HCP_SC_GIF_out = err.returncode
    

if __name__ == '__main__':
    main()
