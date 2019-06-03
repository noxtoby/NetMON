def initialise(base_parser, subparsers):
  import argparse
  from mrtrix3 import app
  parser = subparsers.add_parser('gif2', author='Neil P. Oxtoby (noxtoby@gmail.com), via freesurfer.py by Robert E. Smith (robert.smith@florey.edu.au)', synopsis='Generate the 5TT image based on a GIF2 parcellation image', parents=[base_parser])
  parser.add_argument('input',  help='The input GIF2 (Neuromorphometrics) parcellation image (usually containing \'NeuroMorph_Parcellation\' in its name)')
  parser.add_argument('output', help='The output 5TT image')
  options = parser.add_argument_group('Options specific to the \'gif2\' algorithm')
  options.add_argument('-lut', help='Manually provide path to the lookup table on which the input parcellation image is based (e.g. GIF2ColourLUT.txt)')



def checkOutputPaths():
  pass



def getInputs():
  import os, shutil
  from mrtrix3 import app, path, run
  run.command('mrconvert ' + path.fromUser(app.args.input, True) + ' ' + path.toTemp('input.mif', True))
  if app.args.lut:
    run.function(shutil.copyfile, path.fromUser(app.args.lut, False), path.toTemp('LUT.txt', False))



def execute():
  import os, sys
  from mrtrix3 import app, path, run

  lut_input_path = 'LUT.txt'
  if not os.path.exists('LUT.txt'):
    gif2_home = os.environ.get('GIFDB_HOME', '')
    if not gif2_home:
      app.error('Environment variable GIFDB_HOME is not set. Please set this variable manually, or provide script with path to file GIF2ColourLUT.txt using -lut option.')
    lut_input_path = os.path.join(gif_home, 'GIF2ColourLUT.txt')
    if not os.path.isfile(lut_input_path):
      app.error('Could not find GIF2 lookup table file (expected location: ' + lut_input_path + '), and none provided using -lut')

  if app.args.sgm_amyg_hipp:
    lut_output_file_name = 'GIF22ACT_sgm_amyg_hipp.txt'
  else:
    lut_output_file_name = 'GIF22ACT.txt'
  lut_output_path = os.path.join(path.sharedDataPath(), path.scriptSubDirName(), lut_output_file_name);
  if not os.path.isfile(lut_output_path):
    app.error('Could not find lookup table file for converting GIF2 parcellation output to tissues (expected location: ' + lut_output_path + ')')

  # Initial conversion from GIF2 parcellation to five principal tissue types
  run.command('labelconvert input.mif ' + lut_input_path + ' ' + lut_output_path + ' indices.mif')

  # Use mrcrop to reduce file size
  if app.args.nocrop:
    image = 'indices.mif'
  else:
    image = 'indices_cropped.mif'
    run.command('mrthreshold indices.mif - -abs 0.5 | mrcrop indices.mif ' + image + ' -mask -')

  # Convert into the 5TT format for ACT
  run.command('mrcalc ' + image + ' 1 -eq cgm.mif')
  run.command('mrcalc ' + image + ' 2 -eq sgm.mif')
  run.command('mrcalc ' + image + ' 3 -eq  wm.mif')
  run.command('mrcalc ' + image + ' 4 -eq csf.mif')
  run.command('mrcalc ' + image + ' 5 -eq path.mif')

  run.command('mrcat cgm.mif sgm.mif wm.mif csf.mif path.mif - -axis 3 | mrconvert - result.mif -datatype float32')

# Created by Neil Oxtoby on 2018-03-02
# Modified from gif3.py in MRtrix3 v3.0