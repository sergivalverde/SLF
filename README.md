SALEM-LF (SLF)
=============

SALEM-LF (SLF) is a new algorithm designed to refill WM lesions on MRI T1-w images.  For each slice composing the three-dimensional brain image, lesion voxel intensities are replaced by random intensities of a normal distribution generated from the mean white matter intensity of the current slice. The algorithm recomputes the mean white matter signal intensity at each two-dimensional slice with the aim of reproducing more precisely the signal variability between MRI slices, especially in low resolution images.

Currently, SLF is currently available to download both as a ready-to-use Matlab function and as a SPM8 and SPM12 library. The useful resources section contains all available software source-codes as zip files. 

Please check the Requirements tab before installing it.

Please browse the SALEM webpage for more information (http://eia.udg.edu/salem/slfToolbox).

## Installation as a SPM library

 Installing SLFtoolbox as a SPM library is also straightforward. Download the toolbox, and extract it into a folder. These are the required steps to run SLF inside the SPM bundle.
 
 *  On SPM, all external libraries are installed into the toolbox folder. This folder is located inside the main SPM directory. Just find where SPM8 was installed and open the toolbox folder.
 *  Download the toolbox, and extract it into the spm/toolbox folder. Open MATLAB as usual, and type spm in the command-prompt. This will open the main SPM8 menu. Hence, on the Toolbox list, select the SLFtoolbox.
 *  The SPM batch editor will be opened containing the SLFtoolbox program.

The program requires two input images:

* A T1-w image to refill. If more of one image is selected, the same order must be followed for all T1-w images and their corresponding brain and lesion masks.
* A binary lesion mask image. If more of one image is selected, these images must follow the same order entered for the T1-w images to refill.

* Optionally. A binary brainmask image. If more of one image is selected, these images must follow the same order entered for the T1-w images to refill.

* The program is able to do bias correction and skull-stripping, via SPM8/SPM12 internal procedures. If no binary brain mask is provided, the program will attempt to compute the binary brainmask before the lesion-filling process. 

The resulted lesion-filled image will be saved in the same directory of the input T1-w image. Optionally, also bias corrected images and binray brain mask can be saved.
