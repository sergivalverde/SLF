function SLFtoolbox = tbx_cfg_SLFtoolbox
% SLF refills white matter lesions with the mean and half of the standard
%     deviation of the WM tissue of the input image.
%
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright 2014 Sergi Valverde/Arnau Oliver/ Xavier Lladó
%   $Revision: 1.2$     $Date: 06/05/14$
%       - 1.1: added support for ANALYZE and log
%	- 1.2: alternative processing to Image Processing toolbox  
%===============================================================================

%add to path the toolbox and used functions
if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','SLFtoolbox')); end
% image processing


% ---------------------------------------------------------------------
% WM LESION-FILLING on T1-w images
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% T1 image 
% ---------------------------------------------------------------------
data_T1 = cfg_files;
data_T1.tag  = 'data_T1';
data_T1.name = 'T1 volume to refill';
data_T1.help = {'Select T1-w scans for processing. When more than one scan is processed, this assumes that there is one scan for each subject. The same order of subjects must be followed between T1-w input scans and their corresponding brain masks, and lesion masks.'};
data_T1.filter = 'image';
data_T1.ufilter = '.*';
data_T1.num     = [1 Inf];

% ---------------------------------------------------------------------
% Lesion image
% ---------------------------------------------------------------------
data_LesionMask = cfg_files;
data_LesionMask.tag  = 'data_LesionMask';
data_LesionMask.name = 'Lesion mask binary volume';
data_LesionMask.help = {'Select binary lesion masks for processing. When more than one scan is processed, this assumes that there is one scan for each subject. Lesion masks must be in the same space than thes T1-w image. The same order of subjects must be followed between T1-w input scans and their corresponding brain masks, and lesion masks.'};
data_LesionMask.filter = 'image';
data_LesionMask.ufilter = '.*';
data_LesionMask.num     = [1 Inf];



% ---------------------------------------------------------------------
% biasreg Bias regularisation
% ---------------------------------------------------------------------
biasreg         = cfg_menu;
biasreg.tag     = 'biasreg';
biasreg.name    = 'Bias regularisation';
biasreg.help    = {
                   'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
                   ''
                   'An important issue relates to the distinction between intensity variations that arise because of bias artifact due to the physics of MR scanning, and those that arise due to different tissue properties.  The objective is to model the latter by different tissue classes, while modelling the former with a bias field. We know a priori that intensity variations due to MR physics tend to be spatially smooth, whereas those due to different tissue types tend to contain more high frequency information. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large values for the intensity non-uniformity parameters. This regularisation can be placed within a Bayesian context, whereby the penalty incurred is the negative logarithm of a prior probability for any particular pattern of non-uniformity.'
                   'Knowing what works best should be a matter of empirical exploration.  For example, if your data has very little intensity non-uniformity artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'
                   }';
biasreg.labels = {
                  'no regularisation (0)'
                  'extremely light regularisation (0.00001)'
                  'very light regularisation (0.0001) (SPM8 default)'
                  'light regularisation (0.001)'
                  'medium regularisation (0.01) '
                  'heavy regularisation (0.1)'
                  'very heavy regularisation (1)'
                  'extremely heavy regularisation (10)'
                  }';
biasreg.values = {
                  0
                  1e-05
                  0.0001
                  0.001
                  0.01
                  0.1
                  1
                  10
                  }';
biasreg.val    = {0};

% ---------------------------------------------------------------------
% biasfwhm Bias FWHM
% ---------------------------------------------------------------------
biasfwhm         = cfg_menu;
biasfwhm.tag     = 'biasfwhm';
biasfwhm.name    = 'Bias FWHM';
biasfwhm.help    = {'FWHM of Gaussian smoothness of bias. If your intensity non-uniformity is very smooth, then choose a large FWHM. This will prevent the algorithm from trying to model out intensity variation due to different tissue types. The model for intensity non-uniformity is one of i.i.d. Gaussian noise that has been smoothed by some amount, before taking the exponential. Note also that smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity non-uniformities.'};
biasfwhm.labels = {
                   '30mm cutoff'
                   '40mm cutoff'
                   '50mm cutoff'
                   '60mm cutoff (SPM8 default)'
                   '70mm cutoff'
                   '80mm cutoff'
                   '90mm cutoff'
                   '100mm cutoff'
                   '110mm cutoff'
                   '120mm cutoff'
                   '130mm cutoff'
                   '140mm cutoff'
                   '150mm cutoff'
                   'No correction'
                   }';
biasfwhm.values = {
                   30
                   40
                   50
                   60
                   70
                   80
                   90
                   100
                   110
                   120
                   130
                   140
                   150
                   Inf
                   }';
biasfwhm.val    = {Inf};

% ---------------------------------------------------------------------
% write Save Bias Corrected
% ---------------------------------------------------------------------
write         = cfg_menu;
write.tag     = 'write';
write.name    = 'Save Intermediate Bias Corrected';
write.help    = {'This is the option to save a bias corrected version of your images from this channel, or/and the estimated bias field. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.  The bias corrected version should have more uniform intensities within the different types of tissues.'};
write.labels = {
                'Save Nothing'
                'Save Bias Corrected'
                'Save Bias Field'
                'Save Field and Corrected'
                }';
write.values = {
                [0 0]
                [0 1]
                [1 0]
                [1 1]
                }';
write.val    = {[0 0]};



% ---------------------------------------------------------------------
% bias processing menu
% ---------------------------------------------------------------------
biasproc         = cfg_branch;
biasproc.tag     = 'biasproc';
biasproc.name    = 'Bias processing';
biasproc.val     = {biasreg biasfwhm write};
biasproc.help    = {'Bias processing options for input images which are not intensity corrected. By default, no regularisation is performed on input images.'
                    ''
                    'Bias processing is based on the same method included in the SPM8 tissue segmentation. If needed, please choose the bias field parameters according to the input image noise.'
                    ''
                    'Default values in SPM8 are:'
                    '- Bias regularisation:  very light regularisation 0.0001'
                    '- Bias FWHM:            60mm cutoff'
                    }';
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
% brain mask option
% ---------------------------------------------------------------------
brainmask         = cfg_menu;
brainmask.tag     = 'brainmask';
brainmask.name    = 'Perform skull-stripping';
brainmask.help    = {'Please indicate if skull stripping process have to be performed by the method. By default, the method expects an input binary brainmask'}';
    brainmask.labels = {
                  'do not perform skull-stripping (0)'
                  'perform internal skull-stripping (1)'
                  }';
brainmask.values = {
                  0
                  1                  
                  }';
brainmask.val    = {0};

% ---------------------------------------------------------------------
% thresholding segmentation output
% ---------------------------------------------------------------------
brainthresh         = cfg_entry;
brainthresh.tag     = 'brainthresh';
brainthresh.name    = 'Threshold to binarize the resulting brainmask';
brainthresh.help    = {'Value to threshold probabilistic SPM8 tissue segmentation output used to generate the resulting binary brainmask.'
                       'If the resulting skull-stripping mask tends to include non-brain parts, try to increase the value. By default, t=0.5 is considered sufficient for satisfactory results'
                     }';
brainthresh.num = [1 1];
brainthresh.val    = {0.5};


% ---------------------------------------------------------------------
% write Save Bias Corrected

% ---------------------------------------------------------------------
writebrainmask         = cfg_menu;
writebrainmask.tag     = 'writebrainmask';
writebrainmask.name    = 'Save computed brain mask';
writebrainmask.help    = {'This is the option to save the binary brainmask of the input T1-w to refill.'}';
writebrainmask.labels = {
                'Save Nothing'
                'Save Brain Mask'
                }';
writebrainmask.values = {
                0
                1
                }';
writebrainmask.val    = {0};



% ---------------------------------------------------------------------
% external brainmask
% ---------------------------------------------------------------------
data_BrainMask = cfg_files;
data_BrainMask.tag  = 'data_BrainMask';
data_BrainMask.name = 'Binary brainmask volume';
data_BrainMask.help = {'Indicate the corresponant binary brainmask associated with the T1-w input image.'
                        'When more than one scan is processed, this assumes that there is one binary brainmask for each subject. The same order of subjects must be followed between T1-w input scans and their corresponding brain masks.'
    };
data_BrainMask.filter = 'image';
data_BrainMask.ufilter = '.*';
data_BrainMask.num     = [0 Inf];
data_BrainMask.val = {''};

% ---------------------------------------------------------------------
% brain mask menu
% ---------------------------------------------------------------------
brainmaskproc         = cfg_branch;
brainmaskproc.tag     = 'brainmaskproc';
brainmaskproc.name    = 'Skull stripping';
brainmaskproc.val     = {brainmask writebrainmask brainthresh data_BrainMask};
brainmaskproc.help    = {'Skull stripping options for input images which do not include a binary brainmask. By default, no brainmask extraction is performed.'
                        ''
                        'Skull-stripping is based on the SPM8 tissue segmentation output. SPM8 probability masks are thresholded to generate a binary brainmask.'
                        'Please indicate if internal skull-stripping have to be performed, if you want to save the resulting binary mask, and the binary threshold.'
                        'If the resulting skull-stripping mask tends to include non-brain parts, try to increase the value. By default, t=0.5 is considered sufficient for satisfactory results'
                        }';
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
% Lesion Filling process
% ---------------------------------------------------------------------
lesionFilling         = cfg_exbranch;
lesionFilling.tag     = 'lesionFilling';
lesionFilling.name    = 'WM lesion-filling toolbox';
lesionFilling.val     = {data_T1 data_LesionMask biasproc brainmaskproc};
lesionFilling.help    = {'This toolbox can be useful to refill WM lesions with signal intensities similar to normal WM.'
                          ' Basically, for each slice composing the three-dimensional image, lesion voxel intensities are replaced by random intensities of a normal distribution generated from the mean normal WM intensity of the current slice.'
                          ' WM lesions are masked out from the T1-w image using the provided lesion mask, in order to avoid the influence of artificial lesions on tissue distributions. The resulting image is used to estimate the probability of each voxel to be classified as CSF, GM, and NAWM, by segmenting tissue with a Fuzzy-C-means approach. From the obtained tissue segmentation output,  we compute the three-dimensional NAWM mask from the image voxels with the highest probability to pertain to the WM cluster.'
                          ' Finally, the lesion-filling process is achieved as follows: for each slice composing the three-dimensional image, we compute the mean and standard deviation of the signal intensity of NAWM tissue. These values are used to generate a normal distribution with mean equal to the computed NAWM mean intensity and standard deviation equal to half of the computed NAWM standard deviation. Lesion voxels intensities from the current image slice are replaced by random values of the generated distribution. The procedure is repeated until all image slices are completed.'
                      
                          ''
                          ' The proposed algorithm requires at least two input images:'
                          '  - T1-w image'
                          '  - the corresponding  lesion mask.'
                          'If the corresponding brain mask is not provided, the method will compute it automatically based on the SPM8 segmentation output.'
                          'If input images are not intensity inhomogeneity corrected, the method can correct them using the internal SPM8 procedure.'
                          ''
     };
                          
%lesionFilling.prog    = @SLF_lesion_filling;
lesionFilling.prog    = @spm_local_preproc;
lesionFilling.vout    = @vout;

SLFtoolbox  = cfg_choice;
SLFtoolbox.name = 'SLFtoolbox';
SLFtoolbox.tag  = 'SLFtoolbox';

% menu d'entrada inicial 
SLFtoolbox.values = {lesionFilling};

end


% ------------------------------------------------------------------------
% lesion filling process is encapsulated inside the spm_local_preproc
% function
% ------------------------------------------------------------------------

function varargout = spm_local_preproc(job)
    
    % start log write for debug purposes
    global logfile;
    [pthT1, ~, ~] = spm_fileparts(job.data_T1{1});
    logfile = [pthT1,'/log_execution','.txt'];
    % start log file
    write_log(logfile, '', 'reset');
    write_log(logfile,'Starting Filling function SLF_spm8_lesion_filling_vout','');
    write_log(logfile,'-----------------------------------------------','');
    write_log(logfile, 'Chosen parameters:', '');
    write_log(logfile, ['T1: ', job.data_T1{1}], '');
    write_log(logfile, ['Lesion Mask: ', job.data_LesionMask{1}], '');
    write_log(logfile, ['biasreg: ', num2str(job.biasproc.biasreg)], '');
    write_log(logfile, ['biasfwhm: ', num2str(job.biasproc.biasfwhm)], '');
    write_log(logfile, ['write: ', num2str(job.biasproc.write)], '');
    write_log(logfile, ['brainmask: ', num2str(job.brainmaskproc.brainmask)], '');
    write_log(logfile, ['brainthresh: ', num2str(job.brainmaskproc.brainthresh)], '');
    write_log(logfile, ['writebrainmask: ', num2str(job.brainmaskproc.writebrainmask)], '');
    %write_log(logfile, ['data_BrainMask: ', job.brainmaskproc.data_BrainMask{1}], '');
    write_log(logfile,'-----------------------------------------------','');
    
    % run lesion filling process
    varargout{1} = SLF_spm8_lesion_filling_vout(job);
end

% ------------------------------------------------------------------------
% function vout:
%   generates the dependence files for the next process in the SPM batch
% 
% ------------------------------------------------------------------------

function dep = vout(job)
    % This depends on job contents, which may not be present when virtual
    % outputs are calculated.

    cdep = cfg_dep;

    for i=1:numel(job.data_T1)
        cdep(end +1) = cfg_dep;
        cdep(end).sname      = sprintf(['filled image ',num2str(i)]);
        cdep(end).src_output = substruct('.','filled','()',{i},'.','t1','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end

    % out dependencies
    dep = cdep;
end

