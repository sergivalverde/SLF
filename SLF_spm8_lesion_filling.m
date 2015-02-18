function SLF_spm8_lesion_filling(job) 
%SLF refills white matter lesions with the mean and half of the standard
%deviation of the WM tissue of the input image.
%
% SLF(input_scan, output_scan, lesion_mask) saves a NIFTI image on the [output_scan]
%   location provided as input. The input_scans is already an skull-stripped image.
%   All input parameters should be existing paths to NIFTI images. For instance:
%
%       input_scan = '/dir/to/input/image/foo.nii'
%       refilled_scan= '/dir/to/output/image/output.nii'
%       lesion_mask= '/dir/to/lesionmask/image/foomask.nii'
%
%       SLF_lesion_filling(input_scan, refilled_scan, lesion_mask) 
% 
%   would save the refilled_image image on '/dir/to/output/image/output.nii'
%
% SLF(input_scan, output_scan, lesion_mask, brainmask) saves a NIFTI image on the [output_scan]
%   location provided as input. All input parameters should be existing paths
%   to NIFTI images. For instance:
%
%       input_scan = '/dir/to/input/image/foo.nii'
%       refilled_scan= '/dir/to/output/image/output.nii'
%       lesion_mask= '/dir/to/lesionmask/image/foomask.nii'
%       brain_mask =  '/dir/to/brainmask/image/brainmask.nii'
%
%       SLF_lesion_filling(input_scan, refilled_scan, lesion_mask, brain_mask) 
% 
%   would save the refilled_image image on '/dir/to/output/image/output.nii'
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright 2014 Sergi Valverde/Arnau Oliver/ Xavier Lladó
%   $Revision: 1.0$     $Date: 12/03/14$ 
%===============================================================================


% check input parameters. The number of T1_images should be equal to the
% number of brain masks and lesion masks. 
t1_elements = size(job.data_T1,1);
brainmask_elements = size(job.brainmaskproc.data_BrainMask,1);
lesion_elements = size(job.data_LesionMask,1);
if brainmask_elements > 0
    if t1_elements ~= brainmask_elements
        error('Looks like the number of T1 input images and brain masks is different');
    end
end

if t1_elements ~= lesion_elements
    error('Looks like the number of T1 input images and lesion masks is different');
end

% lesion filling is performed for all input images. We assume that T1 images and
% brain and lesion masks are sorted with the same order.

for im=1:t1_elements
    
    % ----------------------------------------------------------------- 
    % input data 
    % -----------------------------------------------------------------
    
    %T1 image
    [pthT1, nameT1, extT1] = spm_fileparts(job.data_T1{im});
    input_scan = fullfile(pthT1, [nameT1, extT1]);
    im_nii = load_untouch_nii(input_scan); input_image = im_nii.img;
    
    % out name
    refilled_scan = fullfile(pthT1, [nameT1,'_slf_filled', extT1]);
    disp(['Processing image: ', nameT1]);
    
    
    % lesion image
    [pth_lesionmask, name_lesionmask, ext_lesionmask] = spm_fileparts(job.data_LesionMask{im});  
    lesion_mask = fullfile(pth_lesionmask, [name_lesionmask, ext_lesionmask]);
    ls_nii = load_untouch_nii(lesion_mask); lesionmask = ls_nii.img;
    
    % brainmask image
    if job.brainmaskproc.brainmask == 0
        [pth_brainmask, name_brainmask, ext_brainmask] = spm_fileparts(job.brainmaskproc.data_BrainMask{im});
        brain_mask = fullfile(pth_brainmask, [name_brainmask, ext_brainmask]);
        bm_nii = load_untouch_nii(brain_mask);   brainmask = bm_nii.img;
        current_img = double(input_image) .* double(brainmask);
        im_nii.img = current_img;
        % save temporal file with skull stripping
        orig_name = nameT1;
        nameT1 = strcat(nameT1,'_brain_tmp');
        save_untouch_nii(im_nii, fullfile(pthT1, [nameT1,extT1]));
    end
    
   
    preprocess = job.brainmaskproc.brainmask + job.biasproc.biasreg > 0;
    
    if preprocess
        
        % ----------------------------------------------------------------- 
        % SPM 8 PREPROCESSING 
        % ----------------------------------------------------------------- 
        
        % if external brainmask, we load the tmp brainmasked image
        if job.brainmaskproc.brainmask == 0
            matlabbatch{1}.spm.tools.preproc8.channel.vols = {fullfile(pthT1, [nameT1, extT1])};
        else
            matlabbatch{1}.spm.tools.preproc8.channel.vols = {job.data_T1{im}};
        end
        matlabbatch{1}.spm.tools.preproc8.channel.biasreg = job.biasproc.biasreg;
        matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = job.biasproc.biasfwhm;
        matlabbatch{1}.spm.tools.preproc8.channel.write = job.biasproc.write;
        matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spm('dir') '/toolbox/Seg/TPM.nii,1']};
        matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
        matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spm('dir') '/toolbox/Seg/TPM.nii,2']};
        matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
        matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spm('dir') '/toolbox/Seg/TPM.nii,3']};
        matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];	
        matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spm('dir') '/toolbox/Seg/TPM.nii,4']};
        matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[spm('dir') '/toolbox/Seg/TPM.nii,5']};
        matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[spm('dir') '/toolbox/Seg/TPM.nii,6']};
        matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
        matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
        matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
        matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
        matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0]; 

        spm_jobman('run',matlabbatch);

    
        
        % if internal brainmask, we save the resulting brainmask
        if job.brainmaskproc.brainmask == 1
            c1 = load_untouch_nii([pthT1, '/c1', nameT1, extT1]);
            c2 = load_untouch_nii([pthT1, '/c2', nameT1, extT1]);
            c3 = load_untouch_nii([pthT1, '/c3', nameT1, extT1]);
            
            % normalize input images
            c1.img = c1.img ./ max(c1.img(:));
            c2.img = c2.img ./ max(c2.img(:));
            c3.img = c3.img ./ max(c3.img(:));
            
            % binary threshold
            thresh = 0.4;
            brainmask = ((c1.img > thresh) + (c2.img > thresh) + (c3.img > thresh)) > 0;
           
            % saving brainmask
            if job.brainmaskproc.writebrainmask
                im_nii.img = brainmask;
                save_untouch_nii(im_nii, fullfile(pthT1, [nameT1,'_brainmask',extT1]));
            end
        end
    end
    
    
    % loading image if internal regularisation
    if job.biasproc.biasreg > 0
        im_reg = load_untouch_nii(fullfile(pthT1,['m',nameT1,extT1]));
        current_img = im_reg.img;
        current_img = double(current_img) .* double(brainmask);
        
    end
    
    
    % -----------------------------------------------------------------
    % Tissue estimation using a Fuzzy C-Means approach
    % - 3 classes are estimated.
    % - Matlab implementation based on the method of Mahdi Amiri 
    %   http://ce.sharif.edu/~m_amiri/
    % -----------------------------------------------------------------

    % lesion voxels are removed from the input mask before segmentation
    brainmask = and(brainmask, not(lesionmask));
    current_img = double(current_img) .* double(brainmask);


    % image segmentation pre-filtering to remove possible non-brain parts such a
    % eyes, neck, etc.

    m_int = mean(current_img(brainmask));
    s_int = std(current_img(brainmask));
    outlier_voxels = current_img < (m_int + 3*s_int);
    constr_brainmask = and(brainmask, outlier_voxels);

    % tissue estimation
    [centers, members, ~] = tissue_estimation(current_img(constr_brainmask),3,[2.0,100,1e-3,0,0]); 
    % reassign labels to classes
    [~, index] = sort(centers);
    tissue = zeros(size(members));
    tissue(1,:) = members(index(1),:);
    tissue(2,:) = members(index(2),:);
    tissue(3,:) = members(index(3),:);
    % segmentation output
    max_tissue = max(tissue);
    segmentation = (tissue == repmat(max_tissue,3,1));
    wm_tissue = zeros(size(current_img));
    wm_tissue(constr_brainmask) = segmentation(3,:);

    % option to save the FCM segmentation (for debug purposes)

    %csf_tissue = zeros(size(current_img));
    %gm_tissue = zeros(size(current_img));
    %csf_tissue(constr_brainmask) = segmentation(1,:);
    %gm_tissue(constr_brainmask) = segmentation(2,:);
    %out_seg = csf_tissue + gm_tissue.*2 + wm_tissue.*3;
    % tmp_fcm = im_nii;
    % tmp_fcm.img = out_seg;
    % save_untouch_nii(tmp_fcm , ['tmp_fcm_seg_',date,'.nii']);


    % -----------------------------------------------------------------
    % lesion refilling. Each lesion is refilled acording to the mean
    % and half of the standard deviation of WM slice intensity.
    % -----------------------------------------------------------------

    slices = size(current_img, 3);  % <- be careful with that (z direction)
    for i=1:slices
        % current slice 
        slice_seg = wm_tissue(:,:,i);
        slice_lesion = logical(lesionmask(:,:,i));
        slice_image = input_image(:,:,i);
        % refill intensities
        sel_voxels = slice_image(slice_lesion);
        m_slice_int = median(double(slice_image(slice_seg == 1))); 
        s_slice_int = std(double(slice_image(slice_seg == 1)));
        slice_image(slice_lesion) = normrnd(double(m_slice_int), s_slice_int/2, size(sel_voxels));
        % save refilled lesions
        input_image(:,:,i) = slice_image;
    end

    % ---------------------------------------------------------------------
    % Save the results 
    %   - output image has the same heading than input image but lesion
    %   voxels have been refilled
    % ---------------------------------------------------------------------
    
    
    filled_image = im_nii;
    filled_image.img = input_image;
    save_untouch_nii(filled_image, refilled_scan);
    
    
    
end
end












