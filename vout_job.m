function vout = vout_job(job)
%SLF refills white matter lesions with the mean and half of the standard
%deviation of the WM tissue of the input image.
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


    % ******************************************************************
    % this function creates a new structure with the out lesion-filled  
    % files. This structure will be used by the vout function inside the 
    % tbx_cfg_SLFtoolbox to pass the output files to the next SPM process
    % if necessary
    % ******************************************************************
    
    % find all filled images 
    filled = struct('t1',{});

    t1_elements = size(job.data_T1, 1);
    parts = cell(t1_elements, 4);

    % load image paths
    for j=1:t1_elements
        [parts{j,:}] = spm_fileparts(job.data_T1{j});
    end

    % save the correspondent lesion-filled image location
     for i=1:t1_elements
        filled(i).t1{1} = fullfile(parts{i,1},[parts{i,2},'_slf_filled.nii']);
     end
    
     vout  = struct('filled',filled);

end

