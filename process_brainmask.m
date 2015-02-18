function [bm] = process_brainmask(input_image, numn)
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
%
%   input:
%       - input image containing binary mask of the brain
%   	- number of neighbors
%   output:
%       - processed brainmask
%   
%===============================================================================

    % size components
    [x,y,s] = size(input_image);
    bm = zeros(x,y,s);

    % for each slice, we remove outlier voxels.
    % quite rudimentary
    for i=1:s
        bm_s = input_image(:,:,i);
        bm_dmp = zeros(size(bm_s));
        for j=numn+1:x-numn-1
            for k=numn+1:y-numn-1
                neig = sum(sum(bm_s(j-numn:j+numn,k-numn:k+numn)));
                % if more that half of the neighbor voxels are filled, the
                % current voxel is also filled 
                if neig > (((numn*2+1)^2)/2)
                    bm_dmp(j,k) = 1;
                end
            end
        end
        bm(:,:,i) = bm_dmp;
    end

end
