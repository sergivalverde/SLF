function write_log(log_file, message, status)
% SLF refills white matter lesions with the mean and half of the standard
% deviation of the WM tissue of the input image.
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
%   $Revision: 1.2$     $Date: 06/05/14$
%       - 1.1: added support for ANALYZE and log
%	- 1.2: alternative processing to Image Processing toolbox
%===============================================================================    

    % if status is reset, start a new file
    if strcmp(status, 'reset')
       fid = fopen(log_file,'wt');
       fprintf(fid, '*****************  SLFToolbox Log: %s  *******************\n\n', fix_clock(fix(clock)));
    else
        fid = fopen(log_file,'at');
        fprintf(fid, '%s %s   \n', fix_clock(fix(clock)), message);
    end
end

% function to format the clock
function ftime = fix_clock(clock)

  ftime = [ num2str(clock(3)),'/',num2str(clock(2)),'/', num2str(clock(1)),'-',...
            num2str(clock(4)),':',num2str(clock(5)),':',num2str(clock(6)), ' '];

end
