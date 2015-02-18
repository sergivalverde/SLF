function dist = euclidean_distance_vectors (V, X)
% euclidean_distance_vectors: Compute the euclidean distance between two vectors
% Part of the SLF toolbox
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  Seexs the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright 2014 Sergi Valverde/Arnau Oliver/ Xavier Lladó
%   $Revision: 1.0$     $Date: 12/03/14$ 
%===============================================================================

c = size(V, 1);
p = size(V, 2);
n = size(X, 1);

dist = zeros(c, n);

% fill the output matrix
if p > 1,
    for k = 1:c,
	dist(k, :) = sqrt( sum(((X-ones(n, 1)*V(k, :)).^2)') );
    end
else
    % 1-D data
    for k = 1:c,
	dist(k, :) = abs(V(k)-X)';
    end
end
