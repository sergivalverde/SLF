function [V, U, E] = fuzzy_next_step(X, V, c, m)
% fuzzy_ext_step: next step on the fuzzy-C-means algorithm
% Part of the SLF toolbox
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



    n = size (X, 1);
    p = size (X, 2);

    dist = euclidean_distance_vectors(V, X); % fill the distance matrix

    % calculate new U, suppose m != 1
    tmp = dist.^(-2/(m-1));
    U = tmp./(ones(c, 1)*sum(tmp));

    % Correct the situation of "singularity" (one of the data points is
    % exactly the same as one of the cluster centers).

    si = find (tmp == Inf);
    U(si) = 1;
    if (size (si, 1) ~= 0)
        display ('FCMC, Warning: Singularity occured and corrected.');
    end

    % Check constraint
    tmp = find ((sum (U) - ones (1, n)) > 0.0001);
    if (size(tmp,2) ~= 0)
        display ('FCMC, Warning: Constraint for U is not hold.');
    end

    V_old = V;
    mf = (U).^m; % MF matrix after exponential modification
    V = mf*X./((ones(p, 1)*sum(mf'))'); % new center

    E = norm (V - V_old, 1);

