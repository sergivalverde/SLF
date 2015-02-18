function [V, U, E] = tissue_estimation (X, c, options, init_V)
% tissue estimation: Fuzzy C-Means Clustering with initial prefiltering
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




if nargin < 2,
	error('Too few input arguments!');
end

if nargin > 4,
	error('Too many input arguments!');
end

% Change the following to set default options
default_options = [2;	% weighting exponent (m)
		100;	% max. number of iteration
		1e-3;	% termination threshold
		1;      % info display during iteration 
        0];     % use provided init_V 

if nargin == 2,
	options = default_options;
else
	% If "options" is not fully specified, pad it with default values.
	if length(options) < 5,
		tmp = default_options;
		tmp(1:length(options)) = options;
		options = tmp;
	end
	% If some entries of "options" are nan's, replace them with defaults.
	nan_index = find(isnan(options)==1);
	options(nan_index) = default_options(nan_index);
end

m = options(1);	        	% Weighting exponent
max_iter = options(2);		% Max. iteration
term_thr = options(3);		% Termination threshold
display = options(4);		% Display info or not
use_init_V = options(5);    % use provided init_V

if m <= 1,
    error('The weighting exponent should be greater than 1!');
end


n = size(X, 1);
p = size(X, 2);

E = zeros(max_iter, 1);	% Array for termination measure values

if use_init_V,
    V = init_V;
else
    V = rand(c, p);   
end

%U = zeros (c, n);

% Main loop
for i = 1:max_iter,
	[V, U, E(i)] = fuzzy_next_step(X, V, c, m);

    if display, 
		fprintf('Iteration count = %d, Termination measure value = %f\n', i, E(i));
	end

    % check termination condition
	if E(i) <= term_thr, break; end,
end

iter_n = i;	% Actual number of iterations 
E(iter_n+1:max_iter) = [];
