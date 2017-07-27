% [raster, binunit] = ASDFToSparse(asdf)
%
%    asdf - {n_neu + 2, 1} ASDF data.
%    Optional args:
%    'row'/'col' ('column') - in the raster that is returned, if 'row' then
%    each neuron's spike train will be a row, if 'col', then it will be a
%    column. If 'row' the raster will have a number of rows equal to the
%    number of neurons, otherwise it will have a number of columns equal to
%    the number of neurons.
%    binsize - the size of the desired time bin in the raster if not the
%    same as the one specified in the asdf file. If the differences in
%    spike timing are less than the bin size, the returned raster will
%    have how many times a given neuron spiked in each bin.
%
% Returns :
%    raster - (n_neu, duration) time raster expressed as a sparse matrix,
%    optionally in row or column (default) major format
%    binunit - the unit of time in the original asdf argument
%
% Description :
%    This function converts the time raster of ASDF to a raster in the form of a
% 	 dense or sparse matrix depending on its size. If the given (optional)
%    bin size is larger than the bin size of the ASDF file, the resulting
%    raster will contain the number of spikes in each bin. The resulting
%    raster starts where the first spike starts, so the first spike time
%    can be arbitrary.

%==============================================================================
% Copyright (c) 2017, The Trustees of Indiana University
% All rights reserved.
%
% Authors:
% Zach Tosi (ztosi@indiana.edu),
% Michael Hansen (mihansen@indiana.edu),
% Shinya Ito (itos@indiana.edu)
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%   1. Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%   2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%   3. Neither the name of Indiana University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================

function [raster, binunit] = ASDFToRaster(asdf, varargin)
%tic;
narginchk(1, 3);
if ~isempty(asdf{end})
    info = asdf{end};
    n_neu = info(1);
    binunit = asdf{end - 1};
else 
    n_neu = length(asdf)-2;
    binunit = 0.25;
end
% change to column vectors if need be
asdf = cellfun(@(x) reshape(x, length(x), 1), ...
    asdf, 'UniformOutput', false);
colMajor = true;
if isempty(varargin)
    binsize = binunit;
else
    for ii=1:length(varargin)
        if isnumeric(varargin{ii})
            binsize = varargin{ii};
            if binsize < binunit
                error('Cannot have custom bin size less than asdf time bin.');
            end
        else
            if strcmp(varargin{ii}, 'row') || ...
                    strcmp(varargin{ii}, 'row major')
                colMajor = false;
            elseif strcmp(varargin{ii}, 'column') || ...
                    strcmp(varargin{ii}, 'col')  || ...
                    strcmp(varargin{ii}, 'column major')
                colMajor = true;
            else
                error('Unknown input');
            end
        end
    end
end
nemp = find(cellfun(@isempty, asdf(1:end-2)));
nemp = [setdiff(1:n_neu, nemp) n_neu+1 n_neu+2];
asdf = asdf(nemp);
n_neu = length(nemp)-2;

minVal = min(cellfun(@min, asdf(1:end-2)));

% very simple check of validity of ASDF
if n_neu ~= size(asdf,1) - 2
    error('Invalid n_neu information is contained in this ASDF');
end

% Get rid of unnecessary metadata
asdf = asdf(1:(length(asdf)-2));

% Convert the spike times to bin indices based on the time unit for each bin
% and the specified size of each bin in the resulting raster
convToBinInd = @(x) ceil(((x-minVal+binsize) ./ binsize)) ;
asdf = cellfun(convToBinInd, asdf, 'UniformOutput', 0);
if binsize > binunit
    bins = cellfun(@unique, asdf, 'UniformOutput', false);
    numSpks = cellfun(@(x) diff(find(diff([0; x; 0]) ~= 0)), ...
        asdf, 'UniformOutput', false);
    V = cell2mat(numSpks);
    I = cell2mat(bins);
    J = cell2mat(cellfun(@(x,y) ones(length(x),1).*y, bins, ...
        num2cell((1:n_neu)'), 'UniformOutput', 0 ));
else
    I = cell2mat(asdf);
    V = ones(size(I));
    J = cell2mat(cellfun(@(x,y) ones(length(x),1).*y, asdf, ...
        num2cell((1:n_neu)'), 'UniformOutput', 0 ));
end

if colMajor
    raster = sparse(I, J, V);
else
    raster = sparse(J, I, V);
end

%toc;
end
