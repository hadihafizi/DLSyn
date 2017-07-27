% [ str_avas, szes, lens, rast ] = findAvas(asdf, ts, varargin)
%
% Returns :
%    str_avas - a cell array containing all avalanches
%    szes - the size of each avalanche in number of spikes
%    lens - the duration of each avalanch
%    rast - a spike raster in row major ordering over the time range
%    specified, if that was specified
%
% Required args:
%    asdf - spike trains to find avalanches in, in another spike data
%    format (asdf)
%    ts - the time bin to use for finding avalanches
%
% Optional args:
%     Example: [...] = findAvas(asdf, ts, 'offset', a, 'range', b);
%     Finds avalanches only in a specified region of the raster between
%     times (in ts) a and a+b (both inclusive)
%
%     Example: [...] = findAvas(asdf, ts, 'threshold', theta);
%     Avalanches are defined as periods of activity greater than some
%     threshold, by default it is mean activity in each time bin/2, but
%     this can be set manually and this is how
%
%     Arguments can be entered in any order so long as argument name
%     preceeds argument.
%
% Description :
%    A high-performance, vectorized function which takes an asdf cell array
%    and a time bin and finds neuronal avalanches present in the spike trains
%    stored in the asdf when binned in the given time bin. 
%    It returns the size, duration, and avalanches
%    themselves, defined as all the times each neuron spiked in a given
%    avalanche. This can be performed over the whole asdf or some range of
%    time bins. The threshold activity level of avalanches can also be set.
%
% DEPENDENTS :
%    ASDFToRaster.m

%==============================================================================
% Copyright (c) 2017, The Trustees of Indiana University
% All rights reserved.
%
% Authors:
% Zach Tosi (ztosi@indiana.edu),
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

function [ str_avas, szes, lens, rast ] = findAvas( asdf, ts, varargin)
%tic;
narginchk(2, 8);

% Defaults...
offset = 0;
% range needs size of rast first...
customThresh = false;
defaultRange = true;

% assign values if specified
for i=1:2:length(varargin)
    switch varargin{i}
        case 'offset'
            offset = varargin{i+1};
            defaultRange=false;
        case 'range'
            range = varargin{i+1};
            defaultRange=false;
        case 'threshold'
            thresh = varargin{i+1};
            customThresh = true;
        otherwise
            error('Unknown input.');
    end
end

% operations are optimized when the raster is a sparse matrix such that
% each column represents a time bin since str_avas requires us to take sub
% matrices of the raster from one time index to another. See CSC sparse
% data format for why this is the case...
rast=ASDFToRaster(asdf, ts, 'row');
if defaultRange
    range = size(rast,2);
end
% get only specified range
if ~defaultRange
    rast=rast(:, int32((offset+1):(range+offset)));
end

% find number of spikes in each time bin
pop_fir = full(sum(rast));
if ~customThresh
    thresh = mean(nonzeros(pop_fir))/2;
end
% put 1s where a time bin is within a valid avalache, 0 otherwise
evts = pop_fir>thresh;

% pad and find where avalanches start or end based on where 0s change to
% 1s and vice versa
act_change = diff([0 evts 0]);

% 1s indicate an avalanche began in the given time bin
starts = find(act_change == 1);
% -1s indicate an avalanche ended in the previous time bin
ends = find(act_change== -1) - 1;

% durations...
lens = ends-starts;
str_avas = cell(length(lens), 1);
szes = zeros(size(lens));

for ii=1:length(lens)
    % select sub arrays of valid avalanches
    str_avas{ii} = rast(:,starts(ii):ends(ii));
    % find the number of spikes in each valid avalanche
    szes(ii) = nnz(str_avas{ii});
end

%toc;
end

