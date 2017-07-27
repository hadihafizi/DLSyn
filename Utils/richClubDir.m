%% Calculates the rich-club coefficient, normalized rich-club coefficient, and the
% significance of the normalized rich-club coefficient for a directed 
% graph, based on an arbitrary set of properties assigned to each vertex.
% In other words: produces values answering the question: are vertices with
% "x" quantified property more (or less) connected to each other than would
% be expected by chance? Where "x" is completely arbitrary.
%
% CITATION:
%   V. Colizza, A. Flamini, M. Serano, and A. Vespignani,
%       "Detecting rich-club ordering in complex networks" Nature physics,
%       vol. 2, no. 2, pp. 110-115,2006
%
% DEPENDENTS:
%
%   dir_generate_srand: 
%       For generating directed degree preserving rewirings
%       AUTHOR: Sergei Maslov: Brookhaven National Laboratory
%
%   semistd: 
%       finds the lower and upper standard deviations for some set of data,
%       here for results the analysis on the null models
%       AUTHOR: Abdulaziz Alhouti: George Washington University
%
%   shadedErrorBar:
%       Uses shaded area to express standard error instead of errorbars for
%       null model means & std devs.
%       AUTHOR:  Rob Campbell: Basel University
%
%   rc_dir: 
%       Performs the actual calculation for the original matrix
%       and each null model
%       AUTHOR: Zach Tosi: Indiana University
%
%
% COMPATIBILITY: Not tested below R2015b, but likely compatible with earlier versions
%
%   nrm_rccs = richClubDir(mat, riParam) calculates the normalized
%   rich-club coefficients in the order of the sorted richness parameter
%   (riParam) for each element in riParam. nrm_rccs will be the same length
%   as riParam, and indicate the normalized rich club coefficient for all
%   members of sorted riParam (ascending) with the value of the riParam at
%   that sorted index or greater.
%
%   [nrm_rrcs, rccs] = richClubDir(mat, riParam) also provides the raw rich
%   club coefficients.
%
%   [nrm_rccs, rccs, h, pValues] = richClubDir(mat, riParam) provides h
%   which indicates the result of the hypothesis test as well as the
%   p-values themselves. Sigificance level is 1/# of null models.
%
%   [...] = richClubDir(mat, riParam, null) allows specification of the
%   number of null models to generate for determining the normalized rccs
%   and the significance values if "null" is a scalar. "null" can also be a
%   3D matrix with the same number of rows/columns as "mat" and with a
%   number of elements along the 3rd dim being the number of null models
%   provided. That is, the user can supply their own null-models either for
%   speed (say if many riParams are being tested for the same graph) or if
%   they want to use a different definition of a null model from a degree
%   preserving rewire (the default). Must be the 3rd argument.
%
%   [...] = richClubDir(mat, riParam, null, 'PARAM1', val1, ...) specifies
%   one or more of the following name/value pairs:
%   
%   Parameter       Value
%   'p-value'       How hypothesis testing is conducted: "empirical" which
%                   calculates p as the largest possible p value as a
%                   proportion of the number of rccs in the null model
%                   greater (if rcc is greater than the null mean) or less
%                   (if rcc is less than the null mean rcc) than the rcc at
%                   that point for the given graph. Ex. if 100 null models
%                   are specified/provided an entry in pValues of 0.01
%                   indicates p < 0.01 meaning that the node represented by
%                   that entry had a rcc greater/less than all rccs for
%                   that node calculated in the null models. The other
%                   option "analytical" assigns rccs z-scores based on the
%                   distribution of rccs for each element produced by the
%                   null models. This method uses different upper and lower
%                   standard deviations, but on either side of the mean
%                   still assumes normalcy of the distribution of null
%                   models. "empricial" is the default.
%
%   'ShowFigure'    Whether or not to show a figure of the results which
%                   displays the normalized rccs, raw rccs, the rccs of the
%                   null models with standard shaded error bars and
%                   highlights regions where h == 1 indicating rejection of
%                   the null hypothesis at the default level. true/false;
%                   default is false--plots to the current axis.
%   
%
%%
% AUTHORS: Zach Tosi
% EMAIL: ztosi@indiana.edu
% March 2016
% Febuary 2017

%%==============================================================================
% Copyright (c) 2016, The Trustees of Indiana University
% All rights reserved.
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
%
function [ nrm_rccs, rccs, h, pValues ] = richClubDir( adjMat, riParam, varargin )
    narginchk(2,7);
    nargoutchk(0, 4);
    [m, n] = size(adjMat);
    type = 'empirical';
    %% Error Checking
    if m~=n
        error('Input matrix is not square');
    end
    if length(riParam) ~= m
        error('Vector of values for each node being tested for richness must have length equal to # of rows (and columns) of the adjacency matrix.')
    end
    
    %Defaults:
    showFigs = false;
    NUM_NULL = 100;
    nullMods = [];
    
    % Sort out inputs...
    adjMat = uint8(adjMat~=0);
    if ~isempty(varargin)
        if isscalar(varargin{1})
            NUM_NULL = varargin{1};
        else
            nullMods = uint8(varargin{1}~=0);
            dims = size(nullMods);
            NUM_NULL = dims(3);
            if dims(1:2) ~= [m ,n]
                error('Null model have the same dimensions as the given matrix.');
            end    
        end
        
        for i=2:2:length(varargin)
            switch varargin{i}
                case 'p-value'
                    type = varargin{i+1};
                    if ~(strcmpi(type, 'empirical') || strcmpi(type, 'analytical'))
                        error('Invalid p-value option. Options are "empirical" or "analytical"');
                    end
                case 'ShowFigure'
                    showFigs = varargin{i+1};
                    if ~isa(showFigs, 'logical')
                        error('ShowFigure option must be boolean.');
                    end
                otherwise
                    error('Invalid input.');        
            end
        end
    end
    
    %% Calculate Rccs for each unique riParam, for the original matrix
    [n, ~] = size(adjMat);
    maxVal = max(riParam);
    [nrm_rccs] = rc_dir(adjMat, riParam); % Perform Rich-Club exam on original
    rccs = nrm_rccs;

    %% Now for the null models
    % If no null models were supplied generate them...
    if isempty(nullMods)
        disp('Generating Null Models...');
        nullMods = zeros(m,n, NUM_NULL, 'uint8');
        parfor i = 1:NUM_NULL
            nullMods(:,:,i) = uint8(dir_generate_srand(adjMat));
        end
    end
    % Calculate the Rccs for the null models
    rccShuff = zeros(NUM_NULL,n);
    parfor i = 1:NUM_NULL
        [rccShuff(i, :)] = rc_dir(nullMods(:,:,i), riParam);
    end
    
    %% Find out how significant the results are and determine the normalized Rccs from the null models
    % This is the order the Rccs are return in, so we need to sort this to
    % make graphs
    sortedRIP = sort(riParam, 'ascend');       
    meanRccSh = mean(rccShuff);
    gr = rccs >= meanRccSh;
    les = rccs < meanRccSh;
    [stdRccShL, stdRccShU] = semistd(rccShuff);  
    if strcmpi(type, 'empirical')
        pValues = zeros(size(rccs));
        pValues(gr) = sum(bsxfun(@ge, rccShuff(:, gr), rccs(gr)));
        pValues(les) = sum(bsxfun(@le, rccShuff(:,les), rccs(les)));
        h = pValues < 1/NUM_NULL;
        pValues = (pValues ./ NUM_NULL) + 1/NUM_NULL;
    else
        pValues(gr) = abs(rccs(gr) - meanRccSh(gr)) ./ stdRccShU(gr);
        pValues(les) = abs(rccs(les) - meanRccSh(les)) ./ stdRccShL(les);
        pValues = normcdf(-pValues, 0, 1);
        pValues(isnan(pValues)) = 1; % Means there was 0 standard dev
        h = pValues < 1/NUM_NULL;
    end

    nrm_rccs = nrm_rccs ./ meanRccSh;
    
    %% Optionally plot the results, showing where the results are significant
    if showFigs
        figure;
        mx = 1:uint32(maxVal/20):uint32(maxVal+maxVal/10);
        maxn = max(nrm_rccs);
        sig = find(h);
        if ~isempty(sig)
            
            if sig(end) ~= length(pValues)
                sig = [sig length(pValues)];
            end
            
            inds = diff(sig);
            %inds
            if sum(inds>1) == 0
                inds = [1 length(sig)];
            else
                if inds(1) == 1
                    inds = [0 find(inds>1)];
                else
                    inds = find(inds>1);
                end
            end
            %inds
            pairs = zeros(length(inds)-1, 2);
            for i = 1:(length(inds)-1)
                pairs(i,:) = [sortedRIP(sig(inds(i)+1)) ...
                    sortedRIP(sig(inds(i+1)))];
            end
            %pairs
            
            X = zeros(4, length(inds)-1);
            Y = zeros(4, length(inds)-1);
            Y(2,:) = maxn + maxn/10;
            Y(3,:) = maxn + maxn/10;
            for i = 1:length(pairs(:,1))
                X(1, i) = pairs(i, 1);
                X(3, i) = pairs(i, 2);
                X(2, i) = pairs(i, 1);
                X(4, i) = pairs(i, 2);
            end
            patchSaturation=0.15;
            faceAlpha=patchSaturation;
            patchColor=[1, .84, 0];
            %set(gcf,'renderer','openGL');
            
            patch(X, Y, patchColor,...
                'edgecolor','none',...
                'facealpha',faceAlpha);
            hold on;
        end
       
        rc = plot(sortedRIP, rccs, 'r', 'LineWidth', 2);
        hold on;
        seb = shadedErrorBar(sortedRIP, ...
                meanRccSh, ...
                 [stdRccShU; stdRccShL], {'Color', [0, .45, .74], ...
                 'LineWidth', 2}, 1);
        nrc = plot(sortedRIP, nrm_rccs, 'k', 'LineWidth', 2);
        plot(mx, ones(1, numel(mx)), '--k', 'LineWidth', 2);
        xlim([0 maxVal+maxVal/20]);
        ylim([0 maxn+maxn/20]);
        dummy = patch([0 0 0 0], [0 0 0 0], [1 .97 .75]);
        legend([nrc rc seb.mainLine dummy], '\phi_{norm}',...
            '\phi', '\phi_{null}', ['p < ' num2str(1/NUM_NULL)], ...
            'Location', 'northwest');
        hold off;
    end
    hold off;
    
end

