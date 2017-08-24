function [layerMem] = makeLayersFromMat(adj_mat, varargin)

%narginchk(1, 3);
%nargoutchk(1,2);

[m, n] = size(adj_mat);
if m~=n
    % For now, just to simplify life...
    error('Matrix must be square');
end
%adj_mat = double(adj_mat~=0);
% Replace the above with the below for a weighted version
%adj_mat = abs(adj_mat);

indegs = sum(adj_mat)';
outdegs = sum(adj_mat,2);
recNess = (indegs - outdegs) ./ (indegs + outdegs);
adj_matT = adj_mat';
[~, recNI] = sort(recNess, 'ascend'); % TODO make sure to include all with 0 in-deg

if nargin == 1
    initL1 = recNI(1:int32(size(adj_mat,1)/5));
    %initL1 = randperm(size(adj_mat,1), int32(size(adj_mat,1)/5));
end

lscores = zeros(n, 2); % Layer, layer score
lscores(initL1,1) = 1;
lscores(initL1,2) = indegs(initL1)./2;
visited = zeros(n,1);
visited(initL1) = 1;
currLay = 1;
maxit = 100;


for ii=1:100
    if ii > 1
        visited = zeros(n,1);
        visited(lscores(:,1)==currLay) = 1;
        visited(indegs==0) = 1;
    end
    %lscores(:,2) = 0;
    stLay = lscores(:,1);
    itcnt = 0;
    while any(~visited) && itcnt < maxit
        layIni = double(lscores(:,1) == currLay);
        %size(layIni)
        LNextScore = adj_matT * layIni;
        visited = visited | (double(adj_matT~=0)*layIni);
        mvL = (lscores(:,2) < LNextScore) | (rand(size(LNextScore)) < ((100-ii)./10000)) ;
        nextLay = lscores(:,1) == (currLay+1) | mvL;
        penalty = sum(adj_mat(mvL, nextLay),2);
        LNextScore(mvL) = LNextScore(mvL)-penalty;
        mvL = (lscores(:,2) < (LNextScore));
        %visited = visited & ~mvL;
        lscores(mvL,2) = LNextScore(mvL);
        lscores(mvL,1) = currLay+1;
        
        %lscores(mvL,2) = lscores(mvL,2)-sum(adj_mat(mvL, nextLay),2);
        currLay = currLay+1;
        %any(~visited)
        %sum(mvL) ~= 0
        itcnt = itcnt+1;
    end
    disp(sum(stLay ~= lscores(:,1)));
    lscores(:,1) = lscores(:,1) - min(lscores(:,1))+1;
    currLay = max(lscores(:,1));
    visited = zeros(n,1);
    visited(lscores(:,1)==currLay) = 1;
    visited(indegs==0) = 1;
    %maxLay =  max(lscores(:,1));
    %lscores(:,2)=0;
    itcnt=0;
    while any(~visited) && itcnt < maxit
        layIni = double(lscores(:,1) == currLay);
        %size(layIni)
        LNextScore = adj_mat * layIni;
        visited = visited | (double(adj_mat~=0)*layIni);
        mvL1 = (lscores(:,2) < LNextScore) | (rand(size(LNextScore)) < ((100-ii)./10000)) ;
        %visited = visited & ~mvL;
        nextLay = lscores(:,1) == (currLay-1) | mvL1;
        penalty = sum(adj_mat(nextLay, mvL1))';
        LNextScore(mvL1) = LNextScore(mvL1)-penalty;
        mvL = (lscores(:,2) < LNextScore);
        %         if ~any(mvL) && any(mvL1) %no movers after penalty
        %             lscores(lscores(:,1)<currLay) = lscores(lscores(:,1)<currLay);
        %         end
        
        lscores(mvL,2) = LNextScore(mvL);
        lscores(mvL,1) = currLay-1;
        currLay = currLay-1;
        
        %any(~visited)
        %sum(mvL) ~= 0
        itcnt = itcnt+1;
    end
    %lscores(:,2)=0;
    lscores(:,1) = lscores(:,1) - min(lscores(:,1))+1;
end



layerMem = lscores;