function [ rcc ] = rc_dir( mat, riParam, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[n, m] = size(mat);
if n ~= m
    error('Matrix must be square.');
end
if ~isvector(riParam)
    error('Variable in question for each neuron must be a vector.');
end
if length(riParam) ~= n
    error('Cannot be more than one variable per vertex/node. Parameter vector is not the same size as dims 1 or 2 of the supplied matrix.')
end
binMat = mat ~= 0; % Transform to adjacency matrix
rcc = zeros(1, n);
j = 1;
[sortedRIP, Ind] = sort(riParam, 'ascend'); % Smallest -> largest of the parameter
binMat = binMat(Ind,Ind);
for i = 1:numel(riParam)
    gr =  sortedRIP >= sortedRIP(i); % Vector of all nodes with variable val geq to the value in qestion
    num = sum(gr);
    if num < 2 % Default to a rich club coefficient of 1 for the richest node in the network.
        rcc(j) = 1;
        continue;
    end  
    rcc(j) = nnz(binMat(gr, gr)) / (num*(num-1));
    j = j + 1;
end

end

