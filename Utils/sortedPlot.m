function sortedPlot( x, y, varargin)
%% A utility function for plotting x against y where x is sorted in 
% ascending order first and y is arranged arranged according to the same 
% ordering. Any arguments that can be passed to the normal plot function
% can be passed to this function.
    [xp, I] = sort(x);
    plot(xp, y(I), varargin);

end

