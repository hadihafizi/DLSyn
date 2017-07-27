function [p] = unityline(ax, varargin)
    hold on;
    xl = ax.XLim;
    yl = ax.YLim;
    
    minV = min([xl;yl]);
    maxV = max([xl;yl]);
    
    if isempty(varargin)
       p = plot([minV(1) maxV(2)], [minV(1) maxV(2)], 'k--'); 
    else
       p = plot([minV(1) maxV(2)], [minV(1) maxV(2)], varargin); 
    end

    hold off;

end