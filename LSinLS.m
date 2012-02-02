function l = LSinLS( LS1, LS2 )
% LSINLS to indicate whether LS1 is fully contained inside LS2 or not
l = false;
switch size(LS1.T,2)
    
    case 0 % 1D problem - 0D level set
        if (max(LS1.X)<max(LS2.X))&&(min(LS1.X)>min(LS2.X))
            l = true;
        end
        
    case 1 % 2D problem - 1D level set
        x1 = LS1.X(LS1.T(:,1),1);
        y1 = LS1.X(LS1.T(:,1),2);
        x2 = LS2.X(LS2.T(:,1),1);
        y2 = LS2.X(LS2.T(:,1),2);
        if all(inpolygon( x1, y1, x2, y2 ))
            l = false;
        end
        
end
