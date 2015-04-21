%**************************************************************************
% FIX TITLE STRING
% ----------------
% This function fixes a string before it is used as a title for a graph.
% It is known that an underscore '_' is used in the title string to set the
% next character to lowercase. But, in many cases we don't need this
% feature, we just want to print the underscore. The correct way is to
% insert the '\_' pair instead of '_'.
% This function is doing exactly this. It can handle any number of
% underscores in the string.
%
% Written by Alecsander Eitan, Qualcomm Israel ltd. September 2004.
% 
%**************************************************************************
function s1 = fix_title_string(s1)

k=1;
while (k==1)
    j=-1;
    n=length(s1);
    if (s1(1)=='_')
        loc=1;
    else
        loc=0;
    end
    if (loc==0)
        for (i=2:n)
            if ((s1(i)=='_') && (s1(i-1)~='\'))
                loc=i;
            end
        end
    end
    if (loc>0)
        s2 = [s1(1:(loc-1)) '\_' s1((loc+1):end)];
        s1 = s2;
        k=1;
    else
        k=0;
    end
end
        