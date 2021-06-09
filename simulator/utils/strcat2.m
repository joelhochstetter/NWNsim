function outstring = strcat2(cellstrings, fmt)
%{
    On the cell array cellstrings combines cellstr into a single character
    array. For numeric types applies num2str with given format (fmt).
    No format given behaves as default num2str
    Returns a character array outrstring


    Written by Joel Hochstetter
%}
    if nargin < 2
        fmt = -1;
    end
    
    outstring = '';
    for i = 1:numel(cellstrings)
        curr = cellstrings{i};
        if isnumeric(curr)
            if fmt == -1
                curr = num2str(curr);
            else
                curr = num2str(curr, fmt);
            end
            outstring = strcat(outstring, curr);
        else
            outstring = strcat(outstring, curr);
        end
    end
    
end