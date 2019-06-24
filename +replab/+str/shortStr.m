function [s overLimit] = shortStr(obj, limit)
% Returns a one-line description of the given object, that fits within the given width limit
%
%   obj: Object to pretty print
% limit: maximum column size before switching to crude string description
%        optional parameter with default value 80
%
% Returns the string 's' and a Boolean 'overLimit' that states whether the output has been shortened
% because it ran over limit.
%   
% Thank you Matlab for being a regular language, we totally don't have to consider a thousand
% particular cases below.
    if nargin < 2
        limit = 80;
    end
    
    try
        s = obj.shortStr;
    catch
        if isscalar(obj) % prints scalars using the relevant method
            if  isnumeric(obj)
                s = num2str(obj);
            elseif isa(obj, 'sym')
                s = char(obj);
            else
                s = class(obj);
            end
        elseif iscell(obj) % prints cells
            if isequal(size(obj), [0 0])
                s = '{}'
            elseif numel(obj) == 0
                s = [replab.str.printSize(size(obj)) ' empty cell array'];
            if isrow(obj) && numel(obj) < limit
                s = ['{' strjoin(cellfun(@(x) replab.str.shortStr(x), obj, 'uniform', 0), ', ') '}'];
            else iscolumn(obj)
                [s ol] = replab.str.shortStr(obj');
                if ol > 
                s = [replab.str.shortStr(obj') ''''];
                elseif ismatrix(obj) && numel(obj) < limit
            elseif iscolumn(obj)
                
                
        elseif isrow(obj) && numel(obj) < limit
            s = ['[' strjoin(arrayfun(@(x) replab.str.shortStr(x), obj, 'uniform', 0), ', ') ']'];
        elseif iscolumn(obj)
            s = [replab.str.shortStr(obj') ''''];
        elseif ismatrix(obj) && numel(obj) < limit
            rowFun = @(r) strjoin(arrayfun(@(c) replab.str.shortStr(obj(r, c)), 1:size(obj, 2), 'uniform', 0), ', ');
            s = ['[' strjoin(arrayfun(rowFun, 1:size(obj, 1), 'uniform', 0), '; ') ']'];
        else
            s = class(obj);
        end
        if length(s) > limit
            overLimit = true;
            if ismatrix(obj)
                s = [replab.str.printSize(size(obj))  ' ' class(obj)];
            else
                s = class(obj);
            end
        else
            overLimit = false;
        end
    end
end
