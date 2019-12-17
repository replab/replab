function [lawName argFuns] = collectTypes(laws, methodName)
% Splits the given "methodName" according to the law naming convention
%
% (see `replab.Laws` documentation)
%
% Returns
% - "lawName", the function name transformed from camel case
% - "isRandom", whether the law check method needs any random elements
% - "argFuns", a cell array of function handles that return random elements
%   whose types correspond to the method name "signature"
    parts = strsplit(methodName, '_');
    argFuns = {};
    lastPart = parts{end};
    if isstrprop(lastPart(1), 'upper')
        % has random arguments
        nameParts = parts(2:end-1);
        desc = [parts{end} '#']; % add a canary at the end
        i = 1;
        while true
            switch desc(i)
              case '#'
                break
              case 'Z'
                i = i + 1;
                digits = '';
                while isstrprop(desc(i), 'digit')
                    digits = [digits desc(i)];
                    i = i + 1;
                end
                assert(length(digits) > 0, 'Error in natural number description');
                f = @() randi(str2num(digits));
              case 'N'
                i = i + 1;
                digits = '';
                while isstrprop(desc(i), 'digit')
                    digits = [digits desc(i)];
                    i = i + 1;
                end
                assert(length(digits) > 0, 'Error in natural number description');
                if digits(1) == '0'
                    assert(length(digits) > 1, 'Error in natural number description');
                    f = @() randi([0 str2num(digits(2:end))]);
                else
                    f = @() randi(str2num(digits));
                end
              otherwise
                try
                    domain = laws.(desc(i));
                catch
                    error(sprintf('Type %s is unknown', desc(i)));
                end
                f = @() domain.sample;
                i = i + 1;
            end
            argFuns{end+1} = f;
        end
    else
        % has not any random arguments
        nameParts = parts(2:end);
    end
    lawName = strjoin(nameParts, ' ');
end
