function [lawName isRandom argFuns] = collectTypes(laws, methodName)
    parts = strsplit(methodName, '_');
    argFuns = {};
    lastPart = parts{end};
    if isstrprop(lastPart(1), 'upper')
        isRandom = true;
        % has random arguments
        nameParts = parts(2:end-1);
        desc = [parts{end} '#']; % add a canary at the end
        i = 1;
        while true
            switch desc(i)
              case '#'
                break
              case 'D'
                err = sprintf('The object %s needs to be a domain', str(laws));
                assert(isa(laws, 'replab.cat.Domain'), err);
                f = @() laws.sample;
                i = i + 1;
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
        isRandom = false;
        % has random arguments
        nameParts = parts(2:end);
    end
    lawName = strjoin(nameParts, ' ');
end
