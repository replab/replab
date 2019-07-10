function msg = message(errorDesc, context, assertArgs, testVarNames, testVarValues)
% Formats a failed law check error message
    if ~isequal(context, '')
        msg = [errorDesc char(10) context];
    else
        msg = errorDesc;
    end
    for i = 1:length(assertArgs)
        assertArgs{i} = replab.shortStr(assertArgs{i});
    end
    msg = sprintf(msg, assertArgs{:});
    msg = [msg char(10) 'with variables' char(10)];
    maxLength = max(cellfun(@(x) length(x), testVarNames));
    for i = 1:length(testVarNames)
        vn = testVarNames{i};
        vn = [blanks(maxLength - length(vn)) vn ':' char(10)];
        vv = replab.shortStr(testVarValues{i});
        msg = [msg sprintf('%s\n%s\n\n', vn, vv)];
    end
end
