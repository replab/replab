function flags = parseFlags(s, errFun)
% Parses the flags given in a ``doctest: XXX`` comment
%
% Args:
%   s (charstring): Flag string to parse
%   errFun (function_handle): Function called before an error is thrown;
%                             can be used to display the error context
%
% Returns:
%   struct: Flags structure
    errId = 'replab:docTestParseError';
    if nargin < 2
        errFun = @() [];
    end
    s = strtrim(s);
    flags = struct;
    while true
        token = regexp(s, '^((\+|-)?[A-Za-z][A-Za-z_0-9]*)', 'tokens', 'once');
        if isempty(token)
            errFun();
            error(errId, 'Invalid syntax in flags');
        end
        if iscell(token)
            token = token{1};
        end
        s = s(length(token)+1:end);
        if token(1) == '+'
            k = token(2:end);
            v = true;
        elseif token(1) == '-'
            k = token(2:end);
            v = false;
        else
            if s(1) ~= '('
                errFun();
                error(errId, sprintf('Unrecognized argument list for %s token', token));
            end
            s = s(2:end);
            ind = find(s == ')', 1);
            if isempty(ind)
                errFun();
                error(errId, sprintf('Argument list of %s should be closed by )', token));
            end
            args = s(1:ind-1);
            s = s(ind+1:end);
            k = token;
            v = strsplit(args, ',');
        end
        if isfield(flags, k)
            errFun();
            error(errId, sprintf('Flag %s defined twice', k));
        end
        flags.(k) = v;
        s = strtrim(s);
        if isempty(s)
            break
        end
        if s(1) ~= ','
            errFun();
            error(errId, sprintf('Missing flag separator , after %s', token));
        end
        s = s(2:end);
        s = strtrim(s);
    end
end
