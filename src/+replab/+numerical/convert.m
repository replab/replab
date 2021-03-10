function X = convert(X, type)
% Converts an internal image to the requested type, if possible
%
% Raises:
%   An error if the conversion is impossible.
%
% Args:
%   rho (double(\*,\*) or intval(\*,\*) or cyclotomic(\*,\*)): Matrix to convert
%   type ({'double', 'double/sparse', 'exact'}): Requested type
%
% Returns:
%   double(\*,\*) or cyclotomic(\*,\*): Converted matrix
    if strcmp(type, 'exact')
        if isa(X, 'double')
            X = replab.cyclotomic(X);
        end
        return
    end
    if ~isa(X, 'double')
        X = double(X); % can be sparse
    end
    if strcmp(type, 'double')
        X = full(X);
    end
end
