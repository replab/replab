function X = convert(X, outputType)
% Converts an internal image to the requested type, if possible
%
% Raises:
%   An error if the conversion is impossible.
%
% Args:
%   rho (double(\*,\*) or intval(\*,\*) or cyclotomic(\*,\*)): Matrix to convert
%   outputType ('double', 'double/sparse', 'intval', 'exact'): Requested type
%
% Returns:
%   double(\*,\*) or cyclotomic(\*,\*) or intval(\*,\*): Converted matrix
    inputType = class(X);
    if strcmp(inputType, 'replab.cyclotomic')
        inputType = 'exact';
    end
    switch [inputType ' -> ' outputType]
        % from exact
      case 'exact -> exact'
        X = X;
      case 'exact -> intval'
        X = X.intval;
      case {'exact -> double', 'exact -> double/sparse'}
        X = X.double;
        % from intval
      case 'intval -> exact'
        assert(all(all(rad(X) == 0)), 'Input matrix needs to be exact');
        X = replab.cyclotomic.fromDoubles(mid(X));
      case 'intval -> intval'
        X = X;
      case 'intval -> double/sparse'
        X = mid(X);
      case 'intval -> double'
        X = full(mid(X));
        % from double (double/sparse is never used as input)
      case 'double -> exact'
        assert(all(all(X == round(X))), 'Input matrix needs to have integer coefficients');
        X = replab.cyclotomic.fromDoubles(X);
      case 'double -> intval'
        assert(all(all(X == round(X))), 'Input matrix needs to have integer coefficients');
        X = intval(X);
      case 'double -> double/sparse';
        X = X;
      case 'double -> double'
        X = full(X);
    end
end
