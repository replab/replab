function X = convert(X, type)
% Converts an internal image to the requested type, if possible
%
% Only the following conversions are possible:
% * ``type`` to ``type``: no change
% * ``cyclotomic`` to ``intval`` or ``double`: by approximation
% * ``intval`` to ``double``: by computing the interval midpoint of each coefficient
% * ``double`` to ``intval`` or ``cyclotomic``: only when the double matrix has integer entries
%
% Raises:
%   An error if the conversion is impossible.
%
% Args:
%   rho (double(\*,\*) or intval(\*,\*) or cyclotomic(\*,\*)): Matrix to convert
%   type ({'double', 'intval', 'cyclotomic'}): Requested type
%
% Returns:
%   double(\*,\*) or intval(\*,\*) or cyclotomic(\*,\*): Converted matrix
    switch class(X)
      case 'replab.cyclotomic'
        switch type
          case 'intval'
            [approx error] = X.doubleApproximation;
            X = midrad(approx, error);
          case 'double'
            X = X.double;
          case 'cyclotomic'
          otherwise
            error('Unknown type %s', type);
        end
      case 'intval'
        switch type
          case 'cyclotomic'
            r = round(mid(X));
            if all(all(rad(X) == 0)) && isreal(X) && all(all(mid(X) == r))
                X = replab.cyclotomic.fromDoubles(r);
            else
                error('Cannot convert a non-integer interval matrix to an exact matrix');
            end
          case 'double'
            X = mid(X);
          case 'intval'
          otherwise
            error('Unknown type %s', type);
        end
      case 'double'
        switch type
          case 'cyclotomic'
            r = round(X);
            if isreal(X) && all(all(X == r))
                X = replab.cyclotomic.fromDoubles(r);
            else
                error('Cannot convert a floating-point matrix to an exact matrix');
            end
          case 'intval'
            r = round(X);
            if isreal(X) && all(all(X == r))
                X = intval(r);
            else
                error('Cannot convert a floating-point matrix to an interval matrix');
            end
          case 'double'
          otherwise
            error('Unknown type %s', type);
        end
      otherwise
        error('Unknown parameter type %s', class(X));
    end
end
