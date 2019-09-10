classdef Normalized < replab.Str
% Represents a vector vec*factor, where 'vec' is nice to print
    properties (SetAccess = protected)
        vec;
        factor;
    end
    methods
        function self = Normalized(vec, factor)
            assert(factor > 0, 'Factor must be positive');
            self.vec = vec;
            self.factor = factor;
        end
        function s = shortStr(self, maxColumns)
            f = self.factor;
            invsquare = 1/f^2;
            if abs(invsquare - round(invsquare)) < replab.Settings.doubleEigTol
                n = round(invsquare);
                if n == 1
                    rhs = '';
                elseif sqrt(n) == round(sqrt(n))
                    n1 = sqrt(n);
                    rhs = sprintf('/%d', n1);
                else
                    rhs = sprintf('/sqrt(%d)', n);
                end
            else
                rhs = sprintf('* %e', 1/f);
            end
            maxL = maxColumns - length(rhs);
            lhs = replab.shortStr(self.vec, maxL);
            s = [lhs rhs];
        end
    end
end
