classdef TrivialRep < replab.Rep
% Describes copies of the real or complex trivial representation of a group

    methods

        function self = TrivialRep(group, field, dimension)
            assert(isa(group, 'replab.CompactGroup'));
            args = {'isUnitary', true, 'isIrreducible', dimension == 1, ...
                            'trivialDimension', dimension, 'frobeniusSchurIndicator', dimension};
            if strcmp(field, 'R') && dimension == 1
                args = horzcat(args, {'isDivisionAlgebraCanonical' true});
            end
            self@replab.Rep(group, field, dimension, args{:});
        end

    end

    methods (Access = protected)

        function rho = image_exact(self, g)
            rho = replab.cyclotomic.eye(self.dimension);
        end

        function rho = image_double_sparse(self, g)
            rho = speye(self.dimension);
        end

        function e = computeErrorBound(self)
            e = 0;
        end

        function c = computeConditionNumberEstimate(self)
            c = 1;
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = headerStr@replab.Rep(self); % logic in parent class
        end

        % Rep

        function b = isExact(self)
            b = true;
        end

        function complexRep = complexification(self)
            assert(self.overR, 'Representation should be real to start with');
            complexRep = replab.rep.TrivialRep(self.group, 'C', self.dimension);
        end

    end

end
