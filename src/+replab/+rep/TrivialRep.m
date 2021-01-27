classdef TrivialRep < replab.Rep
% Describes copies of the real or complex trivial representation of a group

    properties (SetAccess = protected)
        isExactValue % (logical): What to return for `.isExact`
    end

    methods

        function self = TrivialRep(group, field, dimension, isExactValue)
            if nargin < 4 || isempty(isExactValue);
                isExactValue = true;
            end
            assert(isa(group, 'replab.CompactGroup'));
            args = {'isUnitary', true, 'isIrreducible', dimension == 1, ...
                            'trivialDimension', dimension, 'frobeniusSchurIndicator', dimension};
            if strcmp(field, 'R') && dimension == 1
                args = horzcat(args, {'isDivisionAlgebraCanonical' true});
            end
            self@replab.Rep(group, field, dimension, args{:});
            self.isExactValue = isExactValue;
        end

    end

    methods (Access = protected)

        function b = computeIsUnitary(self)
            b = true;
        end

        function r = computeDouble(self)
            r = replab.rep.TrivialRep(self.group, self.field, self.dimension, false);
        end

        function rho = image_exact(self, g)
            assert(self.isExact);
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
            b = self.isExactValue;
        end

        function p = invariantBlocks(self)
            blocks = arrayfun(@(i) {i}, 1:self.dimension, 'uniform', 0);
            % each coordinate in its own block
            p = replab.Partition.fromBlocks(blocks);
        end

        function complexRep = complexification(self)
            assert(self.overR, 'Representation should be real to start with');
            complexRep = replab.rep.TrivialRep(self.group, 'C', self.dimension);
        end

    end

end
