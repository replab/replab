classdef EncodedRep < replab.Rep
% Encoding of a real/complex representation into a complex/real representation
%
% The following types are available:
%
% * 'R^d -> C^d': no transformation of coefficients, the images are real
% * 'C^d -> R^d': no transformation of coefficients, the images are real
% * 'R^2d -> C^d': parent is a real complex-type representation with divisionAlgebraName set
% * 'C^d -> R^2d': this will be a real complex-type representation with divisionAlgebraName set
% * 'R^4d -> C^2d': parent is a real quaternion-type representation
% * 'C^2d -> R^4d': parent is a complex quaternion-type representation

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Real representation being transformed
        type % ('R^d -> C^d', 'R^2d -> C^d', 'C^d -> R^2d', 'R^4d -> C^2d', 'C^2d -> R^4d'): Type of encoding
    end

    methods

        function self = EncodedRep(parent, type)
            switch type
              case 'R^d -> C^d'
                assert(parent.overR);
                field = 'C';
                dimension = parent.dimension;
                divisionAlgebraName = '';
              case 'C^d -> R^d'
                assert(parent.overC);
                field = 'R';
                dimension = parent.dimension;
                divisionAlgebraName = '';
              case 'R^2d -> C^d'
                assert(parent.overR);
                assert(strcmp(parent.divisionAlgebraName, 'complex'));
                field = 'C';
                dimension = parent.dimension/2;
                divisionAlgebraName = '';
                error('TODO');
              case 'C^d -> R^2d'
                assert(parent.overC);
                field = 'R';
                dimension = parent.dimension*2;
                divisionAlgebraName = 'complex';
                error('TODO');
              case 'R^4d -> C^2d'
                assert(parent.overR);
                assert(strcmp(parent.divisionAlgebraName, 'quaternion.rep'));
                field = 'C';
                dimension = parent.dimension/2;
                divisionAlgebraName = '';
                error('TODO');
              case 'C^2d -> R^4d'
                assert(parent.overC);
                field = 'R';
                dimension = parent.dimension*2;
                divisionAlgebraName = 'quaternion.rep';
                error('TODO');
            end
            self@replab.Rep(parent.group, field, dimension, 'isUnitary', parent.isUnitary, 'divisionAlgebraName', divisionAlgebraName);
            self.parent = parent;
            self.type = type;
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_complexifiedTrivial(self, options)
        % Encoded trivial is still trivial
            if isa(self.parent, 'replab.rep.TrivialRep')
                res = replab.rep.TrivialRep(self.parent.group, self.field, self.dimension);
            else
                res = [];
            end
        end

    end

    methods (Access = protected)

        % Rep

        function c = decomposeTerm(self)
            c = {self.parent};
        end

        function r = composeTerm(self, newParents)
            r = replab.rep.EncodingRep(newParents{1}. type);
        end

        function rho = image_double_sparse(self, g)
            switch self.type
              case 'R^d -> C^d'
                rho = self.parent.image(g, 'double/sparse');
              case 'C^d -> R^d'
                % kill imaginary error part if present
                rho = real(self.parent.image(g, 'double/sparse'));
              otherwise
                error('TODO');
            end
        end

        function rho = image_exact(self, g)
            switch self.type
              case 'R^d -> C^d'
                rho = self.parent.image(g, 'exact');
              case 'C^d -> R^d'
                % should be exactly real
                rho = self.parent.image(g, 'exact');
              otherwise
                error('TODO');
            end
        end

        function e = computeErrorBound(self)
            switch self.type
              case 'R^d -> C^d'
                e = self.parent.errorBound;
              otherwise
                e = inf;
            end
        end

        function c = computeConditionNumberEstimate(self)
            switch self.type
              case {'R^d -> C^d', 'C^d -> R^d'}
                c = self.parent.conditionNumberEstimate;
              otherwise
                error('TODO');
            end
        end

        function k = computeKernel(self)
            k = self.parent.kernel;
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function M = matrixRowAction_double_sparse(self, g, M)
            switch self.type
              case 'R^d -> C^d'
                  M = self.parent.matrixRowAction(g, real(M)) + ...
                      self.parent.matrixRowAction(g, imag(M)) * 1i;
              case 'C^d -> R^d'
                M = real(self.parent.matrixRowAction(g, M));
              otherwise
                M = matrixRowAction_double_sparse@replab.Rep(self, g, M);
            end
        end

        function M = matrixColAction_double_sparse(self, g, M)
            switch self.type
              case 'R^d -> C^d'
                M = self.parent.matrixColAction(g, real(M)) + ...
                    self.parent.matrixColAction(g, imag(M)) * 1i;
              case 'C^d -> R^d'
                M = real(self.parent.matrixColAction(g, M));
              otherwise
                M = matrixColAction_double_sparse@replab.Rep(self, g, M);
            end
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = ['Encoding of representation: ' self.type];
        end

        % Rep

        function b = isExact(self)
            b = self.parent.isExact;
        end

        function p = invariantBlocks(self)
            switch self.type
              case {'R^d -> C^d', 'C^d -> R^d'}
                p = self.parent.invariantBlocks;
              otherwise
                error('TODO');
            end
        end

        function b = hasTorusImage(self)
            switch self.type
              case {'R^d -> C^d', 'C^d -> R^d'}
                b = self.parent.hasTorusImage;
              otherwise
                b = false;
            end
        end

        function [torusMap, torusInjection, torusProjection] = torusImage(self)
            switch self.type
              case {'R^d -> C^d', 'C^d -> R^d'}
                [torusMap, torusInjection, torusProjection] = self.parent.torusImage;
              otherwise
                error('Unsupported');
            end
        end

    end

end
