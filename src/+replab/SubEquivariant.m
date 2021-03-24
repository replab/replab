classdef SubEquivariant < replab.Equivariant
% Describes a subspace of an equivariant space induced by subspaces of source/target representations
%
% This equivariant space relates to a parent equivariant space as follows.
%
% - ``self.repR`` must be a `.SubRep` of ``self.parent.repR``,
% - ``self.repC`` must be a `.SubRep` of ``self.parent.repC``.

    properties (SetAccess = protected)
        parent % (`.Equivariant`): Parent equivariant space
    end

    methods (Access = protected)

        function self = SubEquivariant(parent, repR, repC, special)
        % Constructs an equivariant space similar to another equivariant space
        %
        % Args:
        %   parent (`.Equivariant`): Equivariant space of ``repR.parent`` and ``repC.parent``
        %   repR (`.SubRep`): A subrepresentation of ``parent.repR``
        %   repC (`.SubRep`): A subrepresentation of ``parent.repC``
        %   special (charstring): Whether the equivariant subspace has special structure
            self@replab.Equivariant(repR, repC, special);
            assert(isa(parent, 'replab.Equivariant'));
            assert(isa(repR, 'replab.SubRep'));
            assert(isa(repC, 'replab.SubRep'));
            self.parent = parent;
        end

    end

    methods % Implementations

        % Domains

        function l = laws(self)
            l = replab.laws.SubEquivariantLaws(self);
        end

        % Equivariant
        function b = isExact(self)
            b = self.repR.isExact && self.repC.isExact && self.parent.isExact;
        end

    end

    methods (Access = protected) % Implementations

        % Equivariant

        function X1 = project_exact(self, X)
            parentX = self.repR.injection('exact') * X * self.repC.projection('exact');
            X1 = self.projectFromParent(parentX, 'exact');
        end

        function [X1, err] = project_double_sparse(self, X)
            parentX = self.repR.injection('double/sparse') * X * self.repC.projection('double/sparse');
            X1 = self.projectFromParent(parentX, 'double/sparse');
            if nargout > 1
                %if self.repR.isSimilarRep && self.repC.isSimilarRep
                %    c = self.repR.injectionConditionNumberEstimate * self.repC.injectionConditionNumberEstimate;
                %    % error in the first conversion
                %    e1 = replab.numerical.norm2UpperBound(parentX) * self.repR.biorthogonalityErrorBound * c;
                %    % error in the projection
                %    e2 = eParent * c;
                %    % error in the second conversion
                %    e3 = replab.numerical.norm2UpperBound(X1) * self.repC.biorthogonalityErrorBound;
                %    err = e1 + e2 + e3;
                %else
                assert(replab.globals.yolo);
                eR = self.repR.errorBound;
                eC = self.repC.errorBound;
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
                sX = replab.numerical.norm2UpperBound(X);
                err = sX*(eR*cC + cR*eC);
                %end
            end
        end

        function X1 = projectFromParent_exact(self, parentX)
            X1 = self.repR.projection('exact') * self.parent.project(parentX, 'exact') * self.repC.injection('exact');
        end

        function [X1, err] = projectFromParent_double_sparse(self, parentX)
            parentX1 = self.parent.project(parentX, 'double/sparse');
            X1 = self.repR.projection('double/sparse') * parentX1 * self.repC.injection('double/sparse');
            if nargout > 1
                assert(replab.globals.yolo);
                eR = self.repR.errorBound;
                eC = self.repC.errorBound;
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
                sX = replab.numerical.norm2UpperBound(X1);
                err = sX*(eR*cC + cR*eC);
            end
        end

    end

    methods % Projection from parent space

        function [X1, err] = projectFromParent(self, X, type)
        % Projects the given matrix, written in the parent representation space, into the equivariant space of the subrepresentation
        %
        % Raises:
        %   An error if ``type`` is ``'exact'`` and `.isExact` is false.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   X1: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     Projected matrix
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned ``X1`` to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
              case 'exact'
                X1 = self.projectFromParent_exact(X);
                err = 0;
              case 'double'
                if nargout > 1
                    [X1, err] = self.projectFromParent_double_sparse(X);
                else
                    X1 = self.projectFromParent_double_sparse(X);
                end
                X1 = full(X1);
              case 'double/sparse'
                if nargout > 1
                    [X1, err] = self.projectFromParent_double_sparse(X);
                else
                    X1 = self.projectFromParent_double_sparse(X);
                end
              otherwise
                error('Unknown type %s', type);
            end
        end

    end

    methods (Static)

        function E = make(repR, repC, varargin)
        % Returns the space of equivariant linear maps between two subrepresentations
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``repC.image(g) * X = X * repR.image(g)``
        %
        % Args:
        %   repR (`+replab.SubRep`): Representation on the target/row space
        %   repC (`+replab.SubRep`): Representation on the source/column space
        %
        % Keyword Args:
        %   special ('commutant', 'hermitian', 'trivialRows', 'trivialCols' or '', optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double'
        %   parent (`.Equivariant`, optional): Equivariant space of ``repR.parent`` and ``repC.parent``, default: ``[]``
        %
        % Returns:
        %   `+replab.SubEquivariant`: The equivariant vector space
            args = struct('special', '', 'type', 'double', 'parent', []);
            args = replab.util.populateStruct(args, varargin);
            parent = args.parent;
            if isempty(parent)
                switch args.special
                  case 'antilinear'
                    parent = repR.parent.antilinearInvariant(args.type);
                  case 'commutant'
                    parent = repR.parent.commutant(args.type);
                  case 'hermitian'
                    parent = repR.parent.hermitianInvariant(args.type);
                  case 'sesquilinear'
                    parent = repC.parent.sesquilinearInvariant(args.type);
                  case 'trivialRows'
                    parent = repC.parent.trivialRowSpace(args.type);
                    d = repC.dimension;
                    D = repC.parent.dimension;
                    injection = sparse(1:d, 1:d, ones(1, d), D, d);
                    if strcmp('type', 'exact')
                        injection = replab.cyclotomic(injection);
                    end
                    projection = injection';
                    repR = parent.repR.subRep(injection, 'projection', projection);
                  case 'trivialCols'
                    parent = repR.parent.trivialColSpace(args.type);
                    d = repR.dimension;
                    D = repR.parent.dimension;
                    injection = sparse(1:d, 1:d, ones(1, d), D, d);
                    if strcmp('type', 'exact')
                        injection = replab.cyclotomic(injection);
                    end
                    projection = injection';
                    repC = parent.repC.subRep(injection, 'projection', projection);
                  case ''
                    parent = repR.parent.equivariantFrom(repC.parent, 'type', args.type);
                  otherwise
                    error('Invalid special structure');
                end
            end
            E = replab.SubEquivariant(parent, repR, repC, args.special);
        end

    end

end
