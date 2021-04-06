classdef equiop < replab.Obj
% Describes a linear or affine map between equivariant spaces
%
% Formally, in MATLAB, both scalars and vectors are matrices, with one or two dimensions equal to one.
% Thus, `.equiop` describes maps between those different types of objects.
%
% The map is defined between two spaces:
%
% * a `.source` equivariant space,
% * a `.target` equivariant space,
%
% both of which must be defined over related groups, which describe the symmetries.
%
% This `.equiop` is equivariant over a given `.group`, and this group must be a subgroup of
% both ``source.group`` and ``target.group``. The subgroup inclusion is characterized by the injection
% maps ``sourceInjection`` and ``targetInjection``. These must satisfy:
%
% * ``sourceInjection.source == group``
% * ``sourceInjection.target == source.group``
% * ``targetInjection.source == group``
% * ``targetInjection.target == target.group``.
%

% The resulting `.equiop` can be called as any MATLAB function, or through the `.apply` method.

    properties (SetAccess = protected)
        group % (`.CompactGroup`): Map equivariant group
        source % (`.Equivariant`): Source equivariant space
        target % (`.Equivariant`): Target equivariant space
        sourceInjection % (`.Morphism`): Morphism from `.group` to ``source.group``
        targetInjection % (`.Morphism`): Morphism from `.group` to ``target.group``
    end

    methods

        function self = equiop(group, source, target, sourceInjection, targetInjection)
            assert(isa(group, 'replab.CompactGroup'));
            assert(isa(source, 'replab.Equivariant'));
            assert(isa(target, 'replab.Equivariant'));
            assert(isa(sourceInjection, 'replab.Morphism'));
            assert(isa(targetInjection, 'replab.Morphism'));
            assert(sourceInjection.source == group && sourceInjection.target == source.group);
            assert(targetInjection.source == group && targetInjection.target == target.group);
            self.group = group;
            self.source = source;
            self.target = target;
            self.sourceInjection = sourceInjection;
            self.targetInjection = targetInjection;
        end

    end

    methods % Operator application

        function Y = apply(self, X)
        % Computes the application of the operator to a matrix
        %
        % The syntax ``E.apply(X)`` is equivalent to the call ``E(X)``.
        %
        % This method works on three argument types:
        %
        % * when called on a ``double`` matrix, it directly applies the function and factorizes the output,
        % * when called on an ``sdpvar`` matrix, it decomposes the sdpvar and applies the user function on the components of the affine combination,
        % * when called on an `.equivar` matrix, it computes the corresponding ``sdpvar``, and applies the user function on the components as well.
        %
        % In all cases it returns an `.equivar` matrix.
        %
        % Args:
        %   X (double(\*,\*) or sdpvar(\*,\*) or `.equivar`): Input matrix
        %
        % Returns:
        %   `.equivar`: Map output
            error('Abstract');
        end

    end

    methods

        function res = restrict(self, injection)
        % Restricts the equiop to be equivariant under a subgroup of its symmetry group
        %
        % Args:
        %   injection (`.Morphism`): Morphism from the subgroup to `.group`
            res = replab.equiop.generic(self.source, self.target, @(X) X, 'sourceInjection', injection.andThen(self.sourceInjection), 'targetInjection', injection.andThen(self.targetInjection), 'supportsSparse', true);
        end

    end

    methods % Composition

        function res = compose(self, applyFirst)
        % Composition of equivariant operators, the argument applied first
        %
        % Args:
        %   applyFirst (`.equiop`): Morphism to apply first
        %
        % Returns:
        %   `.equiop`: The composition of equiops
            assert(self.group == applyFirst.group);
            res = replab.equiop.generic(applyFirst.source, self.target, @(X) self.apply(applyFirst.apply(X)), 'sourceInjection', applyFirst.sourceInjection, 'targetInjection', self.targetInjection, 'supportsSparse', true);
            % we put supportsSparse true because our function accepts sparse matrices (though the two equiops may convert to full)
        end

        function res = andThen(self, applyLast)
        % Composition of equivariant operators, the argument applied last
        %
        % Args:
        %   applyLast (`.equiop`): Morphism to apply second
        %
        % Returns:
        %   `.equiop`: The composition of equiops
            res = applyLast.compose(self);
        end

    end

    methods % Matlab standard methods

        function n = numel(self)
            n = 1;
        end

        function n = numArgumentsFromSubscript(self, s, indexingContext)
            switch s(1).type
              case '()'
                n = 1;
              otherwise
                n = 0;
            end
        end

        function varargout = subsref(self, s)
        % MATLAB subscripted reference
        %
        % Allows the use of an `.equiop` as if it were a function handle
            switch s(1).type
              case '.'
                [varargout{1:nargout}] = builtin('subsref', self, s);
              case '()'
                assert(length(s(1).subs) == 1, 'Affine operator must be applied to a single argument');
                arg = s(1).subs{1};
                varargout{1} = self.apply(arg);
            end
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.equiopLaws(self);
        end

    end

    methods (Static) % Construction

        function E = generic(source, target, f, varargin)
        % Constructs an `.equiop` from a user-provided function
        %
        % The user-provided function is always called on arguments of type ``double``; when it is applied on
        % variables containing sdpvars, the argument will be expanded on a linear combination of optimization
        % variables, and the user-defined function will be called using basis matrices of type ``double``.
        %
        % The user-defined function may or may not support the use of sparse arguments (for example, if it uses
        % ``reshape/permute`` internally). The ``supportsSparse`` parameter may be set accordingly.
        %
        % The user-provided function must be equivariant over a group, and this group must be a subgroup of
        % both ``source.group`` and ``target.group``. The subgroup inclusion is characterized by the injection
        % maps ``sourceInjection`` and ``targetInjection``, which must satisfy
        % ``group = sourceInjection.source = targetInjection.source``.
        %
        % In the case one or both of ``sourceInjection`` and ``targetInjection`` are missing, they are set to
        % the identity isomorphism of ``source.group`` and ``target.group`` respectively.
        %
        % Args:
        %   source (`.Equivariant`): Source equivariant space
        %   target (`.Equivariant`): Target equivariant space
        %   f (function_handle): User-defined linear or affine map
        %
        % Keyword Args:
        %   sourceInjection (`.Morphism`, optional): An injection from the equiop group to ``source.group``
        %   targetInjection (`.Morphism`, optional): An injection from the equiop group to ``target.group``
        %   supportsSparse (logical, optional): Whether the user-defined function ``f`` may be applied on sparse matrices, default: ``false``
        %
        % Returns:
        %   `.equiop`: Super-operator between equivariant spaces
            args = struct('sourceInjection', [], 'targetInjection', [], 'supportsSparse', false);
            args = replab.util.populateStruct(args, varargin);
            if isempty(args.sourceInjection)
                sourceInjection = replab.Isomorphism.identity(source.group);
            else
                sourceInjection = args.sourceInjection;
            end
            if isempty(args.targetInjection)
                targetInjection = replab.Isomorphism.identity(target.group);
            else
                targetInjection = args.targetInjection;
            end
            group = sourceInjection.source;
            assert(targetInjection.source == group);
            E = replab.evar.equiop_generic(group, source, target, sourceInjection, targetInjection, f, args.supportsSparse);
        end

    end

end
