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
% The resulting `.equiop` can be called as any MATLAB function. It can work on three types of arguments:
%
% * when called on a ``double`` matrix, it applies the function and returns an `.equivar` in the target space,
% * when called on an ``sdpvar`` matrix, it decomposes the sdpvar and applies the user function on the components of the affine combination,
% * when called on an `.equivar` matrix, it computes the corresponding sdpvar, and applies the user function on the components as well.

    properties (SetAccess = protected)
        source % (`.Equivariant`): Source equivariant space
        target % (`.Equivariant`): Target equivariant space
    end

    methods

        function self = equiop(source, target)
            assert(isa(source, 'replab.Equivariant'));
            assert(isa(target, 'replab.Equivariant'));
            self.source = source;
            self.target = target;
        end

    end

    methods

        function Y = apply(self, X)
        % Computes the application of the operator to a matrix
        %
        % The syntax ``E.apply(X)`` is equivalent to the call ``E(X)``.
        %
        % Args:
        %   X (double(\*,\*) or sdpvar(\*,\*) or `.equivar`): Input matrix
        %
        % Returns:
        %   `.equivar`: Map output
            error('Abstract');
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

    methods (Static)

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
        % Both the source and the target space must be defined over the same group.
        %
        % Args:
        %   source (`.Equivariant`): Source equivariant space
        %   target (`.Equivariant`): Target equivariant space
        %   f (function_handle): User-defined linear or affine map
        %
        % Keyword Args:
        %   supportsSparse (logical, optional): Whether the user-defined function ``f`` may be applied on sparse matrices, default: false
        %
        % Returns:
        %   `.equiop`: Super-operator between equivariant spaces
            args = struct('supportsSparse', false);
            args = replab.util.populateStruct(args, varargin);
            E = replab.evar.equiop_generic(source, target, f, args.supportsSparse);
        end

        function E = restriction(source, target, mu)
        % Returns a linear operator that maps an equivariant space over a group to the equivariant space over a subgroup
        %
        % Args:
        %   source (`.Equivariant`): Source equivariant space
        %   target (`.Equivariant`): Target equivariant space
        %   mu (`.Morphism`): Morphism that embeds ``target.group`` into ``source.group``
            E = replab.evar.equiop_restriction(source, target, mu);
        end

    end

end
