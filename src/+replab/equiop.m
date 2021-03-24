classdef equiop < replab.Str
% Describes a linear or affine map between equivariant spaces
%
% Formally, in MATLAB, both scalars and vectors are matrices, with one or two dimensions equal to one.
% Thus, `.equiop` describes maps between those different types of objects.
%
% The map is defined by:
%
% * a `.source` equivariant space,
% * a `.target` equivariant space, both must be defined over the same group, which together describe the symmetries,
% * and a user-provided `.f` function handle that should act linearly (or affinely) on its input.
%
% The resulting `.equiop` can be called as any MATLAB function. It can work on three types of arguments:
%
% * when called on a ``double`` matrix, it applies the function and returns an `.equivar` in the target space,
% * when called on an ``sdpvar`` matrix, it decomposes the sdpvar and applies the user function on the components of the affine combination,
% * when called on an `.equivar` matrix, it computes the corresponding sdpvar, and applies the user function on the components as well.
%
% Thus, the user-provided function is always called on arguments of type ``double``. Note that those arguments may be sparse.

    properties (SetAccess = protected)
        source % (`.Equivariant`): Source equivariant space
        target % (`.Equivariant`): Target equivariant space
        f % (function_handle): Affine map
        supportsSparse % (logical): Whether the function `.f` handles sparse arguments
    end

    methods

        function self = equiop(source, target, f, supportsSparse)
        % Constructs an affine map between equivariant spaces
            if nargin < 4 || isempty(supportsSparse)
                supportsSparse = false;
            end
            self.source = source;
            self.target = target;
            self.f = f;
            self.supportsSparse = supportsSparse;
        end

        function n = numArgumentsFromSubscript(self, s, indexingContext)
            n = 1;
        end

        function varargout = subsref(self, s)
            switch s(1).type
              case '.'
                [varargout{1:nargout}] = builtin('subsref', self, s);
              case '()'
                assert(length(s(1).subs) == 1, 'Affine operator must be applied to a single argument');
                arg = s(1).subs{1};
                f = self.f;
                if isa(arg, 'double')
                    if ~self.supportsSparse
                        arg = full(arg);
                    end
                    value = f(arg);
                    varargout{1} = replab.equivar(self.target, 'value', value);
                    return
                end
                if isa(arg, 'sdpvar')
                    sourceValue = arg
                elseif isa(arg, 'replab.equivar')
                    sourceValue = sdpvar(arg);
                else
                    error('Invalid parameter type %s', class(arg));
                end
                dec = self.target.decomposition;
                n1 = dec.repR.nComponents;
                n2 = dec.repC.nComponents;
                blocks = cell(n1, n2);
                arg = getbasematrix(sourceValue, 0);
                if ~self.supportsSparse
                    arg = full(arg);
                end
                v = f(arg);
                for i = 1:n1
                    for j = 1:n2
                        [M, err] = dec.blocks{i,j}.projectAndFactorFromParent(v);
                        replab.msg(2, 'Error in block (%d,%d), constant factor: %6.2f', i, j, err);
                        blocks{i,j} = M;
                    end
                end
                ind = getvariables(sourceValue);
                for k = 1:length(ind)
                    arg = getbasematrix(sourceValue, ind(k));
                    if ~self.supportsSparse
                        arg = full(arg);
                    end
                    v = f(arg);
                    for i = 1:n1
                        for j = 1:n2
                            [M, err] = dec.blocks{i,j}.projectAndFactorFromParent(v);
                            replab.msg(2, 'Error in block (%d,%d), variable %d: %6.2f', i, j, ind(k), err);
                            blocks{i,j} = blocks{i,j} + recover(ind(k))*M;
                        end
                    end
                end
                varargout{1} = replab.equivar(self.target, 'blocks', blocks);
            end
        end

    end

end
