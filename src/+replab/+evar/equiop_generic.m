classdef equiop_generic < replab.equiop

    properties (SetAccess = protected)
        f % (function_handle): User defined linear/affine map
        supportsSparse % (logical): Whether the function `.f` handles sparse arguments
    end

    methods

        function self = equiop_generic(group, source, target, sourceInjection, targetInjection, f, supportsSparse)
            assert(isa(f, 'function_handle'));
            self@replab.equiop(group, source, target, sourceInjection, targetInjection);
            self.f = f;
            self.supportsSparse = supportsSparse;
        end

    end

    methods % Implementations

        % equiop

        function Y = apply(self, X)
            f = self.f;
            if isa(X, 'replab.equivar')
                X = value(X);
            end
            if isa(X, 'double')
                if ~self.supportsSparse
                    X = full(X);
                end
                v = f(X);
                Y = replab.equivar(self.target, 'value', v);
                return
            end
            if ~isa(X, 'sdpvar')
                error('Invalid parameter type %s', class(X));
            end
            dec = self.target.decomposition;
            n1 = dec.repR.nComponents;
            n2 = dec.repC.nComponents;
            blocks = cell(n1, n2);
            arg = getbasematrix(X, 0);
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
            ind = getvariables(X);
            for k = 1:length(ind)
                arg = getbasematrix(X, ind(k));
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
            Y = replab.equivar(self.target, 'blocks', blocks);
        end

    end

end
