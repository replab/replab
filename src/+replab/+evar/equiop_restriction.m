classdef equiop_restriction < replab.equiop

    properties (SetAccess = protected)
        mu % (`+replab.Morphism`): Morphism from ``target.group`` to ``source.group``
    end

    methods

        function self = equiop_restriction(source, target, mu)
        % Constructs a map between an equivariant space and an equivariant space defined over a restricted representation
            assert(isa(mu, 'replab.Morphism'));
            self@replab.equiop(source, target);
            self.mu = mu;
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.equiop_restrictionLaws(self)
        end

        % equiop

        function Y = apply(self, X)
            if isa(X, 'replab.equivar')
                X = value(X);
            end
            if ~isa(X, 'sdpvar') && ~isa(X, 'double')
                error('Invalid parameter type %s', class(X));
            end
            dec = self.target.decomposition;
            n1 = dec.repR.nComponents;
            n2 = dec.repC.nComponents;
            blocks = cell(n1, n2);
            arg = getbasematrix(X, 0);
            for i = 1:n1
                for j = 1:n2
                    [M, err] = dec.blocks{i,j}.projectAndFactorFromParent(arg);
                    replab.msg(2, 'Error in block (%d,%d), constant factor: %6.2f', i, j, err);
                    blocks{i,j} = M;
                end
            end
            ind = getvariables(X);
            for k = 1:length(ind)
                arg = getbasematrix(X, ind(k));
                for i = 1:n1
                    for j = 1:n2
                        [M, err] = dec.blocks{i,j}.projectAndFactorFromParent(arg);
                        replab.msg(2, 'Error in block (%d,%d), variable %d: %6.2f', i, j, ind(k), err);
                        blocks{i,j} = blocks{i,j} + recover(ind(k))*M;
                    end
                end
            end
            Y = replab.equivar(self.target, 'blocks', blocks);
        end

    end

end
