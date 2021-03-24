classdef equiop < replab.Str
% Affine map between equivariant spaces

    properties (SetAccess = protected)
        source % (`.Equivariant`): Source equivariant space
        target % (`.Equivariant`): Target equivariant space
        f % (function_handle): Affine map
    end

    methods

        function self = equiop(source, target, f)
        % Constructs an affine map between equivariant spaces
            self.source = source;
            self.target = target;
            self.f = f;
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
                    value = f(arg);
                    varargout{1} = replab.equivar(self.target, 'value', value);
                elseif isa(arg, 'replab.equivar')
                    sourceValue = sdpvar(arg);
                    dec = self.target.decomposition;
                    n1 = dec.repR.nComponents;
                    n2 = dec.repC.nComponents;
                    blocks = cell(n1, n2);
                    v = f(getbasematrix(sourceValue, 0));
                    for i = 1:n1
                        for j = 1:n2
                            [M, err] = dec.blocks{i,j}.projectAndFactorFromParent(v);
                            replab.msg(2, 'Error in block (%d,%d), constant factor: %6.2f', i, j, err);
                            blocks{i,j} = M;
                        end
                    end
                    ind = getvariables(sourceValue);
                    for k = 1:length(ind)
                        v = f(getbasematrix(sourceValue, ind(k)));
                        for i = 1:n1
                            for j = 1:n2
                                [M, err] = dec.blocks{i,j}.projectAndFactorFromParent(v);
                                replab.msg(2, 'Error in block (%d,%d), variable %d: %6.2f', i, j, ind(k), err);
                                blocks{i,j} = blocks{i,j} + recover(ind(k))*M;
                            end
                        end
                    end
                    varargout{1} = replab.equivar(self.target, 'blocks', blocks);
                else
                    error('Invalid operator parameter of class %s', class(arg));
                end
            end
        end

    end

end
