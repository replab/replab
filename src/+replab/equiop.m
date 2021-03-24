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

    end

    methods % Affine map composition

        function res = andThen(self, applyLast)
        end

        function res = compose(self, applyFirst)
        end

    end

    methods % Multiplication

        function mtimes(lhs, rhs)
        end

    end

end
