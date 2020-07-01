classdef MixedRadix < replab.Str
% Converts an integer to and from its tabular representation
%
% The digits are represented using doubles, while the integer is of ``vpi`` type.

    properties
        n % (integer): Number of base elements
        base % (integer(1,\*)): Base
    end

    methods

        function self = MixedRadix(base, isOneBased, isBigEndian)
        % Constructs a mixed radix basis
        %
        % Args:
        %   base (integer(1,\*)): Mixed radix basis
        %   isOneBased (logical): Whether digits start at 1 (used for 1-based indexing), must be true
        %   isBigEndian (logical): Whether the digits are written with the most significant digit first, must be true
            assert(isBigEndian);
            assert(isOneBased);
            self.base = base;
            self.n = length(base);
        end

        function sub = sub2ind(self, ind)
        % Return the tabular representation vector corresponding to an integer
        %
        % Args:
        %   index (vpi): Integer
        %
        % Returns:
        %   integer(1,\*): Digits in the mixed radix basis
            sub = zeros(1, self.n);
            base = self.base;
            ind = ind - 1;
            for i = self.n:-1:1
                r = mod(ind, base(i));
                ind = (ind - r)/base(i);
                sub(i) = double(r) + 1;
            end
        end

        function ind = ind2sub(self, sub)
        % Returns the integer corresponding to the given digits
        %
        % Args:
        %   sub (integer(1,\*)): Digits in the mixed radix basis
        %
        % Returns:
        %   vpi: Integer
            ind = vpi(0);
            base = self.base;
            for i = 1:self.n
                ind = ind * base(i);
                ind = ind + vpi(sub(i) - 1);
            end
            ind = ind + vpi(1);
        end

    end

end
