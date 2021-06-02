classdef MixedRadix < replab.Str
% Converts an integer to and from its tabular representation
%
% The digits are represented using doubles, while the integer is of ``vpi`` type.

    properties (SetAccess = protected)
        n % (integer): Number of base elements
        base % (integer(1,\*)): Base
        isOneBased % (logical): Whether digits start at 1
        isBigEndian % (logical): Whether the subindices are written with the most significant subindex first
    end

    methods

        function self = MixedRadix(base, isOneBased, isBigEndian)
        % Constructs a mixed radix basis
        %
        % Args:
        %   base (integer(1,\*)): Mixed radix basis
        %   isOneBased (logical): Whether digits start at 1 (used for 1-based indexing), must be true
        %   isBigEndian (logical): Whether the digits are written with the most significant digit first, must be true
            self.isOneBased = isOneBased;
            self.isBigEndian = isBigEndian;
            self.base = base;
            self.n = length(base);
        end

        function sub = ind2sub(self, ind)
        % Return the tabular representation vector corresponding to an integer
        %
        % Args:
        %   index (vpi): Integer
        %
        % Returns:
        %   integer(1,\*): Digits in the mixed radix basis
            sub = zeros(1, self.n);
            base = self.base;
            if self.isBigEndian
                base = fliplr(base);
            end
            if self.isOneBased
                ind = ind - 1;
            end
            for i = self.n:-1:1
                r = mod(ind, base(i));
                ind = (ind - r)/base(i);
                if self.isOneBased
                    sub(i) = double(r) + 1;
                end
            end
            if self.isBigEndian
                sub = fliplr(sub);
            end
        end

        function ind = sub2ind(self, sub)
        % Returns the integer corresponding to the given digits
        %
        % Args:
        %   sub (integer(1,\*)): Digits in the mixed radix basis
        %
        % Returns:
        %   vpi: Integer
            if self.isOneBased
                sub = sub - 1;
            end
            ind = vpi(0);
            base = self.base;
            if self.isBigEndian
                base = fliplr(base);
                sub = fliplr(sub);
            end
            for i = 1:self.n
                ind = ind * base(i);
                ind = ind + vpi(sub(i));
            end
            if self.isOneBased
                ind = ind + vpi(1);
            end
        end

    end

end
