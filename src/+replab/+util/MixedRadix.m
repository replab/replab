classdef MixedRadix < replab.Str
% Converts an integer to and from its tabular representation
%
% The digits are represented using doubles, while the integer is of ``vpi`` type.

    properties (SetAccess = protected)
        n % (integer): Number of base elements
        base % (integer(1,\*)): Base
        oneBased % (logical): Whether digits start at 1
        bigEndian % (logical): Whether the subindices are written with the most significant subindex first
    end

    methods

        function self = MixedRadix(base, varargin)
        % Constructs a mixed radix basis
        %
        % Args:
        %   base (integer(1,\*)): Mixed radix basis
        %
        % Keyword Args:
        %   oneBased (logical): Whether digits start at 1 (used for 1-based indexing), default: true
        %   bigEndian (logical): Whether the digits are written with the most significant digit first, default: true
            args = struct('oneBased', true, 'bigEndian', true);
            args = replab.util.populateStruct(args, varargin);
            self.oneBased = args.oneBased;
            self.bigEndian = args.bigEndian;
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
            if self.bigEndian
                base = fliplr(base);
            end
            if self.oneBased
                ind = ind - 1;
            end
            for i = self.n:-1:1
                r = mod(ind, base(i));
                ind = (ind - r)/base(i);
                if self.oneBased
                    sub(i) = double(r) + 1;
                end
            end
            if self.bigEndian
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
            if self.oneBased
                sub = sub - 1;
            end
            ind = vpi(0);
            base = self.base;
            if self.bigEndian
                base = fliplr(base);
                sub = fliplr(sub);
            end
            for i = 1:self.n
                ind = ind * base(i);
                ind = ind + vpi(sub(i));
            end
            if self.oneBased
                ind = ind + vpi(1);
            end
        end

    end

end
