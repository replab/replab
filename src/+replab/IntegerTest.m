classdef IntegerTest < replab.TotalOrder

    methods
        function l = eqv(self, x, y)
            l = x == y;
        end

        function l = sample(self)
            l = randi([-1000,1000]);
        end

        function c = compare(self,x, y)
            c = sign(x - y);
        end
    end

end
