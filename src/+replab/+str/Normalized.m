classdef Normalized < replab.Str
% Represents a vector along with a normalization factor for pretty printing
%
% The vector value is ``vec * factor``, except that `factor` is already represented by a string
% to be concatenated with the vector string representation
    
    properties (SetAccess = protected)
        vec % double row vector: Vector to print
        factor % char: String description of the normalization factor
    end
    
    methods
        
        function self = Normalized(vec, factor)
            assert(isa(factor, 'char'), 'Factor must be  a string');
            self.vec = vec;
            self.factor = factor;
        end
        
        function s = shortStr(self, maxColumns)
            maxL = maxColumns - length(self.factor);
            vecStr = replab.shortStr(self.vec, maxL);
            s = [vecStr self.factor];
        end
        
    end
    
end
