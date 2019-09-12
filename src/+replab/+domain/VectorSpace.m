classdef VectorSpace < replab.Domain
% Base class for vector spaces that can be either real or complex
    
    properties
        field % {'R', 'C'}: Matrices with real (R) or complex (C) coefficients
    end
    
    methods
        
        function b = overR(self)
        % Returns whether this vector space is defined over the real field
        %
        % Returns:
        %   logical: True if the vector space is defined over the real field
            b = isequal(self.field, 'R');
        end
        
        function b = overC(self)
        % Returns whether this vector space is defined over the complex field
        %
        % Returns:
        %   logical: True if this vector space are defined over the complex field
            b = isequal(self.field, 'C');
        end
        
    end
    
end
