classdef Settings
    
    methods (Static)
        
        function value = bsgsFailureProbability(newValue)
            persistent BsgsFailureProbability;
            if nargin == 1
                BsgsFailureProbability = newValue;
            elseif isequal(BsgsFailureProbability, [])
                BsgsFailureProbability = 2^-100;
            end
            value = BsgsFailureProbability;
        end
        
        function value = doubleEigTol(newValue)
            persistent DoubleEigTol;
            if nargin == 1
                DoubleEigTol = newValue;
            elseif isequal(DoubleEigTol, [])
                DoubleEigTol = 1e-10;
            end
            value = DoubleEigTol;
        end
            
    end
    
end
