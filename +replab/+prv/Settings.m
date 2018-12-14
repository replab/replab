classdef Settings
    
    methods (Static)
       
        function value = doubleEigTol(newValue)
            persistent DoubleEigTol;
            if nargin == 1
                DoubleEigTol = newValue;
            elseif isequal(DoubleEigTol, [])
                    DoubleEigTol = 1e-9;
            end
            value = DoubleEigTol;
        end
            
    end
    
end
