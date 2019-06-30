classdef Settings
    
    methods (Static)
        
        function value = strMaxColumns(newValue)
            persistent StrMaxColumns;
            if nargin == 1
                StrMaxColumns = newValue;
            elseif isequal(StrMaxColumns, [])
                StrMaxColumns = 120;
            end
            value = StrMaxColumns;
        end
           
        function value = strMaxRows(newValue)
            persistent StrMaxRows;
            if nargin == 1
                StrMaxRows = newValue;
            elseif isequal(StrMaxRows, [])
                StrMaxRows = 25;
            end
            value = StrMaxRows;
        end
                
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

        function value = doubleSdpTol(newValue)
            persistent DoubleSdpTol;
            if nargin == 1
                DoubleSdpTol = newValue;
            elseif isequal(DoubleSdpTol, [])
                DoubleSdpTol = 1e-6;
            end
            value = DoubleSdpTol;
        end
        
    end
    
end
