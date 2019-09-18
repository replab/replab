classdef Settings
    
    methods (Static)
        
        function value = averagingSamples(newValue)
            persistent AveragingSamples;
            if nargin == 1
                AveragingSamples = newValue;
            elseif isempty(AveragingSamples)
                AveragingSamples = 3;
            end
            value = AveragingSamples;
        end
        
        function value = averagingIterations(newValue)
            persistent AveragingIterations;
            if nargin == 1
                AveragingIterations = newValue;
            elseif isempty(AveragingIterations)
                AveragingIterations = 1000;
            end
            value = AveragingIterations;
        end
        
        function value = useSparse(newValue)
        % Gets/sets the flag whether to use sparse matrices when possible
        %
        % The default value is `false`.
        %
        % Warnings:
        %   This flag should be set once before any computations is made,
        %   as matrix values are cached by RepLAB and those cached values
        %   not updated when this flag changes.
        %
        % Args:
        %   newValue (logical, optional): New flag value
        %
        % Returns:
        %   logical: The current flag value.
            persistent UseSparse;
            if nargin == 1
                UseSparse = newValue;
            elseif isempty(UseSparse)
                UseSparse = false;
            end
            value = UseSparse;
        end

        function value = randomizedSchreierSimsTries(newValue)
        % Gets/sets the number of sifted elements before the BSGS chain is declared complete
        %
        % This is the number of successive failed attempts to generate a new strong generator
        % before deciding the chain is complete in the randomized Schreier-Sims algorithm; 
        % the probability of failure is then less than 1/2^value.
        %
        % Args:
        %   newValue (integer, optional): New value
        %
        % Returns:
        %   integer: The current value.
            persistent RandomizedSchreierSimsTries;
            if nargin == 1
                RandomizedSchreierSimsTries = newValue;
            elseif isempty(RandomizedSchreierSimsTries)
                RandomizedSchreierSimsTries = 1000;
            end
            value = RandomizedSchreierSimsTries;
        end
        
        function value = strMaxColumns(newValue)
            persistent StrMaxColumns;
            if nargin == 1
                StrMaxColumns = newValue;
            elseif isempty(StrMaxColumns)
                StrMaxColumns = 120;
            end
            value = StrMaxColumns;
        end
        
        function value = strMaxRows(newValue)
            persistent StrMaxRows;
            if nargin == 1
                StrMaxRows = newValue;
            elseif isempty(StrMaxRows)
                StrMaxRows = 25;
            end
            value = StrMaxRows;
        end
        
        function value = bsgsFailureProbability(newValue)
            persistent BsgsFailureProbability;
            if nargin == 1
                BsgsFailureProbability = newValue;
            elseif isempty(BsgsFailureProbability)
                BsgsFailureProbability = 2^-100;
            end
            value = BsgsFailureProbability;
        end
        
        function value = doubleEigTol(newValue)
            persistent DoubleEigTol;
            if nargin == 1
                DoubleEigTol = newValue;
            elseif isempty(DoubleEigTol)
                DoubleEigTol = 1e-10;
            end
            value = DoubleEigTol;
        end

        function value = doubleSdpTol(newValue)
            persistent DoubleSdpTol;
            if nargin == 1
                DoubleSdpTol = newValue;
            elseif isempty(DoubleSdpTol)
                DoubleSdpTol = 1e-5;
            end
            value = DoubleSdpTol;
        end
        
    end

end
