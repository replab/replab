classdef TestParameters
% TestParameters
%
% This class stores parameters meant to be accessed by the tests.
%
% See also:
%     TestParameters.withCoverage,
%     TestParameters.fastTestsOnly

    methods (Static)

        function value = withCoverage(newValue)
            % value = withCoverage([newValue])
            % 
            % Sets/tells whether test coverage is active
            %
            % Args:
            %     newValue: boolean, sets the coverage variable if provided
            %         (optional)
            %
            % Returns:
            %     value: whether the coverage variable is true or false

            persistent coverage;
            if isempty(coverage)
                coverage = false;
            end
            if nargin >= 1
                assert(isequal(newValue, false) || isequal(newValue, true), 'newValue should be a boolean');
                coverage = newValue;
            end
            value = coverage;
        end
        
        function value = onlyFastTests(newValue)
            % value = onlyFastTests([newValue])
            % 
            % Sets/tells whether only fast tests should be performed
            %
            % Args:
            %     newValue: boolean, sets the fast test only variable if
            %         provided (optional)
            %
            % Returns:
            %     value: whether the fast test only variable is true or
            %         false

            persistent fastOnly;
            if isempty(fastOnly)
                fastOnly = false;
            end
            if nargin >= 1
                assert(isequal(newValue, false) || isequal(newValue, true), 'newValue should be a boolean');
                fastOnly = newValue;
            end
            value = fastOnly;
        end
        
    end
end
