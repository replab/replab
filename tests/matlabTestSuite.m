classdef matlabTestSuite < matlab.unittest.TestCase
% To call the test suits through matlab's test interface, here without

    methods (Static)
        function value = doctests(newValue)
            persistent includeDoctests;
            if nargin == 1
                includeDoctests = newValue;
            elseif isempty(includeDoctests)
                includeDoctests = true;
            end
            value = includeDoctests;
        end

        function value = notebooks(newValue)
            persistent includeNotebooks;
            if nargin == 1
                includeNotebooks = newValue;
            elseif isempty(includeNotebooks)
                includeNotebooks = true;
            end
            value = includeNotebooks;
        end

        function value = slowtests(newValue)
            persistent includeSlowtests;
            if nargin == 1
                includeSlowtests = newValue;
            elseif isempty(includeSlowtests)
                includeSlowtests = true;
            end
            value = includeSlowtests;
        end
    end
    
    methods (Test)
        function testAll(testCase)
            allTestPass = replab_runtests('withCoverage', false, 'doctests', matlabTestSuite.doctests(), 'notebooks', matlabTestSuite.notebooks(), 'slowtests', matlabTestSuite.slowtests());
            testCase.verifyEqual(allTestPass, true)
        end
    end
    
end
