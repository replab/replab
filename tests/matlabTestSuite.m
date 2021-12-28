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

        function value = notebooktests(newValue)
            persistent includeNotebooktests;
            if nargin == 1
                includeNotebooktests = newValue;
            elseif isempty(includeNotebooktests)
                includeNotebooktests = true;
            end
            value = includeNotebooktests;
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
            allTestPass = replab_runtests('withCoverage', false, 'doctests', matlabTestSuite.doctests(), 'notebooktests', matlabTestSuite.notebooktests(), 'slowtests', matlabTestSuite.slowtests());
            testCase.verifyEqual(allTestPass, true)
        end
    end
    
end
