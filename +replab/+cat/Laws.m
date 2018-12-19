classdef Laws < replab.Str
% Describes a structure that obeys to algebraic laws
%
% Laws are described by methods called law_METHODNAME_TYPES
% where
% - METHODNAME is a snake_case identifier with each word starting
%   with a lower case character,
% - TYPES is a string composed of the concatenation of the possible definitions below
%   - D is an element of the current domain
%   - Zn is an integer from -n to n, represented as a double
%   - Nn is an integer from 1 to n, represented as a double
%   - N0n is an integer from 0 to n, represented as a double
%   - any other letter should be the name of a property of this Laws instance
%     representing a domain
    
    methods
        
        function assertTrue(self, predicate, context)
            if ~isscalar(predicate) || ~islogical(predicate)
                errorDesc = 'input %s is not a logical scalar';
                errorId = 'assertTrue:invalidCondition';
            elseif ~predicate
                errorDesc = 'input %s does not evaluate to true';
                errorId = 'assertTrue:falseCondition';
            else
                return
            end
            
            if nargin < 3
                context = '';
            end
            
            names = evalin('caller', 'who');
            nV = length(names);
            values = cell(1, nV);
            for i = 1:nV
                values{i} = evalin('caller', names{i});
            end
            
            message = replab.cat.lawsMessage(errorDesc, context, {predicate}, names, values);
            if moxunit_util_platform_is_octave()
                error(errorId, '%s', message);
            else
                throwAsCaller(MException(errorId, '%s', message));
            end
        end
                
        function testSuite = lawsAddTestCases(self, testSuite, varargin)
        % Adds law checks as test cases to the given test suite
        %
        % Inputs:
        % testSuite      - MOxUnitTestSuite object
        % varargin       - name/value pairs such as
        %                  {'nRuns', integer} changes the number of runs
        %                  for random checks (default: 100)
            nRuns = 100;
            testName = class(self);
            for i = 1:2:length(varargin)
                argName = varargin{i};
                argValue = varargin{i + 1};
                switch argName
                  case 'name'
                    testName = argValue;
                  case 'nRuns'
                    nRuns = argValue;
                  otherwise
                    error(sprintf('Argument %s unknown', argName));
                end
            end
            % test case name comes from the calling function
            [call_stack, idx] = dbstack('-completenames');
            location = call_stack(idx + 1).file;
            % look for law_propertyName_TYPES methods in self
            MC = metaclass(self);
            M = MC.MethodList;
            nRuns = 100;
            for i = 1:length(M)
                if isa(M, 'cell')
                    % Octave returns a cell array 
                    method = M{i};
                else
                    % Matlab an object array
                    method = M(i);
                end
                name = method.Name;
                if length(name) > 4 && isequal(name(1:4), 'law_')
                    [lawName isRandom argFuns] = replab.cat.collectTypes(self, name);
                    testFun = @() replab.cat.runNTimes(nRuns, self, name, argFuns);
                    caseName = [testName ': ' lawName];
                    testCase = MOxUnitFunctionHandleTestCase(caseName, location, testFun);
                    testSuite = addTest(testSuite, testCase);
                end
            end
        end
        
    end
    
end
