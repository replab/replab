classdef Laws < replab.Str
% Describes the laws that an algebraic structure should obey
%
% Those laws are tested using random instances of the elements in the laws by the testing framework.
%
% Laws are discovered by the `.getTestCases` method, by enumerating the methods of this class.
%
% - Methods that start with a ``law_`` prefix describe a single law. Then the method name has the
%   following structure: ``law_{method name}_{types}``. The part ``{method name}`` is the law name
%   and corresponds to any valid identifier (thus the characters ``A-Z``, ``a-z``, ``0-9`` and ``_`` are allowed).
%   It will be interpreted as the law name, replacing underscore characters by spaces.
%   The ``{types}`` part contains zero or more type identifiers: a type identifier is a letter (``A-Z``, ``a-z``)
%   followed by zero of more digits (``0-9``). This part does not contain any underscore.
%   The number of type identifiers must correspond to the number of arguments of the method (excluding ``self``).
%   Each type identifier must correspond to a property field of this class of type `+replab.Domain`, and
%   those samplable sets will be used to sample the random instances passed to the law check method.
%   Note that a law that does not require any arguments corresponds to a method name ending with an underscore.
%
% - Methods that start with a ``laws_`` prefix must return another `.Laws` instance (see also `+replab.+laws.Collection`).
%   It enables delegation of checks when a tested object has subparts (for example, a `.FiniteGroup` has a `~+replab.FiniteGroup.elements`
%   method of type `.IndexedFamily` that is conveniently checked by `.IndexedFamilyLaws`, see `.FiniteGroupLaws`).
%
% Example:
%    >>> % We build a group from scratch, using function handles,
%    >>> % and verify that it obeys the group laws.
%    >>> n = 10;
%    >>> eqvFun = @(x, y) isequal(x, y);
%    >>> sampleFun = @() randperm(n);
%    >>> composeFun = @(x, y) x(y);
%    >>> identity = 1:n;
%    >>> inverseFun = @(x) arrayfun(@(i) find(x == i), 1:10);
%    >>> S10 = replab.Group.lambda('Permutations of 1..10', eqvFun, sampleFun, composeFun, identity, inverseFun);
%    >>> S10laws = replab.GroupLaws(S10);
%    >>> S10laws.checkSilent
%        1

    properties
        skipSlow = false % (logical): Whether to skip slow tests
    end

    methods

        function thisIsSlow(self)
        % Throws a ``replab:skip`` if `.skipSlow` is true
            if self.skipSlow
                self.skip;
            end
        end

        function self = skippingSlow(self)
        % Sets the ``skipSlow`` property to true and returns the Laws instance
            self.skipSlow = true;
        end

        function assert(self, predicate, context)
        % Assert function with a verbose error message
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
            message = replab.laws.message(errorDesc, context, {predicate}, names, values);
            if replab.compat.isOctave
                error(errorId, '%s', message);
            else
                throwAsCaller(MException(errorId, '%s', message));
            end
        end

        function [testNames testFuns] = getTestCases(self)
        % Enumerates the laws present in this instance by looking for methods with the correct prefix
        %
        % We look for the ``law_`` and ``laws_`` prefixes, see the `.Laws` class description.
            testNames = {};
            testFuns = {};
            M = replab.compat.methodList(metaclass(self));
            [~, I] = sort(cellfun(@(m) m.Name, M, 'uniform', 0));
            M = M(I);
            for i = 1:length(M)
                method = M{i};
                name = method.Name;
                pos = find(name == '_', 1);
                if ~isempty(pos)
                    prefix = name(1:pos);
                    switch prefix
                      case 'laws_'
                        % It is a method that populates test cases from a Laws
                        % subinstance
                        nameParts = strsplit(name, '_');
                        prefix = [strjoin(nameParts(2:end), ' ') '->'];
                        newLaws = self.(name);
                        [newTestNames newTestFuns] = newLaws.getTestCases;
                        newTestNames = cellfun(@(x) [prefix x], newTestNames, 'UniformOutput', false);
                        testNames = horzcat(testNames, newTestNames);
                        testFuns = horzcat(testFuns, newTestFuns);
                      case 'law_'
                        % It is a method that describes a test case
                        [lawName domains] = replab.laws.parseLawMethodName(self, name);
                        if length(domains) > 0
                            % has randomized arguments
                            testFun = @() replab.laws.runNTimes(replab.Laws.nRuns, self, name, domains);
                        else
                            % does not have randomized arguments
                            testFun = @() replab.laws.runNTimes(1, self, name, {});
                        end
                        testNames{end+1} = lawName;
                        testFuns{end+1} = testFun;
                    end
                end
            end
        end

        function res = checkSilent(self)
        % Runs the randomized tests without usign MOxUnit, and returns whether all tests passed
        %
        % Returns:
        %   logical: True if all tests successful
            res = true;
            [~, testFuns] = self.getTestCases;
            for i = 1:length(testFuns)
                f = testFuns{i};
                try
                    f();
                catch
                    res = false;
                end
            end
        end

        function check(self)
        % Runs the randomized tests without using MOxUnit
        %
        % This method is useful from the REPL command line.
        %
        % Example:
        %   >>> S5 = replab.S(5);
        %   >>> L = replab.GroupLaws(S5);
        %   >>> L.check
        %       Checking associativity...
        %       Checking composeAll...
        %       Checking composeN integers...
        %       Checking composeN positive...
        %       Checking composeN zero...
        %       Checking composeWithInverse...
        %       Checking eqv...
        %       Checking identity...
        %       Checking inverse...
        %       Checking inverse compatible with compose...
        %       Checking leftConjugate...
            [testNames testFuns] = self.getTestCases;
            for i = 1:length(testNames)
                disp(sprintf('Checking %s...', testNames{i}));
                f = testFuns{i};
                f();
            end
        end

        function testSuite = addTestCases(self, testSuite, testName)
        % Adds law checks as test cases to the given test suite
        %
        % Args:
        %   testSuite (MOxUnitTestSuite): Adds the test cases of these laws to this
        %   testName (charstring, optional): name of the test
            nRuns = replab.Laws.nRuns;
            if nargin < 3
                testName = class(self);
            end
            % location comes from the calling function
            [call_stack, idx] = dbstack('-completenames');
            location = call_stack(idx + 1).file;
            [testNames testFuns] = self.getTestCases;
            for i = 1:length(testNames)
                caseName = [testName ': ' testNames{i}];
                testCase = MOxUnitFunctionHandleTestCase(caseName, location, testFuns{i});
                testSuite = addTest(testSuite, testCase);
            end
        end

    end

    methods (Static)

        function skip
            errorId = 'replab:skip';
            msg = 'Skipping slow test';
            if replab.compat.isOctave
                error(errorId, msg);
            else
                throwAsCaller(MException(errorId, '%s', msg));
            end
        end

        function out = inexistent(msg)
            errorId = 'replab:inexistent';
            if replab.compat.isOctave
                error(errorId, msg);
            else
                throwAsCaller(MException(errorId, '%s', msg));
            end
            out = [];
        end

        function value = nRuns(newValue)
        % Sets/gets the default number of runs
        %
        % Args:
        %     newValue (integer, optional): Sets the number of runs; if omitted returns the current number of runs
        %
        % Returns:
        %     integer: number of runs
            persistent n
            if isempty(n)
                n = 20;
            end
            if nargin >= 1
                assert((newValue >= 1) && (round(newValue) == newValue), 'newValue should be a positive integer');
                n = newValue;
            end
            value = n;
        end

    end

end
