classdef DocTest < replab.Str
% Describes an executable code sample

    properties (SetAccess = protected)
        statements % (cell(1,\*) of `.DocTestStatement`): Statements parts of this doctest
    end

    methods (Static)

        function d = make(lineNumbers, commands, outputs, flags)
            n = length(lineNumbers);
            statements = cell(1, n);
            for i = 1:n
                statements{i} = replab.infra.doctests.DocTestStatement(lineNumbers(i), ...
                                                                  commands{i}, outputs{i}, flags{i});
            end
            d = replab.infra.doctests.DocTest(statements);
        end

    end

    methods

        function self = DocTest(statements)
            self.statements = statements;
        end

        function l = lineNumbers(self)
            l = cellfun(@(s) s.lineNumber, self.statements);
        end

        function c = commands(self)
            c = cellfun(@(s) s.command, self.statements, 'uniform', 0);
        end

        function o = outputs(self)
            o = cellfun(@(s) s.output, self.statements, 'uniform', 0);
        end

        function f = flags(self)
            f = cellfun(@(s) s.flags, self.statements, 'uniform', 0);
        end

        function b = isSingleLineCommand(self, i)
            b = self.statements{i}.isSingleLineCommand;
        end

        function c = quotedCommand(self, i)
            c = self.statements{i}.quotedCommand;
        end

        function o = quotedOutput(self, i)
            o = self.statements{i}.quotedOutput;
        end

        function n = nStatements(self)
            n = length(self.statements);
        end

        function n = nCommands(self)
            n = length(self.statements);
        end

        function newDT = mapLineNumbers(self, fun)
        % Returns a new `.DocTest` with the line numbers mapped according to the given function
        %
        % Args:
        %   fun (function_handle): Function that transforms the line numbers (integer -> integer)
        %
        % Returns:
        %   `.DocTest`: The transformed doctest
            newStatements = cellfun(@(s) s.mapLineNumber(fun), self.statements, 'uniform', 0);
            newDT = replab.infra.doctests.DocTest(newStatements);
        end

    end

    methods (Static)

        function [ps command output flags lineNumber] = parseCommandOutputPair(ps, errFun)
        % Parses a command/output pair, where the output may be omitted, and each can be multiline
        %
        % Args:
        %   ps (`.ParseState`): Current parse state
        %   errFun (function_handle): Error context display function, see `.parseTests`
        %                             Called as ``errFun(lineNumber)``
        %
        % Returns
        % -------
        %   ps:
        %     `.ParseState`: Updated parse state (or [] if parse unsuccessful)
        %   command:
        %     charstring: Doctest command; if multiline, lines are separated by '\n'
        %   output:
        %     charstring: Doctest output; may be empty, if multiline, lines are separated by '\n'
        %   flags:
        %     struct: Struct containing the doctest flags, if present
        %   lineNumber:
        %     integer: Position of the first command line in the doctest block being parsed
            errId = 'replab:docTestParseError';
            command = [];
            output = [];
            flags = [];
            lineNumber = [];
            [ps newCommand comment lineNumber] = ps.expect('START');
            if isempty(ps)
                % not having a command input/output pair is not an
                % error, just a failed parse
                return
            end
            % but once we start it has to be valid, so anything unexpected is an error
            % (no backtrack)
            command = {newCommand};
            flagsText = regexp(comment, '^\s*doctest:(.+)', 'tokens', 'once');
            if ~isempty(flagsText)
                if iscell(flagsText)
                    flagsText = flagsText{1};
                end
                errFun1 = @() errFun(lineNumber);
                flags = replab.infra.doctests.parseFlags(flagsText, errFun1);
            else
                flags = struct;
            end
            while 1
                [res newCommand comment lineNumber1] = ps.expect('CONT');
                if isempty(res)
                    % No continuation? then we are at the end
                    break
                else
                    ps = res;
                    if ~isempty(regexp(comment, '^\s*doctest:'))
                        errFun(lineNumber1);
                        error(errId, 'Doctest flags can only be present on the first command line');
                    end
                    command{1,end+1} = newCommand;
                end
            end
            output = {};
            while 1
                [res newOutput comment] = ps.expect('OUT');
                if isempty(res)
                    % end of output lines? then we are done
                    break
                else
                    ps = res;
                    output{1,end+1} = newOutput;
                end
            end
        end

        function dt = parse1(dtt, errFun)
        % Parses a doctest from a tokenized doctest block
        %
        % Args:
        %   dtt (`.DocTestTokens`): Tokenized doctest block
        %   errFun (function_handle, optional): Error context display function, see `.parseTests`
        %                                       Called as ``errFun(lineNumber)``
        %
        % Raises:
        %   An error if the parse is unsuccessful
        %
        % Returns:
        %   `.DocTest`: Parsed doctest
            statements = cell(1, 0);
            pos = 1;
            [res, dts] = replab.infra.doctests.DocTestStatement.parse(dtt, pos, errFun);
            while ~isempty(res)
                pos = res;
                statements{1, end+1} = dts;
                [res, dts] = replab.infra.doctests.DocTestStatement.parse(dtt, pos, errFun);
            end
            if dtt.peek(pos) ~= '$'
                errFun(dtt.lineNumbers(pos));
                error('Doctest not ending properly');
            end
            dt = replab.infra.doctests.DocTest(statements);
        end

        function dt = parse(ps, errFun)
        % Parses a doctest, containing multiple command/output pairs
        %
        % Args:
        %   ps (`+replab.+infra.+doctests.ParseState`): Current parse state
        %   errFun (function_handle, optional): Error context display function, see `.parseTests`
        %                                       Called as ``errFun(lineNumber)``
        %
        % Returns
        % -------
        %   dt:
        %     `.DocTest`: The parsed doctest, with line numbers corresponding to
        %                 position in the parsed block
        %
        % Raises:
        %   An error if the parse is unsuccessful
            errId = 'replab:docTestParseError';
            if nargin < 2
                errFun = @(l) fprintf('Error in line %d\n', l);
            end
            commands = {};
            outputs = {};
            flags = {};
            lineNumbers = [];
            while 1
                [res command output flagss lineNumber] = replab.infra.doctests.DocTest.parseCommandOutputPair(ps, errFun);
                if isempty(res)
                    break
                else
                    ps = res;
                    commands{1,end+1} = command;
                    outputs{1,end+1} = output;
                    flags{1,end+1} = flagss;
                    lineNumbers(1,end+1) = lineNumber;
                end
            end
            [res, ~, ~, lineNumber] = ps.expect('EOF');
            if isempty(res)
                errFun(lineNumber);
                error(errId, 'End of file expected');
            end
            dt = replab.infra.doctests.DocTest.make(lineNumbers, commands, outputs, flags);
        end

    end

end
