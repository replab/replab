classdef DocTest < replab.Str
% Describes an executable code sample
    properties
        lineNumbers % row vector of integer: 1-based line numbers of first line of the command
                    %
                    %                        (interpretation of those line numbers depend on context)
        commands % row cell vector of row vector of charstring: commands to be evaluated
        outputs % row cell vector of row vector of charstring: expected output
        flags % row cell vector of struct: flags corresponding to the command evaluation
    end

    methods

        function self = DocTest(lineNumbers, commands, outputs, flags)
            self.lineNumbers = lineNumbers;
            self.commands = commands;
            self.outputs = outputs;
            self.flags = flags;
        end

        function b = isSingleLineCommand(self, i)
            b = length(self.commands{i}) == 1;
        end

        function c = quotedCommand(self, i)
            c = self.commands{i};
            c = cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], c, 'uniform', 0);
            if length(c) == 1
                c = c{1};
            else
                c = ['strjoin({' strjoin(c, ', ') '}, char(10))'];
            end
        end

        function o = quotedOutput(self, i)
            o = self.outputs{i};
            o = ['{' strjoin(cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], o, 'uniform', 0), ', ') '}'];
        end

        function n = nCommands(self)
            n = length(self.commands);
        end

        function newDT = mapLineNumbers(self, fun)
        % Returns a new `.DocTest` with the line numbers mapped according to the given function
        %
        % Args:
        %   fun (function_handle): Function that transforms the line numbers (integer -> integer)
        %
        % Returns:
        %   `.DocTest`: The transformed doctest
            newLineNumbers = arrayfun(fun, self.lineNumbers);
            newDT = replab.infra.doctests.DocTest(newLineNumbers, ...
                                                  self.commands, self.outputs, self.flags);
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

        function dt = parse(ps, errFun)
        % Parses a doctest, containing multiple command/output pairs
        %
        % Args:
        %   ps (`replab.infra.doctests.ParseState`): Current parse state
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
            dt = replab.infra.doctests.DocTest(lineNumbers, commands, outputs, flags);
        end

    end

end
