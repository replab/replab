classdef DocTest < replab.Str
% Describes an executable code sample
    properties
        lineNumbers % row vector of integer: 1-based line numbers of first line of the command
                    %                        (relative to the parsed comment block)
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

        function n = nCommands(self)
            n = length(self.commands);
        end

    end

    methods (Static)

        function [ps command output flags lineNumber] = parseCommandOutputPair(ps, errFun)
        % Parses a command/output pair, where the output may be omitted, and each can be multiline
        %
        % Args:
        %   ps (`.ParseState`): Current parse state
        %   errFun (function_handle): Error function, see `.parseTests`
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
            command = [];
            output = [];
            flags = [];
            lineNumber = [];
            [ps newCommand comment lineNumber] = ps.expect('START');
            if isempty(ps)
                return
            end
            command = {newCommand};
            flagsText = regexp(comment, '^\s*doctest:(.+)', 'tokens', 'once');
            if ~isempty(flagsText)
                if iscell(flagsText)
                    flagsText = flagsText{1};
                end
                flags = replab.infra.doctests.parseFlags(flagsText);
            else
                flags = struct;
            end
            while 1
                [res newCommand comment] = ps.expect('CONT');
                if isempty(res)
                    break
                else
                    ps = res;
                    assert(isempty(regexp(comment, '^\s*doctest:')), ...
                           'Doctest flags can only be present on the first command line');
                    command{1,end+1} = newCommand;
                end
            end
            output = {};
            while 1
                [res newOutput comment] = ps.expect('OUT');
                if isempty(res)
                    break
                else
                    ps = res;
                    output{1,end+1} = newOutput;
                end
            end
        end

        function dt = parseDocTest(doc, ps, lineOffset)
        % Parses a doctest, containing multiple command/output pairs
        %
        % Args:
        %   doc (`replab.infra.Doc`): Doc comment block in which this is contained
        %   ps (`replab.infra.doctests.ParseState`): Current parse state
        %   lineOffset (integer): Line number of the ``Example:`` directive
        %
        % Returns
        % -------
        %   dt:
        %     `.DocTest`: The parsed doctest, or [] is unsuccessful
            commands = {};
            outputs = {};
            flags = {};
            lineNumbers = [];
            while 1
                [res command output flagss lineNumber] = replab.infra.doctests.DocTest.parseCommandOutputPair(ps);
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
            ps = ps.expect('EOF');
            assert(~isempty(ps));
            dt = replab.infra.doctests.DocTest(doc, lineNumbers + lineOffset, commands, outputs, flags);
        end

    end

end
