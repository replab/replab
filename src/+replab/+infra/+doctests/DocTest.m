classdef DocTest < replab.Str

    properties
        parentDoc % `.Doc`: Parent documentation block
        relativeLineNumbers % row vector of integer: 1-based line numbers of first line of the command
        commands % row cell vector of row vector of charstring: commands to be evaluated
        outputs % row cell vector of row vector of charstring: expected output
    end

    methods

        function self = DocTest(parentDoc, relativeLineNumbers, commands, outputs)
            self.parentDoc = parentDoc;
            self.relativeLineNumbers = relativeLineNumbers;
            self.commands = commands;
            self.outputs = outputs;
        end

        function n = nCommands(self)
            n = length(self.commands);
        end

    end

    methods (Static)

        function [ps command output lineNumber] = parseCommandOutputPair(ps)
        % Parses a command/output pair, where the output may be omitted, and each can be multiline
        %
        % Args:
        %   ps (`.ParseState`): Current parse state
        %
        % Returns
        % -------
        %   ps:
        %     `.ParseState`: Updated parse state (or [] if parse unsuccessful)
        %   command:
        %     charstring: Doctest command; if multiline, lines are separated by '\n'
        %   output:
        %     charstring: Doctest output; may be empty, if multiline, lines are separated by '\n'
            [ps newCommand comment lineNumber] = ps.expect('START');
            command = {newCommand};
            if isempty(ps)
                command = {};
                output = {};
                lineNumber = [];
                return
            end
            while 1
                [res newCommand comment] = ps.expect('CONT');
                if isempty(res)
                    break
                else
                    ps = res;
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
            lineNumbers = [];
            while 1
                [res command output lineNumber] = replab.infra.doctests.DocTest.parseCommandOutputPair(ps);
                if isempty(res)
                    break
                else
                    ps = res;
                    commands{1,end+1} = command;
                    outputs{1,end+1} = output;
                    lineNumbers(1,end+1) = lineNumber;
                end
            end
            ps = ps.expect('EOF');
            assert(~isempty(ps));
            dt = replab.infra.doctests.DocTest(doc, lineNumbers + lineOffset, commands, outputs);
        end

    end

end
