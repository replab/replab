classdef DocTestStatement < replab.Str
% Describes a command-output pair in a doctest

    properties (SetAccess = protected)
        lineNumber % (integer): Line number (1-based) of the first line of the commandb
        command % (cell(1,\*) of charstring): Command to be evaluated, may be multiline
        output % (cell(1,\*) of charstring): Expected output, may be multiline
        flags % (struct): Flags corresponding to the command evaluation
    end

    methods

        function self = DocTestStatement(lineNumber, command, output, flags)
            self.lineNumber = lineNumber;
            self.command = command;
            self.output = output;
            self.flags = flags;
        end

        function b = isSingleLineCommand(self)
            b = length(self.command) == 1;
        end

        function c = quotedCommand(self)
            c = self.command;
            c = cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], c, 'uniform', 0);
            if length(c) == 1
                c = c{1};
            else
                c = ['strjoin({' strjoin(c, ', ') '}, char(10))'];
            end
        end

        function o = quotedOutput(self)
            o = self.output;
            o = ['{' strjoin(cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], o, 'uniform', 0), ', ') '}'];
        end

        function newDTS = mapLineNumber(self, fun)
        % Returns a new `.DocTestStatement` with the line number mapped according to the given function
        %
        % Args:
        %   fun (function_handle): Function that transforms the line numbers (integer -> integer)
        %
        % Returns:
        %   `.DocTestStatement`: The transformed doctest
            newLN = fun(self.lineNumber);
            newDTS = replab.infra.doctests.DocTestStatement(newLN, self.command, self.output, self.flags);
        end

    end

end
