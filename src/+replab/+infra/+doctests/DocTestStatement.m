classdef DocTestStatement < replab.Str
% Describes a command-output pair in a doctest

    properties (SetAccess = protected)
        lineNumber % (integer): Line number (1-based) of the first line of the command
        command % (charstring): Command to be evaluated
        commandType % (silent, assign or repl): Type of command
        expectedStandardOutput % (cell(1,\*) of charstring): Lines of expected standard output
        variableNames % (cell(1,\*) of charstring): Variable names
        variableValues % (cell(1,\*) of cell(1,\*) of charstring): Variable values
        output % (cell(1,\*) of charstring): Expected output, may be multiline
        flags % (struct): Flags corresponding to the command evaluation
    end

    methods

        function self = DocTestStatement(lineNumber, command, commandType, output, expectedStandardOutput, variableNames, variableValues, flags)
            assert(ischar(command));
            self.lineNumber = lineNumber;
            self.command = command;
            self.commandType = commandType;
            self.expectedStandardOutput = expectedStandardOutput;
            self.variableNames = variableNames;
            self.variableValues = variableValues;
            self.output = output;
            self.flags = flags;
        end

        function write(self, fid)
        % Writes the statement in a doctest
        %
        % Args:
        %   fid (integer): File handle
            quote = @(s) ['''' strrep(s, '''', '''''') ''''];
            quotes = @(S) ['{' strjoin(cellfun(quote, S, 'uniform', 0), ', ') '}'];
            quotedCommand = quote(self.command);
            switch self.commandType
              case 'silent'
                fprintf(fid, '  replout_ = evalc(%s);\n', quote(self.command));
                fprintf(fid, '  assertEqualStandardOutput(replout_, %s, filename, %d);\n', quotes(self.expectedStandardOutput), self.lineNumber);
              case 'repl'
                fprintf(fid, '  [replout_, replans_] = evalc(%s);\n', quote(self.command));
                fprintf(fid, '  assertEqualStandardOutput(replout_, %s, filename, %d);\n', quotes(self.expectedStandardOutput), self.lineNumber);
                fprintf(fid, '  assertEqualValue(%s, replans_, %s, filename, %d);\n', quote(self.variableNames{1}), quotes(self.variableValues{1}), self.lineNumber);
              case 'assign'
                fprintf(fid, '  replout_ = evalc(%s);\n', quote([self.command ';']));
                fprintf(fid, '  assertEqualStandardOutput(replout_, %s, filename, %d);\n', quotes(self.expectedStandardOutput), self.lineNumber);
                for i = 1:length(self.variableNames)
                    fprintf(fid, '  assertEqualValue(%s, %s, %s, filename, %d);\n', quote(self.variableNames{i}), self.variableNames{i}, quotes(self.variableValues{i}), self.lineNumber);
                end
            end
        end

        function c = quotedCommand(self, addSemicolon)
        % Returns the command sandwiched between single quotes, with quotes escaped
            if nargin < 2
                addSemicolon = false;
            end
            c = self.command;
            if addSemicolon
                c = ['''' strrep(c, '''', '''''') ';'''];
            else
                c = ['''' strrep(c, '''', '''''') ''''];
            end
        end

        function o = quotedOutput(self)
        % Returns the output lines expressed as a Matlab expression of cell type
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
            newDTS = replab.infra.doctests.DocTestStatement(newLN, self.command, self.commandType, self.output, self.expectedStandardOutput, self.variableNames, self.variableValues, self.flags);
        end

    end

    methods (Static)

        function [pos, command, lineNumber1, flags] = parseCommand(dtt, pos, errFun)
            errId = 'replab:docTestParseError';
            [pos, tag, command, comment1, lineNumber1] = dtt.take(pos);
            switch tag
              case '$'
                pos = [];
                flags = [];
                return
              case {'o', 'O'}
                errFun(lineNumber1);
                error('errId', 'A command must start with >>>');
              case 's'
                % do nothing
              case 'S'
                command = command(1:end-3);
                % collect command lines
                while 1
                    [pos, tag1, line, ~, ln] = dtt.take(pos);
                    switch tag1
                      case {'s', 'S'}
                        errFun(ln);
                        error(errId, 'A line ending with ... cannot be followed by a line starting with >>>.');
                      case 'o'
                        command = [command line];
                        break
                      case 'O'
                        command = [command line(1:end-3)];
                      case '$'
                        errFun(ln);
                        error(errId, 'A line ending with ... must be followed by a continuation line');
                    end
                end
            end
            ind = strfind(comment1, 'doctest:');
            switch length(ind)
              case 0
                flags = struct;
              case 1
                flags = replab.infra.doctests.parseFlags(comment1(ind+8:end), @() errFun(lineNumber1));
              otherwise
                errFun(lineNumber1);
                error(errId, 'Invalid flags comment in doctest');
            end
        end

        function [pos, outputLines] = parseOutput(dtt, pos, errFun)
            errId = 'replab:docTestParseError';
            outputLines = cell(1, 0);
            [res, tag, line] = dtt.take(pos);
            while any(tag == 'oO')
                outputLines{1,end+1} = strtrim(line);
                pos = res;
                [res, tag, line] = dtt.take(pos);
            end
        end

        function [t, outs] = identifyCommand(c)
        % Identifies the type of a command, and the variable names in output
        %
        % Args:
        %   c (charstring): Command
        %
        % Returns
        % -------
        %   t: silent, assign, repl
        %     Command type
        %   outs: char(1,\*) of charstring
        %     Names of variables
            if isempty(c) || c(end) == ';'
                t = 'silent';
                outs = cell(1, 0);
            else
                T = regexp(c, '^\[(\s*(\w(\w|\d)*)(\s*,\s*(\w(\w|\d)*))*)\s*\]\s*=[^=]', 'tokens');
                if ~isempty(T)
                    t = 'assign';
                    outs = cellfun(@strtrim, strsplit(T{1}{1}, ','), 'uniform', 0);
                    return
                end
                T = regexp(c, '^(\w(\w|\d)*)\s*=[^=]', 'tokens');
                if ~isempty(T)
                    t = 'assign';
                    outs = {strtrim(T{1}{1})};
                    return
                end
                outs = {'ans'};
                t = 'repl';
            end
        end

        function [pos, dts] = parse(dtt, pos, errFun)
            [pos, command, lineNumber, flags] = replab.infra.doctests.DocTestStatement.parseCommand(dtt, pos, errFun);
            if isempty(pos)
                dts = [];
                return
            end
            [res, output] = replab.infra.doctests.DocTestStatement.parseOutput(dtt, pos, errFun);
            if isempty(res)
                output = '';
            else
                pos = res;
            end
            [commandType, variableNames] = replab.infra.doctests.DocTestStatement.identifyCommand(command);
            switch commandType
              case 'silent'
                expectedStandardOutput = output;
                variableNames = cell(1, 0);
                variableValues = cell(1, 0);
              case 'repl'
                [~, ind] = ismember('ans =', output);
                assert(length(ind) <= 1);
                if isempty(ind)
                    expectedStandardOutput = cell(1, 0);
                    variableValues = {output};
                else
                    expectedStandardOutput = output(1:ind-1);
                    variableValues = {output(ind+1:end)};
                end
              case 'assign'
                [~, ind] = ismember([variableNames{1} ' ='], output);
                if length(variableNames) == 1
                    assert(length(ind) <= 1);
                    if isempty(ind)
                        expectedStandardOutput = cell(1, 0);
                        variableValues = {output};
                    else
                        expectedStandardOutput = output(1:ind-1);
                        variableValues = {output(ind+1:end)};
                    end
                else
                    assert(length(ind) == 1);
                    expectedStandardOutput = output(1:ind-1);
                    prevInd = ind;
                    variableValues = cell(1, length(variableNames));
                    for i = 2:length(variableNames)
                        [~, ind] = ismember([variableNames{i} ' ='], output);
                        assert(length(ind) == 1 && ind > prevInd);
                        variableValues{i-1} = output(prevInd+1:ind-1);
                    end
                    variablesValues{length(variableNames)} = output(ind+1:end);
                end
            end
            dts = replab.infra.doctests.DocTestStatement(lineNumber, command, commandType, output, expectedStandardOutput, variableNames, variableValues, flags);
        end

    end

end
