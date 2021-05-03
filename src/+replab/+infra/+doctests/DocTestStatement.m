classdef DocTestStatement < replab.Str
% Describes a command-output pair in a doctest

    properties (SetAccess = protected)
        lineNumber % (integer): Line number (1-based) of the first line of the command
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

    methods (Static)

        function [pos, command, lineNumber1, flags] = parseCommand(dtt, pos, errFun)
            errId = 'replab:docTestParseError';
            [pos, tag, command, comment1, lineNumber1] = dtt.take(pos);
            switch tag
              case '$'
                pos = [];
                return
              case {'o', 'O'}
                errFun(lineNumber1);
                error('errId', 'A command must start with >>>');
              case 's'
                % do nothing
              case 'S'
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
                        command = [command line];
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

        function [pos, output] = parseOutput(dtt, pos, errFun)
            errId = 'replab:docTestParseError';
            [pos, tag, output, ~, lineNumber1] = dtt.take(pos);
            switch tag
              case '$'
                pos = [];
                return
              case {'s', 'S'}
                pos = [];
                return
              case 'o'
                % do nothing
              case 'O'
                while 1
                    [pos, tag1, line, ~, ln] = dtt.take(pos);
                    switch tag1
                      case {'s', 'S'}
                        errFun(ln);
                        error(errId, 'A line ending with ... cannot be followed by a line starting with >>>.');
                      case 'o'
                        output = [output line];
                        break
                      case 'O'
                        output = [output line];
                      case '$'
                        errFun(ln);
                        error(errId, 'A line ending with ... must be followed by a continuation line');
                    end
                end
            end
        end

        function [pos, dts] = parse(dtt, pos, errFun)
            [pos, command, lineNumber, flags] = replab.infra.doctests.DocTestStatement.parseCommand(dtt, pos, errFun);
            if isempty(pos)
                return
            end
            [res, output] = replab.infra.doctests.DocTestStatement.parseOutput(dtt, pos, errFun);
            if isempty(res)
                output = '';
            else
                pos = res;
            end
            dts = replab.infra.doctests.DocTestStatement(lineNumber, command, output, flags);
        end

    end

end
