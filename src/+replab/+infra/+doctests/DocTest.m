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

        function n = nStatements(self)
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

        function dt = parse(dtt, errFun)
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

    end

end
