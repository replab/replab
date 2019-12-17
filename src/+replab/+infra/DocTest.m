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
        %   ps (+replab.+infra.DocTestParseState): Current parse state
        %
        % Returns
        % -------
        %   ps:
        %     +replab.+infra.DocTestParseState: Updated parse state (or [] if parse unsuccessful)
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
        
        function dt = parseDocTest(ps, lineOffset)
        % Parses a doctest, containing multiple command/output pairs
        %
        % Args:
        %   ps (+replab.+infra.DocTestParseState): Current parse state
        %
        % Returns
        % -------
        %   dt:
        %     +replab.+infra.DocTest: The parsed doctest, or [] is unsuccessful
            commands = {};
            outputs = {};
            lineNumbers = [];
            while 1
                [res command output lineNumber] = replab.infra.DocTest.parseCommandOutputPair(ps);
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
            dt = replab.infra.DocTest(lineNumbers + lineOffset, commands, outputs);
        end
        
        function doctests = parse(doc)
        % Finds and parses the doctests in the documentation of an object
        %
        % Args:
        %   doc (.Doc): Documentation
        %
        % Returns:
        %   row cell array of `.DocTest`: The parsed doctests
        %
        % Raises:
        %   A warning if some of the parses are unsuccesful.
            n = doc.nLines;
            content = cell(1, n);
            indent = zeros(1, n);
            for i = 1:n
                l = doc.line(i);
                if isempty(l)
                    indent(i) = 0;
                    content{i} = '';
                else
                    tokens = regexp(l, '^(\s*)(.*)', 'tokens', 'once');
                    if length(tokens) == 1
                        % octave doesn't produce first token if string
                        % doesn't start with some spaces
                        indent(i) = 0;
                        content{i} = strtrim(tokens{1});
                    else
                        indent(i) = length(tokens{1});
                        content{i} = strtrim(tokens{2});
                    end
                end
            end
            doctests = {};
            i = 1;
            while i <= n
                if isequal(content{i}, 'Example:')
                    j = i + 1;
                    while j <= n && (isempty(content{j}) || indent(j) > indent(i))
                        j = j + 1;
                    end
                    ps = replab.infra.DocTestParseState.fromDocTestBlock(content(i+1:j-1));
                    dt = replab.infra.DocTest.parseDocTest(ps, i);
                    if isempty(dt)
                        warning(sprintf('Error while parsing Example: block at line %d'), i);
                    else
                        doctests{1, end+1} = dt;
                    end
                    i = j;
                end
                i = i + 1;
            end
        end
        
    end

end
