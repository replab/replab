classdef ConsoleLine < handle
% Prints strings on the same console line

    properties
        lineContent % String currently displayed
    end

    methods

        function self = ConsoleLine
            self.lineContent = '';
        end

        function update(self, str)
        % Updates the displayed line
            n1 = length(self.lineContent);
            % truncate the displayed line if terminal is not big enough
            [nR, nC] = replab.compat.terminalSize;
            nC = nC - 3;
            if ~isempty(nC) && length(str) > nC-2
                str = [str(1:nC-5), '...'];
            end
            n2 = length(str);
            backslash = char(8);
            if n1 > 0
                fprintf('%s', repmat(backslash, 1, n1));
            end
            nErased = max(0, n1 - n2);
            fprintf('%s%s%s', str, repmat(' ', 1, nErased), repmat(backslash, 1, nErased));
            self.lineContent = str;
        end

        function log(self, str)
        % Prints the given string (must be line ending terminated), and reprints the current line afterwards
            lastContent = self.lineContent;
            self.update('');
            fprintf('%s', str);
            fprintf('%s', lastContent);
            self.lineContent = lastContent;
        end

        function h = logFunction(self)
        % Returns a function handle that prints things on the console
            h = @(str) self.log(str);
        end

    end

end
