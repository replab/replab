classdef ProgressBar < handle

    properties
        n % integer: Total number of steps
        startTime % date vector: Time of the first step
        times % double(1,:): Seconds elapsed at each step

        barSize % integer: Size of the progress bar in characters
        consoleLine % `.ConsoleLine`: Console state
        largestTimeEstSize % integer: Largest remaining time estimate string seen so far
    end

    methods

        function self = ProgressBar(n)
        % Creates a progress bar
        %
        % Args:
        %   n (integer): Total number of steps
            self.n = n;
            self.consoleLine = replab.infra.repl.ConsoleLine;
            self.barSize = 20; % default value
            self.times = zeros(1, n);
            self.largestTimeEstSize = 0;
        end

        function step(self, i, txt)
        % Updates the progress bar with a 1-based index
        %
        % This method should be called at the beginning of the ``i``-th step
        %
        % Args:
        %   i (integer):
        %   txt (integer, optional): Text corresponding to the current step
            if nargin < 3
                txt = '';
            end
            if i == 1
                self.startTime = clock;
                elapsed = 0;
            else
                elapsed = etime(clock, self.startTime);
            end
            self.times(i) = elapsed;
            if i > 2
                [p, S] = polyfit(1:i, self.times(1:i), 1);
                [last, lastDelta] = polyval(p, self.n, S);
                remaining = last - elapsed;
                timeEst = sprintf('%s remaining', replab.infra.repl.seconds2human(remaining));
            else
                timeEst = '';
            end
            self.largestTimeEstSize = max(self.largestTimeEstSize, length(timeEst));
            timeEstPad = repmat(' ', 1, self.largestTimeEstSize - length(timeEst));
            barFull = round(i/self.n*self.barSize);
            barEmpty = self.barSize - barFull;
            if ~isempty(txt)
                txt = sprintf('   %s', txt);
            end
            str = sprintf('[%s%s] %d/%d, %s%s%s', repmat('#', 1, barFull), repmat('.', 1, barEmpty), ...
                          i, self.n, timeEst, timeEstPad, txt);
            self.consoleLine.update(str);
        end

        function finish(self, finishStr)
        % Finishes the progress bar
        %
        % This method should be called after the last step has been completed
        %
        % Args:
        %   finishStr (charstring, optional): String to replace the progress bar with
        %                                     Default is ``'Task finished'``
            if nargin < 2
                finishStr = [];
            end
            if isempty(finishStr)
                self.consoleLine.update('');
            else
                self.consoleLine.update(finishStr);
                fprintf('\n');
            end
            self.consoleLine = [];
        end

        function log(self, str)
            self.consoleLine.log(str);
        end

        function h = logFunction(self)
        % Returns a function handle that prints things on the console
            h = self.consoleLine.logFunction;
        end

    end

end
