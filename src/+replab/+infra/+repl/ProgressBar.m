classdef ProgressBar < handle

    properties
        n % integer: Total number of steps
        startTime % date vector: Time of the first step
        steps % (integer(1,\*)): Steps with a sample
        times % (double(1,\*)): Seconds elapsed at each step

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
            self.steps = zeros(1, 0);
            self.times = zeros(1, 0);
            self.largestTimeEstSize = 0;
            self.startTime = clock;
        end

        function step(self, i, txt)
        % Updates the progress bar with a 1-based index
        %
        % This method should be called at the beginning of the ``i``-th step.
        %
        % Args:
        %   i (integer): Current step
        %   txt (integer, optional): Text corresponding to the current step
            if nargin < 3
                txt = '';
            end
            elapsed = etime(clock, self.startTime);
            if (length(self.steps) >= 1) && (self.steps(end) >= i)
                % Remove progresses larger than current one
                lastSmaller = find(self.steps >= i, 1, 'first') - 1;
                self.steps = self.steps(1:lastSmaller);
                self.times = self.times(1:lastSmaller);
            end
            self.steps = [self.steps i];
            self.times(1,end+1) = elapsed;
            if length(self.steps) > 1
                %[p, S] = polyfit(double(self.steps)/double(self.n), self.times, 1);
                %[last, lastDelta] = polyval(p, 1, S);
                p = polyfit(double(self.steps)/double(self.n), self.times, 1);
                last = polyval(p, 1);
                remaining = last - elapsed;
                timeEst = sprintf('%s remaining', replab.infra.repl.seconds2human(remaining));
            else
                timeEst = '';
            end
            self.largestTimeEstSize = max(self.largestTimeEstSize, length(timeEst));
            timeEstPad = repmat(' ', 1, self.largestTimeEstSize - length(timeEst));
            barFull = min(round(double(i)/double(self.n)*self.barSize), self.barSize);
            barEmpty = self.barSize - barFull;
            if ~isempty(txt)
                txt = sprintf('   %s', txt);
            end
            str = sprintf('[%s%s] %s/%s, %s%s%s', repmat('#', 1, barFull), repmat('.', 1, barEmpty), ...
                          strtrim(num2str(i)), strtrim(num2str(self.n)), timeEst, timeEstPad, txt);
            self.consoleLine.update(str);
        end

        function stepNoTimeEstimation(self, i, txt)
        % Updates the progress bar with a 1-based index
        %
        % This method should be called at the beginning of the ``i``-th step.
        % It prints the advancement but does not provide a time estimate
        % for completion.
        %
        % Args:
        %   i (integer): Current step
        %   txt (integer, optional): Text corresponding to the current step
            if nargin < 3
                txt = '';
            end
            elapsed = etime(clock, self.startTime);
            if (length(self.steps) >= 1) && (self.steps(end) >= i)
                % Remove progresses larger than current one
                lastSmaller = find(self.steps >= i, 1, 'first') - 1;
                self.steps = self.steps(1:lastSmaller);
                self.times = self.times(1:lastSmaller);
            end
            self.steps = [self.steps i];
            self.times(1,end+1) = elapsed;
            timeEst = '';
            timeEstPad = repmat(' ', 1, self.largestTimeEstSize);
            barFull = min(round(double(i)/double(self.n)*self.barSize), self.barSize);
            barEmpty = self.barSize - barFull;
            if ~isempty(txt)
                txt = sprintf('   %s', txt);
            end
            str = sprintf('[%s%s] %s/%s, %s%s%s', repmat('#', 1, barFull), repmat('.', 1, barEmpty), ...
                          strtrim(num2str(i)), strtrim(num2str(self.n)), timeEst, timeEstPad, txt);
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
