function runNTimes(nRuns, laws, methodName, sets)
% Runs the given law check a given number of times
%
% Args:
%   nRuns (integer): Number of runs
%   laws (`+replab.Laws`): Laws instance
%   methodName (charstring): Method name describing the law
%   sets (cell{1,:} of `+replab.Samplable`): Samplable sets to use for the law parameters
    nArgs = length(sets);
    for i = 1:nRuns
        skipRun = false;
        try
            args = cell(1, nArgs);
            for j = 1:nArgs
                s = sets{j};
                args{j} = s.sample;
            end
            % Run the test
            laws.(methodName)(args{:});
        catch
            err = lasterror;
            switch err.identifier
              case 'replab:inexistent'
                fprintf('?');
                return
              case 'replab:skip'
                fprintf('skipping slow test');
                return
              otherwise
                    rethrow(err);
            end
        end
    end
end
