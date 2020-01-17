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
        args = cell(1, nArgs);
        for j = 1:nArgs
            s = sets{j};
            try
                % We catch "inexistent" errors during generation
                % and skip the test in that case
                args{j} = s.sample;
            catch
                err = lasterror;
                if ~isequal(err.identifier, 'replab:inexistent')
                    rethrow(err);
                else
                    skipRun = true;
                    fprintf('?');
                end
            end
        end
        if ~skipRun
            % Run the test
            laws.(methodName)(args{:});
        end
    end
end
