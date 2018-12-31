function runNTimes(nRuns, obj, methodName, sampleFuns)
    nArgs = length(sampleFuns);
    for i = 1:nRuns
        errored = false;
        args = cell(1, nArgs);
        for j = 1:nArgs
            sampleFun = sampleFuns{j};
            try
                % We catch "inexistent" errors during generation
                % and skip the test in that case
                args{j} = sampleFun();
            catch
                err = lasterror;
                if ~isequal(err.identifier, 'replab:inexistent')
                    rethrow(err);
                else
                    errored = true;
                    fprintf('?');
                end
            end
        end
        if ~errored
            obj.(methodName)(args{:});
        end
    end
end
