function runNTimes(nRuns, obj, methodName, sampleFuns)
    nArgs = length(sampleFuns);
    for i = 1:nRuns
        args = cell(1, nArgs);
        for j = 1:nArgs
            sampleFun = sampleFuns{j};
            args{j} = sampleFun();
        end
        obj.(methodName)(args{:});
    end
end
