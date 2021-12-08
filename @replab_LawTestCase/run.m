function report = run(obj, report)
    if isempty(obj.domains)
        nRuns = 1;
    else
        nRuns = replab.Laws.nRuns;
    end
    fun = @() replab.laws.runNTimes(nRuns, obj.lawsInstance, obj.lawMethodName, obj.domains);
    delegated = MOxUnitFunctionHandleTestCase(getName(obj), getLocation(obj), fun);
    report = run(delegated, report);
end
