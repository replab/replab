function [M, err] = projectAndFactorFromParent_sdpvar(E, value)
    v = getbasematrix(value, 0);
    [M, err] = E.projectAndFactorFromParent(v);
    ind = getvariables(value);
    for k = 1:length(ind)
        v = getbasematrix(value, ind(k));
        [Mk, errk] = E.projectAndFactorFromParent(v);
        M = M + recover(ind(k)) * Mk;
        err(1,end+1) = errk;
    end
end
