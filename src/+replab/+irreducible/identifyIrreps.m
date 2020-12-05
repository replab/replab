function isIrrep = identifyIrreps(rep, subreps, failureProb, safetyFactor1, safetyFactor2)
% Among subrepresentations of a common representation, identify which are irreducible
%
% Args:
%   rep (`+replab.Rep`): Parent representation
%   subreps (cell(1,\*) of `+replab.SubRep`): Subrepresentations of ``rep``
%   safetyFactor1 (double, optional): Factor on the error estimation, below which it is assumed to cluster, default 100
%   safetyFactor2 (double, optional): Factor on the error estimation, above which it is assumed to split, default 1000
%   failureProb (double, optional): For each subrepresentation, upper bound on ``P(res = irrep|subrep is not irrep)``, default 1e-9
%
% Returns:
%   logical(1,\*): Whether the representation has been identified as irreducible
    if nargin < 5 || isempty(safetyFactor2)
        safetyFactor2 = 1000;
    end
    if nargin < 4 || isempty(safetyFactor1)
        safetyFactor1 = 100;
    end
    if nargin < 3 || isempty(failureProb)
        failureProb = 1e-9;
    end
    d = rep.dimension;
    X = randn(d, d);
    if rep.overR
        R = randn(d, d);
        I = randn(d, d);
        [R1, Rerr] = rep.commutant.project(R);
        [I1, Ierr] = rep.commutant.project(I);
        projErr = Rerr + Ierr;
        X = R + 1i * I;
        X1 = R1 + 1i * I1;
    else
        X = randn(d, d) + 1i * randn(d, d);
        [X1, projErr] = rep.commutant.project(X);
    end
    diff = X - X1;
    n = length(subreps);
    isIrrep = false(1, n);
    for i = 1:n
        sub = subreps{i};
        di = sub.dimension;
        I = sub.injection;
        P = sub.projection;
        A = P*X1*I;
        error = abs(trace(P*diff*I));
        sameUB = error * safetyFactor1 + projErr;
        splitLB = error * safetyFactor2 + projErr;
        mu = trace(A)/di;
        Delta = norm(A - mu * eye(di), 'fro');
        fpLB = sqrt(-2*log1p(-failureProb));
        if fpLB < splitLB
            % TODO: iterate the identification procedure to use a higher probability tolerance, but perform several tests
            error('Not enough precision to satisfy the failure probability requirements.');
        end
        if Delta > fpLB
            isIrrep(i) = false;
        elseif Delta > sameUB
            error('Error estimation failure');
        else
            isIrrep(i) = true;
        end
    end
end
