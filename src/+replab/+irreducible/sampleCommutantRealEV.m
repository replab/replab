function [X, err] = sampleCommutantRealEV(rep)
% Samples a matrix that commutes with the given representation and has real eigenvalues
%
% Args:
%   rep (`+replab.Rep`): Representation
%
% Returns
% -------
%   X: double(\*,\*)
%     Commutant sample
%   err: double
%     Estimated error on the commutant averaging
    dom = replab.domain.SelfAdjointMatrices(rep.field, rep.dimension);
    X = dom.sample;
    if ~rep.knownUnitary
        U = rep.unitarize;
        X = U.Ainv * X * U.A;
    end
    [X, err] = rep.commutant.project(X);
end
