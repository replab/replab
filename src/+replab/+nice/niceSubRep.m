function subRep1 = niceSubRep(subRep)
    subRep1 = [];
    if ~isempty(subRep.divisionAlgebraName)
        return
    end
    P = subRep.projector;
    P = replab.numerical.refineProjector(P);
    P = replab.numerical.refineProjector(P);
    % try to make the subrepresentation injection/projection maps real
    tol = replab.globals.doubleEigTol;
    maximumDenominator = 1200;
    if ~replab.numerical.isNonZeroMatrix(imag(P), replab.globals.doubleEigTol)
        P = full(real(P));
        C = replab.numerical.recoverRational(P, tol, maximumDenominator);
        Cdouble = double(C);
        rationalBasis = false;
        if ~isempty(C)
            if subRep.parent.isExact
                C1 = subRep.parent.commutant.project(C, 'exact');
                if all(all(C == C1))
                    rationalBasis = true;
                end
            else
                [Pproj, error1] = subRep.parent.commutant.project(P);
                [Cproj, error2] = subRep.parent.commutant.project(Cdouble);
                if norm(Cproj - Cdouble, 'fro') <= norm(Pproj - P, 'fro') + error1 + error2
                    rationalBasis = true;
                end
            end
            if rationalBasis
                [~,~,jb] = qr(Cdouble, 'vector');
                injection = C(:, jb(1:subRep.dimension));
                subRep1 = subRep.parent.subRep(injection);
                return
            else
                [~,~,jb] = qr(P, 'vector');
                injection = P(:, jb(1:subRep.dimension));
                subREp1 = subRep.parent.subRep(injection, 'forceReal', true);
                return
            end
        end
    end
end
