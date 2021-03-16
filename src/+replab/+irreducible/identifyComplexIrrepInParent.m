function sub1 = identifyComplexIrrepInParent(sub, sample)
% Identifies the type of a complex irreducible representation and shapes its division algebra encoding
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation over C to identify
%   sample (double(\*,\*)): Sample of ``sub.parent.antilinearInvariant``
%
% Returns:
%   `+replab.SubRep`: Subrepresentation with the Frobenius-Schur indicator set, and the division algebra in the canonical form
    d = sub.dimension;
    J = sub.projection * sample * conj(sub.injection);
    lambda = real(trace(conj(J)*J)/d);
    tol = 1e-10;
    if abs(lambda) < tol
        sub.cache('frobeniusSchurIndicator', 0, '==');
        sub1 = sub;
    elseif lambda > 0
        J = J / sqrt(lambda);
        % now we have conj(J)*J = eye(d) represents an antilinear involution
        % rho(g)*J = J*conj(rho(g))
        % we look for A such that inv(A)*J*conj(A) = eye(d), as then
        % the similar representation A*rho(g)*inv(A) would be real:
        % A*rho(g)*inv(A)*J = J*conj(A)*conj(rho(g))*conj(inv(A))
        % rho(g) * inv(A)*J*conj(A) = inv(A)*J*conj(A) * conj(rho(g))
        if sub.isUnitary
            % we have rho(g)*J = J*conj(rho(g))
            % taking the transpose on both sides
            % J.'*rho(g).' = rho(g)'*J.'
            % J.'*conj(rho(ginv)) = rho(ginv)*J.' as rho(ginv) = rho(g)'
            % substituting ginv -> g
            % rho(g)*J.' = J.'*conj(rho(g))
            %
            % thus J and J.' represent both an equivariant antilinear map
            % as such maps are proportional, we have J = alpha * J.', which
            % means either alpha = 0 or J = J.'
            %
            % now, as J = J.' and as J*conj(J) = 1, we have J*conj(J.') = J*J' = eye(d)
            % and thus J is unitary
            % We have then have a Takagi decomposition J = U.'*D*U with D = eye(d)
            J = (J + J.')/2; % make sure J is symmetric
            [U,~] = replab.numerical.takagi(J);
            A = U.';
            % now inv(A)*U.'*U*conj(A) = (inv(U.')*U.')*U*conj(U.') = eye(d)*U*U' = eye(d)
            % as U' = inv(U) by construction of the Takagi decomposition
        else
        end

        sub1 = sub.withUpdatedMaps(sub.injection('double/sparse') * A, A \ sub.projection('double/sparse'), 'frobeniusSchurIndicator', 1, 'divisionAlgebraName', 'R->C');
    else % lambda < 0
        error('Unsupported');
    end
end
