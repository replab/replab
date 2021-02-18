function subs = absoluteSplitInParent_real_unitary(sub, sample)
% Decomposes a real unitary subrepresentation into subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Unitary subrepresentation to split further
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Subrepresentations with their ``.parent`` set to the ``.parent`` of ``sub``
    assert(sub.isUnitary);
    assert(sub.overR);
    tol = replab.globals.doubleEigTol;
    d = sub.dimension;
    % force a dense matrix to have eig behave well
    subI = sub.injection('double/sparse');
    subP = sub.projection('double/sparse');
    S = full(subP * sample * subI);
    % make a symmetric matrix so that the eigenvectors are orthogonal and the decomposition is real
    X = (S + S')/2;
    [U, D] = eig(X);
    D = reshape(diag(D), [1 d]);
    P = replab.Partition.fromApproximateVector(D, tol);
    blocks = P.blocks;
    n = length(blocks);
    subs = cell(1, n);
    for i = 1:n
        divisionAlgebraName = '';
        blk = blocks{i};
        dblk = length(blk);
        if mod(dblk, 2) == 1
            % if the dimension is odd, the irrep has to be of real-type
            basis = U(:, blk);
            % assert(norm(Y, 'fro') <= tol); % TODO: put a real epsilon
        else
            % otherwise we do not know
            % now sample from the skew symmetric space to see if the representation is of
            % complex-type or quaternion-type
            Y = U(:, blk)'*S*U(:, blk);
            Y = (Y - Y')/2;
            % real eigenvalue decomposition of a skew-symmetric matrix,
            % we use the real Schur that returns an orthogonal basis
            % this could be more efficient with a real skew-symmetric EV algo
            [U1, T1] = schur(Y);
            % enforce a standard form on the imaginary EV
            for j = 1:2:dblk
                if T1(j+1,j) < T1(j,j+1)
                    % enforces the [0 -a; a 0] block structure with a > 0
                    U1(:,[j+1 j]) = U1(:,[j j+1]);
                    T1(:,[j+1 j]) = T1(:,[j j+1]);
                    T1([j+1 j],:) = T1([j j+1],:);
                end
            end
            % enforce the form we are looking for
            diagm = diag(T1, 1);
            diagp = diag(T1, -1);
            lambda = sum(diagp(1:2:end) - diagm(1:2:end))/dblk;
            diag1 = zeros(1, dblk-1);
            diag1(1:2:end) = lambda;
            % this is the 2x2-block-diagonal matrix after clean up
            D1 = diag(diag1, -1) - diag(diag1, 1);
            % and the associated error
            E1 = T1 - D1;
            % check that the error is reasonable, TODO: make this robust
            assert(norm(E1, 'fro') <= tol);
            if norm(D1, 'fro') > tol
                % there is content in the diagonal component
                basis = U(:, blk) * U1;
                divisionAlgebraName = 'complex';
            else
                % the diagonal component is noise
                basis = U(:, blk);
            end
        end

        I = subI * basis;
        if sub.mapsAreAdjoint
            % as subP = subI'
            P = I';
        else
            P = basis' * subP;
        end
        subs{i} = sub.parent.subRep(I, 'projection', P, 'isUnitary', true, 'divisionAlgebraName', divisionAlgebraName);
    end
end
