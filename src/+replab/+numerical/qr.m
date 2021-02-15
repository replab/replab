function [Q, R] = qr(A)
% Economy-size QR decomposition
%
% Provides similar results to ``[Q, R] = qr(A, 0)``, but corrects for diagonal signs
    if isa(A, 'replab.H')
        if nnz(A.Y) == 0
            [Q, R] = replab.numerical.qr(full(A.X));
            Q = replab.H(Q);
            R = replab.H(R);
        else
            Z = full(replab.H.encode(A));
            [Q, R] = replab.numerical.qr(Z);
            Q = replab.H.decode(Q);
            R = replab.H.decode(R);
            for i = 1:min(size(R))
                R.X(i,i) = real(R.X(i,i));
                R.Y(i,i) = 0;
            end
        end
    else
        [Q, R] = qr(full(A), 0);
        D = diag(sign(diag(R)));
        Q = Q*D;
        R = D*R;
    end
end
