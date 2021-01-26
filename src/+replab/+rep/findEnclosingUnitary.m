function I = findEnclosingUnitary(X)
% Given an approximate orthogonal or unitary matrix, returns an interval matrix containing a close orthogonal/unitary matrix
%
% We use a first-order approximation of the Newton iteration corresponding to the equation
% `` X' * X - eye(d) == 0 ``, see `<https://en.wikipedia.org/wiki/Orthogonal_matrix#Nearest_orthogonal_matrix>_`
% by computing the Newton iteration in interval arithmetic. Then either:
%
% - the interval result of the iteration does not contain the previous iterate, thus the precision can still be increased and we repeat the Newton iteration,
% - the interval result of the iteration does contain the previous iterate, in which case we assume it would contain the exact result as well.
%
% TODO: prove that this is a formal result
%
% Args:
%   X (double(d,d)): Approximately orthogonal or unitary matrix with condition number less than 3
%
% Returns:
%   intval(d,d): Interval matrix containing a nearby orthogonal or unitary matrix
    X0 = X;
    d = size(X, 1);
    while 1
        Xi = intval(X);
        N = Xi'*Xi;
        P = Xi*N/2;
        I = 2*Xi + P*(N - 3*eye(d)); % result of the Newton iteration
        if any(any(isnan(I)))
            if cond(X0) > 2.9
                error('The condition number of the input matrix should be less than 3.');
            else
                error('Unknown error');
            end
        end
        if in(X, I)
            return
        end
        X = mid(I);
    end
end
