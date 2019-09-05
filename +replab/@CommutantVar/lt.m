function F = lt(X,Y)
    % lesser than constraint
    warning('Strict inequalities cannot be guaranteed. Imposing a non-strict constraint instead.');
    F = (X <= Y);
end
