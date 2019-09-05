function F = gt(X,Y)
    % greater than constraint
    warning('Strict inequalities cannot be guaranteed. Imposing a non-strict constraint instead.');
    F = (X >= Y);
end
