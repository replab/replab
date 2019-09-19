function S = sparse_(varargin)
% Wraps the `sparse` Octave/Matlab function such that small matrices are returned dense
    S = sparse(varargin{:});
    d = replab.Settings.sparseCutoff;
    if prod(size(S)) < d*d
        S = full(S);
    end
end
