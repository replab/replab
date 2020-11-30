function X = bestStorage(X)
    storageCutoff = 1e6;
    densityCutoff = 0.1;
    if prod(size(X)) > storageCutoff
        % we try to minimize memory use there
        if isreal(X)
            denseMemSize = prod(size(X))*8;
            sparseMemSize = 16*nnz(X) + 8*size(X, 2) + 8; % assumes 64 bits
        else
            denseMemSize = prod(size(X))*16;
            sparseMemSize = 24*nnz(X) + 8*size(X, 2) + 8; % assumes 64 bits
        end
        if denseMemSize > sparseMemSize
            X = sparse(X);
        else
            X = full(X);
        end
    else
        % we try to minimize CPU time there
        if nnz(X) / prod(size(X)) < densityCutoff
            X = sparse(X);
        else
            X = full(X);
        end
    end
end
