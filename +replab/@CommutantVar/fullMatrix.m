function M = fullMatrix(self)
% Constructs the full SDP matrix in the natural basis
    M = sdpvar(1);
    M(self.dim^2) = 0;
    M = reshape(M, self.dim*[1 1]);
    co = 0;
    for i = 1:self.nComponents
        dim = self.dimensions1(i);
        switch self.types(i)
            case 'R'
            case 'C'
                dim = dim/2;
            case 'H'
                dim = dim/4;
            otherwise
                error('Unknown type');
        end
        M(co + (1:dim*size(self.blocks{i},1)), co + (1:dim*size(self.blocks{i},1))) = kron(self.blocks{i}, eye(dim));
        co = co + dim*size(self.blocks{i},1);
    end
    M = self.U*M*self.U';
end
