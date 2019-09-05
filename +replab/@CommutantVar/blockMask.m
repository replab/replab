function M = blockMask(self)
% Returns a mask showing the block structure
    M = sparse(self.dim, self.dim);
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
        M(co+[1:dim*size(self.blocks{i},1)], co+[1:dim*size(self.blocks{i},1)]) = kron(ones(size(self.blocks{i})), eye(dim));
        co = co + dim*size(self.blocks{i},1);
    end
end
