function X = trace(self)
% trace operator
    X = 0;
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
        X = X + trace(self.blocks{i})*dim;
    end
end
