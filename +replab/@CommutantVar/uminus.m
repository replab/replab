function X = uminus(self)
    % unary minus operator
    
    X = self.copy;
    for i = 1:X.nComponents
        X.blocks{i} = -X.blocks{i};
    end
end
