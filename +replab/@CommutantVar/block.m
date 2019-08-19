function M = block(self, i)
% Returns the desired block
    if (i < 1) || (i > self.nComponents)
        error('Block number out of bound');
    end

    M = self.blocks{i};
end
