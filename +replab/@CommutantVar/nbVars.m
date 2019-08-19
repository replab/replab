function n = nbVars(self)
% Returns the number of SDP variable used be the object
    n = 0;
    for i = 1:self.nComponents
        n = n + length(getvariables(self.blocks{i}));
    end
end
