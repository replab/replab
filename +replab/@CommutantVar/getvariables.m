function vars = getvariables(self)
% Returns the SDP variable used be the object
    vars = getvariables(self.blocks{1});
    for i = 2:self.nComponents
        vars = [vars, getvariables(self.blocks{i})];
    end
end
