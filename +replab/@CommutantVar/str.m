function s = str(self)
% Nice string representation
    s = ['SDP matrix of size ', num2str(self.dim), 'x', num2str(self.dim), ' with ', num2str(self.nbVars), ' variables.'];
    s = [s, char(10)];
    s = [s, 'Block structure: '];
    for i = 1:self.nComponents
        switch self.types(i)
            case 'R'
                s = [s, num2str(self.dimensions1(i)), '*', num2str(self.multiplicities(i)), 'x', num2str(self.multiplicities(i)), ' + '];
            case 'C'
                s = [s, num2str(self.dimensions1(i)/2), '*', num2str(2*self.multiplicities(i)), 'x', num2str(2*self.multiplicities(i)), ' + '];
            case 'H'
                s = [s, num2str(self.dimensions1(i)/4), '*', num2str(4*self.multiplicities(i)), 'x', num2str(4*self.multiplicities(i)), ' + '];
            otherwise
                error('Unknown type');
        end
    end
    s = s(1:end-3);
end
