function grp = DihedralGroup(self)
    n = self.domainSize;
    if n <= 2
        grp = self.symmetricGroup;
    else
        g1 = [2:n 1];
        g2 = fliplr(1:n);
        grp = self.subgroup({g1 g2});
    end
end
