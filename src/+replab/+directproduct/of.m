function prd = of(factors)
    className = replab.directproduct.commonClass(factors);
    switch className
      case 'replab.NiceFiniteGroup'
        prd = replab.directproduct.OfNiceFiniteGroups(factors);
      case 'replab.FiniteGroup'
        prd = replab.directproduct.OfFiniteGroups(factors);
      case 'replab.CompactGroup'
        prd = replab.directproduct.OfCompactGroups(factors);
    end
end
