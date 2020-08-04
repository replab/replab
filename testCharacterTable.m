G = replab.S(3);
cr = {[1 2 3] [1 3 2] [2 3 1]};
cc = cellfun(@(r) G.conjugacyClass(r), cr, 'uniform', 0);
cn = {'e' '(1 2)' '(1 2 3)'};
rn = {'t' 's' 'S'};
ce = {'1' '1' '1'
      '1' '-1' '1'
      '2' '0' '-1'};
ct = replab.CharacterTable.make(replab.S(3), cc, cn, ce, rn, {[] [] []});
