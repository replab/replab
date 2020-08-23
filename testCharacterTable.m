G = replab.S(3); % group
cr = {[1 2 3] [1 3 2] [2 3 1]}; % conjugacy class representatives (not necessarily lex-minimal)
cc = cellfun(@(r) G.conjugacyClass(r), cr, 'uniform', 0); % compute conjugacy classes
                                                          % this will compute the lex-minimal
                                                          % rep, and the centralizer
cn = {'e' '(1 2)' '(1 2 3)'}; % labels of conjugacy classes (optional)
rn = {'t' 's' 'S'}; % labels of irreps
ce = {'1' '1' '1'   % table of character values
      '1' '-1' '1'
      '2' '0' '-1'};
irreps = {[] [] []}; % concrete realizations of the irreps
ct = replab.CharacterTable.make(replab.S(3), cc, cn, ce, rn, {[] [] []});
