function C = permGrpClosure(G, el)
    if G.contains(el)
        % trivial closure
        C = G;
        return
    end
    n = G.n;
    C = replab.PermGrp([G.generators; el]);
    % G is the current group
    % C is the new group
    Gelements = G.hashedSortedElements;
    Gorder = Gelements.nElements;
    Gwords = G.words;
    elWord = replab.Word.generator(C.nGenerators);
    % if <G>^<el> = <G> then <C> = <G> * <el>
    % this follows closely grp.gi in Gap System
    if G.elementNormalizes(el)
        rep = el;
        repWord = elWord;
        shift = 0;
        newElements = zeros(G.n, 0, 'int32');
        newWords = {};
        while ~G.contains(rep)
            newElements = [newElements zeros(n, Gorder)];
            newWords = horzcat(newWords, cell(1, Gorder));
            for i = 1:Gorder
                % we cannot have duplicates here
                g = Gelements.M(:, i);
                newElements(:, shift+i) = replab.Perm.compose(g, rep);
                newWords{shift+i} = Gwords{i} * repWord;
            end
            shift = shift + Gorder;
            rep = replab.Perm.compose(rep, el);
            repWord = repWord * elWord;
        end
        Celements = Gelements.append(newElements);
        Cwords = horzcat(Gwords, newWords);
        C.setAllElements(Celements, Cwords);
    else
        % otherwise use a Dimino step
        reps = (1:n);
        repWords = {replab.Word.identity};
        Celements = Gelements;
        Cwords = Gwords;
        while size(reps, 2) > 0
            rep = reps(:, 1);
            repWord = repWords{1};
            reps = reps(:, 2:end);
            repWords = repWords(2:end);
            for i = 1:C.nGenerators
                gen = C.generator(i);
                genWord = replab.Word.generator(i);
                rg = replab.Perm.compose(rep, gen);
                rgWord = repWord * genWord;
                if Celements.find(rg) == 0
                    mem = ismember(rg(:)', Celements.M', 'rows');
                    assert(sum(mem) == 0);
                    newElements = zeros(n, Gorder, 'int32');
                    newWords = cell(1, Gorder);
                    for i = 1:Gorder
                        e = Gelements.M(:,i);
                        eWord = Gwords{i};
                        newElements(:,i) = replab.Perm.compose(e, rg);
                        newWords{i} = eWord * rgWord;
                    end
                    reps = [reps rg];
                    repWords = horzcat(repWords, {rgWord});
                    Celements = Celements.append(newElements);
                    Cwords = horzcat(Cwords, newWords);
                end
            end
        end
        C.setAllElements(Celements, Cwords);
    end

end
