classdef Symmetric

    methods (Static)

        %        function ct = realCharacterTable(G)

% $$$         function ct = realCharacterTable(G)
% $$$             n = G.permutationGroup.domainSize;
% $$$             allData = replab.sym.CTData.instance(n);
% $$$             nData = allData(n);
% $$$             irrepNames = cellfun(@(part) strrep(replab.shortStr(replab.sym.findPartitions.conjugatePart(part)), ' ', ''),nData.partitions.list,'UniformOutput', false);
% $$$             [classes, partCell] = replab.sym.findConjClasses(permGroup);
% $$$             classNames = cellfun(@(c) strrep(replab.shortStr(c), ' ', ''), partCell, 'uniform', 0);
% $$$             chars = replab.cyclotomic(permCTValues(n));
% $$$             permCT = replab.ComplexCharacterTable(G.permutationGroup, replab.ConjugacyClasses(permGroup, classes), chars, 'irrepNames', irrepNames, 'classNames', classNames);
% $$$             ct = permCT.imap(G.niceMorphism);
% $$$
% $$$             function CT = permCTValues(n)
% $$$                 nParts = nData.nParts;
% $$$                 CTunNormed = zeros(nParts);
% $$$                 for row = 1:nParts
% $$$                     partition = nData.partitionList{row};
% $$$                     lin = replab.sym.SymPoly.powerToMonomial(partition,n);
% $$$                     CTunNormed(:,row) = lin.coeffs';
% $$$                 end
% $$$                 CT = grahamSchmidt(CTunNormed,@(x,y) round(sum(x.*y.*nData.conjSizes/nData.fact)));
% $$$                 CT = flip(round(CT));
% $$$                 function mat = grahamSchmidt(mat,innerProd)
% $$$                     d = nData.nParts;
% $$$                     for i = d:-1:2
% $$$                         for j = 1:(i-1)
% $$$                             mat(j,:) = mat(j,:)-innerProd(mat(j,:),mat(i,:))*mat(i,:);
% $$$                         end
% $$$                     end
% $$$                 end
% $$$             end
% $$$
% $$$         end

        function G = make(n)
        % Constructs the symmetric group of degree n, as an abstract group with character tables
            assert(n > 2);
            name = sprintf('Symmetric group S(%d) of degree %d', n, n);
            % Permutation realization
            S = [2:n 1];
            T = [2 1 3:n];
            % this is the presentation from page 2100 of
            % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
            relators = {['s^' num2str(n)], 't^2', ['(s*t)^' num2str(n-1)]};
            for j = 2:floor(n/2)
                relators{1,end+1} = sprintf('(t^-1 s^-%d t s^%d)^2', j, j);
            end
            G = replab.AbstractGroup({'s' 't'}, relators, 'permutationGenerators', {S, T}, 'order', replab.util.factorial(n), 'name', sprintf('Symmetric group S(%d)', n), 'inAtlas', true);
            partitions = replab.sym.IntegerPartition.all(n);
            classes = cellfun(@(p) p.conjugacyClass, partitions, 'uniform', 0);
            classes = replab.ConjugacyClasses.sorted(G.permutationGroup, classes);
            irreps = cellfun(@(p) replab.sym.SymmetricSpechtIrrep(G.permutationGroup, p.partition).rep, partitions, 'uniform', 0);
            values = zeros(length(irreps), classes.nClasses);
            for i = 1:length(irreps)
                for j = 1:classes.nClasses
                    values(i, j) = trace(irreps{i}.image(classes.classes{j}.representative));
                end
            end
            values = replab.cyclotomic(values);
            ctR = replab.RealCharacterTable(G.permutationGroup, classes, values, 'irreps', irreps);
            ctC = replab.ComplexCharacterTable.fromRealCharacterTable(ctR);
            G.cache('conjugacyClasses', classes.imap(G.niceMorphism.inverse), 'error');
            G.cache('realCharacterTable', ctR.imap(G.niceMorphism.inverse), 'error');
            G.cache('complexCharacterTable', ctC.imap(G.niceMorphism.inverse), 'error');
        end

        function R = recognize(G)
        % Recognizes if the given group is the symmetric group and provides the generators according to the standard presentation
            R = [];
            [n, r] = replab.util.unfactorial(G.order);
            if r ~= 0
                return
            end
            n = double(n);
            C = G.conjugacyClasses.classes;
            entry = replab.atl.Symmetric.make(n);
            for i = 1:length(C)
                S = C{i};
                s = S.representative;
                if G.elementOrder(s) == n
                    for j = 1:length(C)
                        T = C{j};
                        if G.elementOrder(T.representative) == 2
                            U = T.elements;
                            for k = 1:length(U)
                                t = U{k};
                                if entry.isMorphismByImages(G, 'images', {s t})
                                    if G.subgroup({s, t}).order == G.order
                                        R = entry.isomorphismByImages(G, 'images', {s t});
                                        return
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

    end

end
