classdef PermGrpList < replab.PermGrp
% Implementation of a permutation group that stores an enumeratio
% of all elements
    
    properties (Access = private)
        n; % domain size
        genMat; % m group generators in a int32 nGens x n matrix
                % one generator per row
    end
    
    methods
        
        function self = PermGrpList(genMat)
            self.n = size(genMat, 2);
            self.cat = replab.cat.PermAsGroup(self.n);
            self.genMat = genMat;
        end

        function d = domainSize(self)
            d = self.n;
        end
        
        function b = contains(self, permutation)
            b = self.hashedSortedElements.find(permutation) ~= 0;
        end
        
        function o = order(self)
            o = vpi(length(self.words));
        end
        
        function nG = nGenerators(self)
            nG = size(self.genMat, 1);
        end
        
        function p = generator(self, i)
            p = double(self.genMat(i, :));
        end

        function genMat = generatorsAsMatrix(self)
            genMat = double(self.genMat);
        end
        
        function p = uniformlyRandomElement(self)
            i = randi(double(self.order));
            p = self.element(i);
        end
        
        function p = element(self, i)
            i = double(i);
            p = double(self.hashedSortedElements.M(:, i)');
        end
        
        function E = elements(self)
            E = cell(1, double(self.order));
            for i = 1:double(self.order)
                E{i} = self.element(i);
            end
        end
        
        function E = elementsAsMatrix(self)
            E = double(self.hashedSortedElements.M');
        end
        
        function w = factorization(self, permutation)
            ind = self.hashedSortedElements.find(permutation);
            if ind == 0
                error(['Permutation ' num2str(permutation) ' is not contained in this group']);
            end
            genInd = self.wordRhs(ind);
            if genInd == 0
                w = replab.Word.identity;
            else
                rhs = self.generator(genInd);
                rhsW = replab.Word.generator(genInd);
                lhs = self.cat.composeWithInverse(permutation, rhs);
                w = self.factorization(lhs) * rhsW;
            end
            assert(isequal(self.evaluateWord(w), permutation));
        end
        
        %    w = self.words{ind};
        
        function p = evaluateWord(self, word)
            p = self.cat.evaluateWord(word, self.generators);
        end

    end
    
    methods (Static)
        
        function G = fromGeneratorsAsMatrix(genMat)
            domainSize = size(genMat, 2);
            assert(domainSize > 0);
            genMat = int32(genMat);
            nG = size(genMat, 1);
            for i = 1:nG
                p = genMat(i, :);
                errmsg = 'The identity cannot be part of the generators';
                assert(~replab.Perm.isIdentity(p), errmsg);
            end
            G = replab.prv.PermGrpList(genMat);
        end
        
    end
    
    properties
        hashedSortedElements_ = []; % list of elements
        wordRhs_ = []; % wordRhs_ has size 1 x nElements; wordRhs_(ind) encodes compactly
                       % the word representing el = self.element(ind)
                       %
                       % wordRhs_(ind) = 0  : el is the identity
                       % otherwise el = lhs * self.generator(wordRhs_(ind))
                       %
                       % and lhs is implicitly encoded by lhs = el * rhs.inverse
        words_ = {};   % list of short words corresponding to the above list
    end
    
    methods
        
        function disp(self)
            if ~isequal(self.hashedSortedElements_, [])
                disp(sprintf('Group of order %d with generators', double(self.order)));
            else
                disp('Group with generators');
            end
            disp(num2str(self.generatorsAsMatrix));
        end

        function b = isSubgroup(self, U)
        % Returns true if U is a subgroup of this group
            b = true;
            for i = 1:U.nGenerators
                if ~self.contains(U.generator(i))
                    b = false;
                    return
                end
            end
        end
        
        function b = isNormalizedByElement(self, h)
        % Returns true if this group G is normalized by the element h
            b = true;
            for i = 1:self.nGenerators
                c = self.cat.conjugate(h, self.generator(i));
                if ~self.contains(c)
                    b = false;
                    return
                end
            end
        end
        
        function b = isNormal(self, U)
        % Returns true if this group G normalizes the given group U
        % and false otherwise
        %
        % (i.e. for all g in G and u in U, g * u * gInv is in U)
        %
        % U does not need to be a subgroup of G
            b = true;
            for i = 1:self.nGenerators
                g = self.generator(i);
                for j = 1:U.nGenerators
                    u = U.generator(j);
                    c = self.cat.conjugate(u, g);
                    if ~U.contains(c)
                        b = false;
                        return
                    end
                end
            end
        end
        
    end
    
    methods% (Access = private)
    
    function H = hashedSortedElements(self)
        if isequal(self.hashedSortedElements_, [])
            self.computeAllElements;
        end
        H = self.hashedSortedElements_;
    end
    
    function W = words(self)
        if isequal(self.words_, {})
            self.computeAllElements;
        end
        W = self.words_;
    end
    
    function W = wordRhs(self, elementIndex)
        if isequal(self.wordRhs_, [])
            self.computeAllElements;
        end
        W = self.wordRhs_(elementIndex);
    end
    
    function H = subgroupRemovingLastGenerator(self)
    % Returns the subgroup generated by the first (nG-1) generators
    % of this group, with nG = the number of generators of this group
        assert(self.nGenerators > 0);
        newGenMat = self.genMat(1:end-1, :);
        H = replab.prv.PermGrpList(newGenMat);
    end
    
    function computeAllElementsNormalizing(self, H, genIndex)
    % Computes all elements of this group, given that
    %
    % - "genIndex" is the index of a generator of this group
    % - H does not contains "gen = self.generator(genIndex)"
    % - This group is the closure of H and gen
    % - H.isNormalizedByElement(gen) is true
    % where gen is the generator corresponding to genWord
        genWord = replab.Word.generator(genIndex);
        gen = self.generator(genIndex);
        Hhse = H.hashedSortedElements;
        Horder = double(Hhse.nElements);
        Hwords = H.words;
        H.computeAllElements;
        HwordRhs = H.wordRhs_;

        newElements = zeros(self.n, 0, 'int32');
        newWords = {};
        newWordRhs = [];
        shift = 0;
        rep = gen;
        repWord = genWord;
        while ~H.contains(rep)
            newElements = [newElements zeros(self.n, Horder, 'int32')];
            newWords = horzcat(newWords, cell(1, Horder));
            newWordRhs = horzcat(newWordRhs, genIndex * ones(1, Horder, 'int32'));
            for i = 1:Horder
                % we cannot have duplicates here
                h = Hhse.M(:, i);
                newEl = h(rep); % compose(h, rep)
                newElements(:, shift+i) = h(rep);
                newWords{shift+i} = Hwords{i} * repWord;
            end
            shift = shift + Horder;
            rep = rep(gen); % compose(rep, gen)
            repWord = repWord * genWord;
        end
        [Ghse, indices] = Hhse.append(newElements).lexColSorted;
        words = horzcat(Hwords, newWords);
        words = words(indices);
        wordRhs = horzcat(HwordRhs, newWordRhs);
        wordRhs = wordRhs(indices);
        self.hashedSortedElements_ = Ghse;
        self.words_ = words;
        self.wordRhs_ = wordRhs;
    end
    
    function computeAllElementsDiminoStep(self, H, genIndex)
    % Computes all elements of this group, given that
    %
    % - "genIndex" is the index of a generator of this group
    % - H does not contains "gen = self.generator(genIndex)"
    % - This group is the closure of H and gen
    % where gen is the generator corresponding to genWord
        
        genWord = replab.Word.generator(genIndex);
        gen = self.generator(genIndex);
        Hhse = H.hashedSortedElements;
        Horder = double(Hhse.nElements);
        Hwords = H.words;
        H.computeAllElements;
        HwordRhs = H.wordRhs_;

        % reps are the coset representatives having been already used
        % to generate parts of this group
        % of course, the elements in H are contained in self,
        % so we start the process with the identity
        reps = int32(self.cat.identity)';
        % with the associated words
        repWords = {replab.Word.identity};
        Ghse = Hhse; % G is self
        Gwords = Hwords;
        GwordRhs = HwordRhs;
        while size(reps, 2) > 0 % while there are representatives
            rep = reps(:, 1); % pop a representative
            repWord = repWords{1};
            reps = reps(:, 2:end);
            repWords = repWords(2:end);
            for i = 1:self.nGenerators
                gen = self.generator(i);
                genWord = replab.Word.generator(i);
                rg = self.cat.compose(rep, gen);
                rgWord = repWord * genWord;
                if Ghse.find(rg) == 0
                    % the new elements are written H.element * rep * gen
                    newElements = zeros(self.n, Horder, 'int32');
                    % index of the first new element
                    rgIndex = length(GwordRhs) + 1;
                    newWordRhs = i*ones(1, Horder, 'int32');
                    newWords = cell(1, Horder);
                    for j = 1:Horder
                        e = Hhse.M(:,j);
                        eWord = Hwords{j};
                        newElements(:, j) = e(rg); % compose(e, rg)
                        newWords{j} = eWord * rgWord;
                    end
                    reps = [reps rg];
                    repWords = horzcat(repWords, {rgWord});
                    Ghse = Ghse.append(newElements);
                    Gwords = horzcat(Gwords, newWords);
                    GwordRhs = horzcat(GwordRhs, newWordRhs);
                end
            end
        end
        [Ghse, indices] = Ghse.lexColSorted;
        self.hashedSortedElements_ = Ghse;
        self.words_ = Gwords(indices);
        self.wordRhs_ = GwordRhs(indices);
    end
    
    function computeAllElements(self)
        if self.isTrivial
            % trivial group
            M = int32(self.cat.identity)';
            self.hashedSortedElements_ = replab.prv.HashIntMatrix(M);
            self.words_ = {replab.Word.identity};
            self.wordRhs_ = int32([0]);
        else
            H = self.subgroupRemovingLastGenerator;
            k = self.nGenerators;
            if H.contains(self.generator(k))
                self.hashedSortedElements_ = H.hashedSortedElements;
                self.words_ = H.words;
            elseif H.isNormalizedByElement(self.generator(k))
                self.computeAllElementsNormalizing(H, k);
            else
                self.computeAllElementsDiminoStep(H, k);
            end
        end
    end
    
    end
end
