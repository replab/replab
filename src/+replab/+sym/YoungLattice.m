classdef YoungLattice < replab.Obj
% Young Lattice
%
% Representation of a subgraph of the Young's lattices, which connects Young Diagrams that differ by the
% removal/addition of one box.
% The subgraph is the set of all diagrams contained in a particular diagram.
% The subgraph is multipartite, with different levels' for each k from 1 to n.

    properties (SetAccess = protected)
        partition  % The partition at the top of the graph
        domainSize % Domain size of the partition at the top of the graph
        sets       % (cell(1,n) of `+replab.+perm.Set`): The k'th element holds data for the partitions of level k
        above      % (cell(1,n-1) of integer(\*,\*)): The k'th element is a matrix recording the connections between the partitions level of of k and k+1.
        below      % (cell(1,n-1) of integer(\*,\*)): The k'th element is a matrix recording the connections between the partitions level of of k+1 and k.

        % These are sort of like an adjacency matrix but the rows and columns are two different levels of the graph.
        % The connections are 'coloured' by an integer labeling the box that is removed/added
        % to make the connection. (The boxes are labeled left to right, then
        % top to bottom)

        labels     % (integer(\*,\*)): Matrix representing the labeling of boxes

        %        1 2 3 4
        % E.g.:  5 6 7 0 represents the labeling of 8 = 4+3+1
        %        8 0 0 0

        % nTabs % (cell of integer(1,\*)): The k'th element lists the number of standard tableaux for diagrams in level k.
        %
        % The rows of the following arrays describe standard tableaux which
        % can be represented as unique paths up this subgraph

        standRows % (integer(\*,n)): Lists the rows of the Young Diagram index k is in
        standCols % (integer(\*,n)): Lists the columns of the Young Diagram index k is in
        standTabs % (integer(\*,n)): Lists the label of the Young Diagram index k is in

        %standEigs % integer(:,n): Lists the eigenvalues of class functions associated with partitions found on the tableau's path
        %tabFuncs % integer(:,n): Lists tableaux function (the square of the change of basis coefficent from seminormal to orthogonal for each tableaux
    end

    methods(Static)

        function partFunc = nPartitions(n)
            partFunc = [1 zeros(1,n)];
            nums = 1:n;
            maxSum1 = floor((sqrt(24*(nums)+1)+1)/6);
            maxSum2 = floor((sqrt(24*(nums)+1)-1)/6);
            for i = nums
                j1 = 1:maxSum1(i);
                j2 = 1:maxSum2(i);
                partFunc(i+1) = sum( (-1).^(j1-1) .* partFunc(i-j1.*(3.*j1-1)/2+1) ) + sum( (-1).^(j2-1) .* partFunc(i-j2.*(3.*j2+1)/2+1) );
            end
            partFunc = partFunc(2:n+1);
        end

        function labels = findLabels(part)
            n = numel(part);
            inds = find(part);
            labels = zeros(sum(part),max(inds));
            rowCount = sum(part);
            labelCount = 0;
            for j = inds
                for k = 1:part(j)
                    labels(rowCount,1:j)= ((n-j+1):n)+labelCount;
                    rowCount = rowCount-1;
                    labelCount = labelCount - j;
                end
            end
        end

    end

    methods



        function self = YoungLattice(part,n)
            self.domainSize = n;
            if sum(part) == n &&(part(1) ~= n || numel(part) == 1)
                part = arrayfun(@(k) nnz(part == k),1:n);
            end
            self.partition = part;
            self.labels = replab.sym.YoungLattice.findLabels(part);
            self.createLattice(part);
        end

        function self = createLattice(self,part)
            n = self.domainSize;
            self.sets = cell(1,n);
            for m = 1:n
                self.sets{m} = replab.perm.Set(n);
            end
            self.sets{n}.insert(part');
            self.below = cell(1,n-1);
            self.above = cell(1,n-1);
            self.below(1:(n-1)) = {zeros(n,3)};
            count = 1;
            subPartitions(part,n);
            for k = 1:(n-1)
                in = @(n) nonzeros(self.below{k}(:,n));
                self.below{k} = sparse(in(1),in(2),in(3));
                self.above{k} = self.below{k}';
            end
            function subPartitions(sup,k)
                if k ==1
                    return
                end
                for i = find(sup)
                    sup2 = sup;
                    switch i
                        case 1
                            sup2(i) = sup2(i)-1;
                        otherwise
                            sup2([i-1 i]) =  sup2([i-1 i])+[1 -1];
                    end
                    row = self.sets{k-1}.find(sup2');
                    col = self.sets{k}.find(sup');
                    val = self.labels(sum(sup(i:n)),i);
                    if ~row
                        row = self.sets{k-1}.insert(sup2');
                        self.below{k-1}(count,:) = [row col val];
                        count = count + 1;
                        subPartitions(sup2,k-1);
                    else
                        self.below{k-1}(count,:) = [row col val];
                        count = count+1;
                    end
                end
                for j = 1:(n-1)
                end
            end
        end

        function mat = symAndAntiSym(self)
            if self.partition(end)
                mat = [1 0]';
                return
            elseif self.partition(1) == numel(self.partition)
                mat = [0 1]';
                return
            end
            n = self.domainSize;
            mat = full(logical(self.below{2}));
            for k = 4:n
                mat = mat*double(logical(self.below{k-1}));
            end
        end

        function nTabs = numTableaux(self)
            n = self.domainSize;
            mat = 1;
            nTabs = cell(1,n);
            nTabs{1} = mat;
            for k = 2:n
                mat = mat*logical(self.below{k-1});
                nTabs{k} = full(mat);
            end
        end

        function [standRows, standCols, standTabs] = generateTableaux(self)
            n = self.domainSize;
            part = self.partition;
            d = self.numTableaux{n};
            words = replab.sym.Words(flip(repelem(1:n,part)));
            standRows = zeros(d,n);
            standCols = zeros(d,n);
            standTabs = zeros(d,n);
            tabVec =zeros(1,n-1);
            tabVec(1) = 1;
            count = 1;
            rec(1,1);
            self.standRows = standRows;
            self.standCols = standCols;
            self.standTabs = standTabs;
            function rec(i,column)
                if i == n
                    standTabs(count,:) = tabVec;
                    standRows(count,:) = words.word(tabVec);
                    standCols(count,:) = words.conjWord(tabVec);
                    count = count + 1;
                    return
                end
                for j = find(self.above{i}(:,column))'
                    self.above{i}(j,column);
                    tabVec(i+1) = self.above{i}(j,column);
                    rec(i+1,j);
                end
            end
        end

        function standEigs = generateEigenVals(self)
            rows = self.standRows;
            cols = self.standCols;
            standEigs = rows(:,2:end)-cols(:,2:end);
        end

        function psi = generateTabFun(self)
            n = self.domainSize;
            part = self.partition;
            words = replab.sym.words(flip(repelem(1:n,part)));
            d = self.numTableaux{n};
            psi = zeros(1,d);
            count = 1;
            rec(1,1,[1 zeros(1,sum(part)-1)],zeros(1,n-1));
            function rec(i,column,maxInRows,phiVec)
                if i == n
                    psi(count) = prod(phiVec);
                    count = count + 1;
                    return
                end
                for j = find(self.above{i}(:,column))'
                    maxInRows2 = maxInRows;
                    phiVec2 = phiVec;
                    index = self.above{i}(j,column);
                    row = words.word(index);
                    col = words.conjWord(index);
                    phiVec2(i) =phiFun(maxInRows2,row,col);
                    maxInRows2(row) = index;
                    rec(i+1,j,maxInRows2,phiVec2);
                end
            end

                function phi = phiFun(maxInRow,iRow,iCol)
                    phi = 1;
                    if iRow == 1
                        return
                    end
                    for k = 1:numel(maxInRow)
                        if maxInRow(k)
                            if k == iRow
                                return
                            end
                            axDist = (abs(words.word(maxInRow(k))-iRow)+...
                                abs(words.conjWord(maxInRow(k))-iCol));
                            phi = phi*(1+1/axDist);
                        end
                    end
                end
            end
    end
end
