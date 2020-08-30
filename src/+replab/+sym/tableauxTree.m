function [connectList,transList,offDiagList,diagList] = tableauxTree(part)
    n = sum(part);
    [rows,cols,~] = replab.sym.YoungLattice(part,n).generateTableaux;
    hash = replab.perm.Set(n);
    hash.insert(rows');
    len = size(rows,1);
    swaps = zeros(n-1,n);
    for j = 1:(n-1)
        swaps(j,:) = swap(j,j+1);
    end
    connections = ~(rows(:,1:(n-1))==rows(:,2:n)|cols(:,1:(n-1))==cols(:,2:n));
    connectList = zeros(1,len-1);
    transList = zeros(1,len-1); 
    diagList = zeros(1,len-1);
    offDiagList = zeros(1,len-1); 
    for i = 1:len
        for k = find(connections(i,:))
            ind = hash.find(rows(i,swaps(k,:))');
            if ind ~=1
                if ~connectList(ind-1) 
                    connectList(ind-1) = i;
                    transList(ind-1) = k;
                    axDist = abs(rows(i,k)-rows(i,k+1))+abs(cols(i,k)-cols(i,k+1));
                    diagList(ind-1) = sign(rows(i,k)-rows(i,k+1))/axDist;
                    offDiagList(ind-1) = sqrt(1-1/axDist^2);
                end
            end
        end
    end
    function perm = swap(a,b)
            perm = 1:n;
            [perm(a),perm(b)] = deal(perm(b),perm(a));
    end
end