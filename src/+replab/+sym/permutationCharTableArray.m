function CT = permutationCharTableArray(n)
   allData = replab.sym.CTData.instance(n);
   nData = allData(n);
   nParts = nData.nParts;
   CTunNormed = zeros(nParts);
   for row = 1:nParts
       partition = nData.partitions.list{row}.partition;
       lin = replab.sym.SymPoly.powerToMonomial(partition,n);
       CTunNormed(:,row) = lin.coeffs';
   end
   CT = round(grahamSchmidt(CTunNormed,@(x,y) round(sum(x.*y./nData.innerProdDenom))));
   CT = fliplr(CT);
   function mat = grahamSchmidt(mat,innerProd)
       d = nData.nParts;
       for i = 1:d-1
           for j = i+1:d
               mat(j,:) = mat(j,:)-innerProd(mat(j,:),mat(i,:))*mat(i,:);
           end
       end
   end
end

