function CGTest
    S5 = replab.S(5);
    [~,mat1,irreps1]  =replab.sym.findAllCGCoeffs(S5,[4 1],[3 1 1],0);
    %Testing of floating point orthogonal decomp
    rep1 = kron(S5.irrep([4 1],'orthogonal'),S5.irrep([3 1 1],'orthogonal'));
    assert(isdiag(removeZeroErrors(mat1'*rep1.commutant.sample*mat1)));
    f = intPartList({[4 1],[3 2],[3 1 1],[2 2 1],[2 1 1 1]});
    isequal(f,irreps1);
    assert(isequal(f,irreps1));
    assert(norm(mat1*mat1'-eye(rep1.dimension))<1e-14);


    %Test of slightly higher dimensions and general orthogonal decomp
    S3 = replab.S(3);
    rep2 = kron(S3.naturalRep,S3.naturalRep,S3.naturalRep);
    [~,mat2,irreps2]  =replab.sym.findAllRepCoeffs(rep2,0);
    assert(nnz(removeZeroErrors(mat2'*rep2.commutant.sample*mat2))==203);
    f = intPartList({[3],[2 1],[1 1 1]});
    assert(isequal(f,irreps2));
    assert(norm(mat2*mat2'-eye(rep2.dimension))<1e-14);


    %Test of symbolic general decomp
    S6 = replab.S(6);
    rep2 = kron(S6.naturalRep);
    [~,mat3,~]  =replab.sym.findAllRepCoeffs(rep2,1);
    numerTest = [1    -1    -6    -3    -2    -3;...
     1    -1    -6    -3    -2     3;...
     1    -1    -6    -3     4     0;...
     1    -1    -6     9     0     0;...
     1    -1    24     0     0     0;...
     1     1     0     0     0     0];
 denomTest = [1     5    25    10     5     5;
     1     5    25    10     5     5;...
     1     5    25    10     5     1;...
     1     5    25    10     1     1;...
     1     5    25     1     1     1;...
     1     1     1     1     1     1];
 [n,d] = rat(mat3);
 assert(isequal(n,numerTest));
 assert(isequal(d,denomTest));

 %Symbolic test of CG decomp
[~,mat4,~]  = replab.sym.findAllCGCoeffs(replab.S(4),[2 2],[2 2],1);
numerTest =  [0    -3     0     3;...
 1     0     3     0;...
-1     0     3     0;...
 0     1     0     4];
denomTest = [1     4     1     1;...
 1     1     4     1;...
 1     1     4     1;...
     1     1     1     1];
 [n,d] = rat(mat4);
 assert(isequal(n,numerTest));
 assert(isequal(d,denomTest));
end

function cleanMat = removeZeroErrors(mat)
    dim = size(mat,1);
    cleanMat = zeros(dim,dim);
    cleanMat(abs(mat)>1e-10) = mat(abs(mat)>1e-10);
end

function partList = intPartList(cellOfParts)
    partList = cellfun(@(part) replab.sym.IntegerPartition(part),cellOfParts,'UniformOutput',0);
end