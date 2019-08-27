%Symmetric Group Representations

partition = [1 1];

SA = []; %Specht rep of (1 2)
SB = []; %Specht rep of (1 2 ... n)
NA = []; %semi-normal rep of (1 2)
NB = []; %semi-normal rep of (1 2 ... n)
OA = []; %orthogonal rep of (1 2)
OB = []; %orthogonal rep of (1 2 ... n)

if isequal(partition,[1]) %trivial
    
    SA = [1];
    SB = [1];
    NA = [1];
    NB = [1];
    OA = [1];
    OB = [1];
    
%S2
    
elseif isequal(partition,[2]) %trivial
    SA = [1];
    SB = [1];
    NA = [1];
    NB = [1];
    OA = [1];
    OB = [1];
    
elseif isequal(partition,[1 1]) %sign
    SA = [-1];
    SB = [1];
    NA = [-1];
    NB = [1];
    OA = [-1];
    OB = [1];
    
%S3
    
elseif isequal(partition,[3]) %trivial
    SA = [1];
    SB = [1];
    NA = [1];
    NB = [1];
    OA = [1];
    OB = [1];
    
elseif isequal(partition,[2 1]) %standard
    SA = [0 -1; -1 0];
    SB = [-1 1; -1 0];
    NA = [1 0; 0 -1];
    NB = [-1/2 3/2; -1/2 -1/2];
    OA = [1 0; 0 -1];
    OB = [-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2];
    
%S4
    
elseif isequal(partition,[4]) %trivial
    SA = [1];
    SB = [1];
    NA = [1];
    NB = [1];
    OA = [1];
    OB = [1];
    
elseif isequal(partition,[3 1])%standard
    SA = [0 -1 0; -1 0 0; 0 0 1];
    SB = [-1 1 -1; -1 0 0; 0 -1 0];
    NA = [1 0 0; 0 1 0; 0 0 -1];
    NB = [-1/3 4/3 0; -1/3 -1/6 3/2; -1/3 -1/6 -1/2];
    OA = [1 0 0; 0 1 0; 0 0 -1];
    OB = [-1/3 sqrt(2)/3 0; -sqrt(2)/3 -1/6 sqrt(3)/2; -sqrt(6)/3 -sqrt(3)/6 -1/2];
    
elseif isequal(partition,[2 1 1]) %product of standard & sign
    SA = [-1 0 0; 0 0 -1; 0 -1 0];
    SB = [1 -1 0; 1 0 -1; 1 0 0];
    NA = [1 0 0; 0 -1 0; 0 0 -1];
    NB = [1/2 -1/2 2; 1/2 1/6 -2/3; 0 2/3 1/3];
    OA = [1 0 0; 0 -1 0; 0 0 -1];
    OB = [1/2 -sqrt(3)/6 sqrt(6)/3; sqrt(3)/2 1/6 -sqrt(2)/3; 0 2*sqrt(2)/3 1/3];
    
elseif isequal(partition,[2 2]) %2D irrep of S4
    SA = [1 0; 1 -1];
    SB = [-1 1; 0 1];
    NA = [1 0; 0 -1];
    NB = [-1/2 -3/2; -1/2 1/2];
    OA = [1 0; 0 -1];
    OB = [-1/2 sqrt(3)/2; -sqrt(3)/2 1/2];
    
elseif isequal(partition,[1 1 1 1])%sign
    SA = [-1];
    SB = [-1];
    NA = [-1];
    NB = [-1];
    OA = [-1];
    OB = [-1];
    
end