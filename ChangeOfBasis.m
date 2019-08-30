%FINDING MATRICES FROM SPECHT TO ORTHOGONAL REPRESENTATION

%S3 STANDARD

%generators
a = [0 -1; -1 0]; %(1 2)
b = [-1 1; -1 0]; %(1 2 3)

ab = a*b; %(2 3)
b2 = b*b; %(1 3 2)
ab2 = a*b*b; %(1 3)
id = a*a;

R = (1/6)*(a'*a + b'*b + ab'*ab + b2'*b2 + ab2'*ab2 + id);

A_S3_S = chol(R); %change of basis matrix for S3 standard rep Specht to orthogonal 

%%
%S4 STANDARD REP SPECHT
S4 = replab.S(4);
images = {};
a = [-1 1 -1; -1 0 0; 0 -1 0]; %image of (1 2 3 4)
b = [0 -1 0; -1 0 0; 0 0 1]; %image of (1 2)

images = {};

images{1,1} = b^2; %e
images{1,2} = a^2*b*a^2; %(3 4)
images{1,3} = a*b*a^3; %(2 3) 
images{1,4} = a*b; %(2 3 4)
images{1,5} = a^3*b; %(2 4 3)
images{1,6} = a*b*a*b*a^3*b*a^3; %(2 4)
images{1,7} = b; %(1 2)
images{1,8} = a^2*b*a^2*b; %(1 2)(3 4)
images{1,9} = b*a*b*a^3; %(1 2 3)
images{1,10} = a; %(1 2 3 4)
images{1,11} = a*b*a; %(1 2 4 3)
images{1,12} = a^2*b*a^3; %(1 2 4)
images{1,13} = a*b*a^3*b; %(1 3 2)
images{1,14} = b*a*b; %(1 3 4 2)
images{1,15} = b*a*b*a^3*b; %(1 3)
images{1,16} = a*b; %(1 3 4)
images{1,17} = a^2; %(1 3)(2 4)
images{1,18} = a*b*a^2*b*a*b; %(1 3 2 4)
images{1,19} = a^2*b*a^3*b*a^3*b; %(1 4 3 2)
images{1,20} = a^2*b*a*b*a^3; %(1 4 2)
images{1,21} = b*a^3; %(1 4 3)
images{1,22} = a^2*b*a^3*b; %(1 4)
images{1,23} = a^2*b; %(1 4 2 3)
images{1,24} = a*b*a^2*b*a; %(1 4)(2 3)

R = zeros(3);
for i = 1:24
    R = R + images{1,i}'*images{1,i};
end
R = R/24;
A_S4_S = chol(R); %change of basis matrix S4 standard rep Specht to orthogonal 

%%
%S4 2D REP SPECHT
S4 = replab.S(4);
images = {};
a = [-1 1; 0 1]; %image of (1 2 3 4)
b = [1 0; 1 -1]; %image of (1 2)

images = {};

images{1,1} = b^2; %e
images{1,2} = a^2*b*a^2; %(3 4)
images{1,3} = a*b*a^3; %(2 3) 
images{1,4} = a*b; %(2 3 4)
images{1,5} = a^3*b; %(2 4 3)
images{1,6} = a*b*a*b*a^3*b*a^3; %(2 4)
images{1,7} = b; %(1 2)
images{1,8} = a^2*b*a^2*b; %(1 2)(3 4)
images{1,9} = b*a*b*a^3; %(1 2 3)
images{1,10} = a; %(1 2 3 4)
images{1,11} = a*b*a; %(1 2 4 3)
images{1,12} = a^2*b*a^3; %(1 2 4)
images{1,13} = a*b*a^3*b; %(1 3 2)
images{1,14} = b*a*b; %(1 3 4 2)
images{1,15} = b*a*b*a^3*b; %(1 3)
images{1,16} = a*b; %(1 3 4)
images{1,17} = a^2; %(1 3)(2 4)
images{1,18} = a*b*a^2*b*a*b; %(1 3 2 4)
images{1,19} = a^2*b*a^3*b*a^3*b; %(1 4 3 2)
images{1,20} = a^2*b*a*b*a^3; %(1 4 2)
images{1,21} = b*a^3; %(1 4 3)
images{1,22} = a^2*b*a^3*b; %(1 4)
images{1,23} = a^2*b; %(1 4 2 3)
images{1,24} = a*b*a^2*b*a; %(1 4)(2 3)

R = zeros(2);
for i = 1:24
    R = R + images{1,i}'*images{1,i};
end
R = R/24;
A_S4_T = chol(R); %change of basis matrix S4 2D rep Specht to orthogonal 

%testing that the matrices generated using this change of basis is
%orthogonal
for j = 1:24
    U = A_S4_T*images{1,j}*inv(A_S4_T);
    U*U'
end

%%
%S4 PRODUCT OF STANDARD AND SIGN REP SPECHT
S4 = replab.S(4);
images = {};
a = [1 -1 0; 1 0 -1; 1 0 0]; %image of (1 2 3 4)
b = [-1 0 0; 0 0 -1; 0 -1 0]; %image of (1 2)

images = {};

images{1,1} = b^2; %e
images{1,2} = a^2*b*a^2; %(3 4)
images{1,3} = a*b*a^3; %(2 3) 
images{1,4} = a*b; %(2 3 4)
images{1,5} = a^3*b; %(2 4 3)
images{1,6} = a*b*a*b*a^3*b*a^3; %(2 4)
images{1,7} = b; %(1 2)
images{1,8} = a^2*b*a^2*b; %(1 2)(3 4)
images{1,9} = b*a*b*a^3; %(1 2 3)
images{1,10} = a; %(1 2 3 4)
images{1,11} = a*b*a; %(1 2 4 3)
images{1,12} = a^2*b*a^3; %(1 2 4)
images{1,13} = a*b*a^3*b; %(1 3 2)
images{1,14} = b*a*b; %(1 3 4 2)
images{1,15} = b*a*b*a^3*b; %(1 3)
images{1,16} = a*b; %(1 3 4)
images{1,17} = a^2; %(1 3)(2 4)
images{1,18} = a*b*a^2*b*a*b; %(1 3 2 4)
images{1,19} = a^2*b*a^3*b*a^3*b; %(1 4 3 2)
images{1,20} = a^2*b*a*b*a^3; %(1 4 2)
images{1,21} = b*a^3; %(1 4 3)
images{1,22} = a^2*b*a^3*b; %(1 4)
images{1,23} = a^2*b; %(1 4 2 3)
images{1,24} = a*b*a^2*b*a; %(1 4)(2 3)

R = zeros(3);
for i = 1:24
    R = R + images{1,i}'*images{1,i};
end
R = R/24;
A_S4_Ss = chol(R); %change of basis matrix S4 product of standard and sign rep Specht to orthogonal 

%testing that the matrices generated using this change of basis is
%orthogonal
for j = 1:24
    U = A_S4_Ss*images{1,j}*inv(A_S4_Ss);
    U*U'
end

%%
%S5 STANDARD REP SPECHT
S5 = replab.S(5);
images = {};
S = replab.NURepByImages(S5,'R', 4, {[-1 1 -1 1; -1 0 0 0; 0 -1 0 0; 0 0 -1 0] [0 -1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1]},{inv([-1 1 -1 1; -1 0 0 0; 0 -1 0 0; 0 0 -1 0]) inv([0 -1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1])});
for i = 1:120
    images{1,i} = image(S,S5.elements.at(i));
end

R = zeros(4);
for j = 1:120
    R = R + images{1,j}'*images{1,j};
end
R = R/120;
A_S5_S = chol(R);

%%
%S5 STANDARD x SIGN REP SPECHT
S5 = replab.S(5);
images = {};
Ss = replab.NURepByImages(S5, 'R', 4, {[-1 1 0 0; -1 0 1 0; -1 0 0 1; -1 0 0 0] [-1 0 0 0; 0 -1 0 0; 0 0 0 -1; 0 0 -1 0]}, {inv([-1 1 0 0; -1 0 1 0; -1 0 0 1; -1 0 0 0]) inv([-1 0 0 0; 0 -1 0 0; 0 0 0 -1; 0 0 -1 0])});
for i = 1:120
    images{1,i} = image(Ss,S5.elements.at(i));
end

R = zeros(4);
for j = 1:120
    R = R + images{1,j}'*images{1,j};
end
R = R/120;
A_S5_Ss = chol(R);

%%
%S5 5D (3 + 2) REP SPECHT
S5 = replab.S(5);
images = {};
f = replab.NURepByImages(S5, 'R', 5, {[-1 1 -1 0 0; 0 0 0 1 -1; 0 0 1 0 -1; 1 0 0 0 0; 1 0 1 0 0] [1 0 0 0 0; 0 0 0 -1 0; -1 0 0 0 -1; 0 -1 0 0 0; -1 0 -1 0 0]}, {inv([-1 1 -1 0 0; 0 0 0 1 -1; 0 0 1 0 -1; 1 0 0 0 0; 1 0 1 0 0]) inv([1 0 0 0 0; 0 0 0 -1 0; -1 0 0 0 -1; 0 -1 0 0 0; -1 0 -1 0 0])});
for i = 1:120
    images{1,i} = image(f,S5.elements.at(i));
end

R = zeros(5);
for j = 1:120
    R = R + images{1,j}'*images{1,j};
end
R = R/120;
A_S5_f = chol(R);

%%
%S5 5D (2 + 2 + 1) REP SPECHT
S5 = replab.S(5);
images = {};
F = replab.NURepByImages(S5, 'R', 5, {[1 -1 -1 1 0; 0 -1 -1 0 1; 1 -1 0 0 0; 0 -1 0 0 0; 1 -1 -1 0 0] [0 0 -1 0 0; 0 0 0 -1 0; -1 0 0 0 0; 0 -1 0 0 0; 0 0 0 0 -1]}, {inv([1 -1 -1 1 0; 0 -1 -1 0 1; 1 -1 0 0 0; 0 -1 0 0 0; 1 -1 -1 0 0]) inv([0 0 -1 0 0; 0 0 0 -1 0; -1 0 0 0 0; 0 -1 0 0 0; 0 0 0 0 -1])});
for i = 1:120
    images{1,i} = image(F,S5.elements.at(i));
end

R = zeros(5);
for j = 1:120
    R = R + images{1,j}'*images{1,j};
end
R = R/120;
A_S5_F = chol(R);

%%
%S5 EXTERIOR PRODUCT OF STANDARD REP SPECHT
S5 = replab.S(5);
images = {};
E = replab.NURepByImages(S5, 'R', 6, {[1 -1 1 0 0 0; 1 0 0 -1 1 0; 0 1 0 -1 0 1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 1 0 0] [-1 0 0 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 0 0 1]}, {inv([1 -1 1 0 0 0; 1 0 0 -1 1 0; 0 1 0 -1 0 1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 1 0 0]) inv([-1 0 0 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 0 0 1])});
for i = 1:120
    images{1,i} = image(E, S5.elements.at(i));
end

R = zeros(6);
for j = 1:120
    R = R + images{1,j}'*images{1,j};
end
R = R/120;
A_S5_E = chol(R);