G = replab.AbstractGroup.parsePresentation('<a, x | a^8 = x^2 = 1, x a x^-1 = a^-1 >');
sqrt2 = replab.cyclotomic.sqrt(2);
img_a = [1/sqrt2 -1/sqrt2
         1/sqrt2 1/sqrt2];
img_x = [1 0; 0 -1];
rep = G.repByImages('R', 2, {img_a img_x});
