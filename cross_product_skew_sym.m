function S = cross_product_skew_sym(vec)
%    
% Produces a skew-symmetric matrix
% representing the vector cross product.
%
% INPUT: vec = 3-by-1 vector.
%
%
%
% OUTPUT: S = 3-by-3 skew-symmetric matrix 
%             representing the vector cross product.

v1 = vec(1);
v2 = vec(2);
v3 = vec(3);

S = [0, -v3, v2;
    v3, 0, -v1;
    -v2, v1, 0];


end