function interleaved=interleaveCols(a,b)
%interleaves columns of matricies or cell arrays of a and b
%a and b must be the same size!
% % create two matrices.
%a = [11 13; 12 14]
%b = [21 23; 22 24]
%col_interleave=interleaveCols(a,b);
% col_interleave =
% 
%     11    21    13    23
%     12    22    14    24

a = a';
b = b';
interleaved = reshape([a(:) b(:)]',2*size(a,1), [])';



end