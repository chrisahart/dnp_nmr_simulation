clear
spins = [1/2 1/2 1/2];
nsp = length(spins);
Ix=sop(spins,'eex');
Iy=sop(spins,'eey');

% A=[1, 2, 3, 4;
%     5, 6, 7, 8;
%     9, 10, 11, 12;
%     13, 14, 15, 16];
% 
% B=[1, 5, 9, 13;
%     2, 6, 10, 14;
%     3, 7, 11, 15;
%     4, 8, 12, 16];
% 
% Ix = kron(A, B);
% Iy = kron(B, A);

tic
for ii = 0:1E6
    temp = Ix * Iy;
    %temp = kron(Ix, Iy);
    % temp = expm(-1i * Ix * 1E4);
end
toc
    