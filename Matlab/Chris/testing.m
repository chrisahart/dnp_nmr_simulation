clear
spins = [1/2 1/2 1/2];
nsp = length(spins);
Ix=sop(spins,'eex');
Iy=sop(spins,'eey');

tic
for ii = 0:1E4
    %temp = Ix * Iy;
    %temp = kron(Ix, Iy);
    temp = expm(-1i * Ix * 1E4);
end
toc
    