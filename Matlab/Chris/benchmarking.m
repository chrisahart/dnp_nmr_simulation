size = 80;

for a=1:size
    for b=1:size
        A(a, b) = 1E-9 + (a * b);
    end
    
end

tic
for test=1:1E4
    t=expm(-1i*A*1E-8);
end
toc
