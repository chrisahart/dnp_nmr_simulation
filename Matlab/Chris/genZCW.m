function[alpha,beta,gamma]=genZCW(M)
%%%%double c[3];
c(1)=1.;
c(2)=2.;
c(3)=1.;
N=getNumberZCW(M);
g2=getNumberZCW(M-2);
for m=1:N
    beta(m)=180/pi*acos(c(1)*(c(2)*mod(m/N,1.)-1.));
    alpha(m)=(180/pi)*2.*pi*(mod((m*g2/N),1.))/c(3); 
    gamma(m)=0.;
end

