function[sum]=getNumberZCW(M)
gM=5;
gMminus1=3; 
sum=5; 
for j=0:1:M
    sum=(gM+gMminus1); 
    gMminus1=gM;
    gM=sum;
end