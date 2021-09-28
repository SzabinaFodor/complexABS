function x=generatecomplex(n)
rand('twister', n);
x=zeros( n, 1 );
for k=1:n
    x(k)=rand+1i*rand;
end
end