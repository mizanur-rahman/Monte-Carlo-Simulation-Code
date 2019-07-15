for a=1:100000
    v(a)=randn(1,1)*2+10;
end
hist(v,30)