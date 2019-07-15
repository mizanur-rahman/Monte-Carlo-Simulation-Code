s=[9.802 9.774;5 5];
se=[.037 .019;1 1];
for i=1:length(s)
    w(:,i)=1./(se(:,i).*se(:,i));
  
end
n=w.*s;
nu=sum(n');
den=sum(w');
mean=nu./den;
error=1./sqrt(den);