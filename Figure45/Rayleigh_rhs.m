function rhs=Rayleigh_rhs(t,V,dummy,alpha,n,K,A,on)

v = V(1:n);w=V(n+1:end);
coupling=zeros(n,1);

for j=1:n
  C=0;
  for jj=1:n
    C=C+on*A(j,jj)*(-1-(v(j) - v(jj)).^2).*(w(j)-w(jj));
  end
  coupling(j,1)=C;
end
dv = w;
dw = 1./alpha.*(w - w.^3/3 - v+(K/n)*coupling);
rhs=[dv;dw];