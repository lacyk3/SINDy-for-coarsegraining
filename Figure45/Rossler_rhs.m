function rhs=Rossler_rhs(t,V,dummy,alpha,n,K,A,on)
a = alpha(1); b= alpha(2); c = alpha(3);
x = V(1:n);y=V(n+1:2*n);z=V(2*n+1:end);
coup=zeros(n,3);

for j=1:n
  C=0;
  for jj=1:n
    C=C+on*A(j,jj)*[sin(x(jj)-x(j)), sin(y(jj)-y(j)), sin(z(jj)-z(j))];
  end
  coup(j,:)=C;
end

dx = -y - z + on*(K/n)*coup(:,1);

dy = x + a*y + on*(K/n)*coup(:,2);

dz = b + z.*(x-c) + on*(K/n)*coup(:,3);

rhs=[dx;dy;dz];