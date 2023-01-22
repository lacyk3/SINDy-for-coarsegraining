function rhs=FHN_rhs(t,V,dummy,alpha,n,K,A,on)
a1=alpha(1);a2=alpha(2);a3=alpha(3); b=alpha(4);c=alpha(5);
v = V(1:n);w=V(n+1:end);
coupling_v=zeros(n,1);
coupling_w = zeros(n,1);


for j=1:n
  Cv=0;
  for jj=1:n
    Cv=Cv+on*A(j,jj)*(v(j)-v(jj));
  end
  coupling_v(j,1)=Cv;
end
dv = a3*v.^3+a2*v.^2+a1*v - w + (K/n)*coupling_v;
dw = c*v-b*w;
rhs=[dv;dw];