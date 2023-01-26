function out = HetNet_rhs(t,x,dummy,nK,omega,KK,nF,alpha, KF,A)
n = nK + nF;
a1=alpha(1);a2=alpha(2);a3=alpha(3); b=alpha(4);c=alpha(5);
theta = x(1:nK); v = x(nK+1:n); w=x(n+1:end);

coupling_theta = zeros(nK,1);
coupling_v=zeros(nF,1);

for j = 1:nK
    Ct = 0;
    for jj = 1:nK
        Ct = Ct + A(jj,j)*sin(theta(jj)-theta(j));
    end
    
    for jj = 1:nF
        Ct = Ct + A(jj,j)*sin(v(jj)-theta(j));
    end
    
    coupling_theta(j,1) = Ct;  
end
dt = omega + (KK/n)*coupling_theta;

for j=1:nF
  Cv = 0;
  
  for jj = 1:nK
    Cv = Cv +A(jj,j)*(v(j) - sin(theta(jj)));
  end
  
  for jj=1:nF
    Cv=Cv+ A(jj,j)*(v(j)-v(jj));
  end
  coupling_v(j,1)=Cv;
  
end
dv = a3*v.^3+a2*v.^2+a1*v - w + (KF/n)*coupling_v;
dw = c*v-b*w;
out=[dt;dv;dw];