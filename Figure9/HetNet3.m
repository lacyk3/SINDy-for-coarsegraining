function rhs = HetNet3(t,X,dummy,nK,omega,KK,nF,alpha, KF,nR, params, KR, A, on)


n = nK + nF + nR;

a1=alpha(1);a2=alpha(2);a3=alpha(3); b=alpha(4);c=alpha(5);
aR = params(1);bR = params(2); cR = params(3);

theta = X(1:nK); v = X(nK+1:nK+nF); w=X(nK+nF+1:nK+2*nF);
x = X(nK+2*nF+1:nK+2*nF+nR); y = X(nK+2*nF+nR+1:nK+2*nF+2*nR); 
z =X(nK+2*nF+2*nR+1:nK+2*nF+3*nR);

coupling_theta = zeros(nK,1);
coupling_v=zeros(nF,1);
coupling_x= zeros(nR,1);
coupling_y = zeros(nR,1);

for j = 1:nK
    Ct = 0;
    for jj = 1:nK
        Ct = Ct + A(jj,j)*sin(theta(jj)-theta(j));
    end
    
    if on
    for jj = 1:nF
        Ct = Ct + A(jj,j)*sin(v(jj)-theta(j));
    end
    
    for jj = 1:nR
        Ct = Ct + A(jj,j)*sin(y(jj) - 10*cos(theta(j)));
    end
    end
    coupling_theta(j,1) = Ct;  
end
dt = omega + (KK/n)*coupling_theta;

for j=1:nF
  Cv = 0;
  
  for jj=1:nF
    Cv=Cv+ A(jj,j)*(v(j)-v(jj));
  end
  
  if on
  for jj = 1:nK
    Cv = Cv +A(jj,j)*(v(j) - cos(theta(jj)));
  end
%   
%   for jj = 1:nR
%     Cv = Cv + A(jj,j)*sin(v(j) - y(jj)/10);
%   end
  end
  coupling_v(j,1)=Cv;
  
end
dv = a3*v.^3+a2*v.^2+a1*v - w + (KF/n)*coupling_v;
dw = c*v-b*w;

for j = 1:nR    
  Cx = 0;
  Cy = 0; 

  if on
  for jj = 1:nK
    Cx = Cx +A(jj,j)*sin(x(j) - theta(jj));
    Cy = Cy +A(jj,j)*sin(y(j) - theta(jj));
  end
  
  for jj=1:nF
    Cx=Cx+ A(jj,j)*sin(v(jj)*10-x(j));
    Cy=Cy+ A(jj,j)*sin(v(jj)*10-y(j));
  end
  end
  
  for jj = 1:nR
    Cx = Cx + 3*A(jj,j)*(x(jj) - x(j));
    Cy = Cy + 3*A(jj,j)*(y(jj) - y(j));
  end
  
  coupling_x(j,1) = Cx;
  coupling_y(j,1) = Cy;
end

dx = -y - z + (KR/n)*coupling_x;
dy = x + aR*y + (KR/n)*coupling_y;
dz = bR + z.*(x-cR);

rhs=[dt;dv;dw;0.2*dx;0.2*dy;0.2*dz];