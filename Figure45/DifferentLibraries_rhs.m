function rhs = DifferentLibraries_rhs(t,y,~,coeffs,n1,p1,u1,p2,u2)
v = y(1:n1).'; w = y(n1+1:end).';
Lib1 = poolData(v,n1,p1,u1);
Lib2 = poolData(w,numel(y)-n1,p2,u2);
Library = [Lib1, Lib2(:,2:end)];
rhs = (Library*coeffs)';