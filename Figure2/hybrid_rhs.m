function dy = hybrid_rhs(t,y,~,weights,partition)
% -- apply function library
yPool = Library(y);

% -- determine which dynamic region y lies in
m = length(y);
p = 0;
for i = 1:size(partition,1)
    if y(1) >= partition(i,1) && y(1) <= partition(i,1)+ partition(i,3)...
        && y(2) >= partition(i,2) && y(2) <= partition(i,2)+ partition(i,4)
    p = i;
    end
end
% -- select appropriate weights
if p > 0 
    ahat = weights(:,(p-1)*m+1:p*m);
else
    ahat = weights(:,end-(m-1):end);
end

dy = (yPool*ahat)';