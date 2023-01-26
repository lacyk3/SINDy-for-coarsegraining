function A = Library(y)
X = y(1); Y = y(2);
A = [ X.^0, X, Y,  X.^2, Y.^2,...
                Y./((1./(Y-1).^2)+(2/3)*(1./(Y-1))-(2/3)*(1./(Y+2))),...
                Y./((1./(Y+1).^2)-(2/3)*(1./(Y+1))+(2/3)*(1./(Y-2))),...
                log(abs(X-1)), log(abs(X+1)),  ...
                1./(X-1), 1./(X+1), X./(1-X.^2),...
                Y./(1-Y.^2)];

end