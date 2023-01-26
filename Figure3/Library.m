function A = Library(y)
X = y(1); Y = y(2);
A = [X.^0,X,X.^2,X.^3, X./(1-X.^2), X.*(1+X.^2)./(1-X.^2).^3,...
    Y.^0, Y, Y.*X, Y.*X.^2, Y.*X.^3, Y.*X.^4, Y.*X.^5];
end
