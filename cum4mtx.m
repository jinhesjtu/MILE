function R = cum4mtx(x,y,X,Y);
% function R = cum4mtx(x,y,X,Y); 
% calculate cumulant matrix where X,Y is in matrix form.

[Sen,Snap] = size(X);

XX = ones(Sen,1)*x.*X; YY = ones(Sen,1)*y.*Y;
R = Snap*XX*YY'  - x*y'*X*Y'  -  (X*x.') * (conj(y)*conj(Y).') - ((conj(Y)*x.') * ((conj(y))*(X).') ).';