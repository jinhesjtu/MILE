%% Calculate Root Mean Squared Error for each row.
%% output B = sqrt(mean(abs(A-a)).^2).
function B = rmse(A,a);
[row,col]=size(A);
if length(a) == 1;
    a = a*ones(1,col);
end

if row == 1;
    B = sqrt(mean(abs(A - a).^2));
elseif col == 1 
    B = sqrt(mean(abs(A - a).^2));
else
    for num = 1:col
        B(:,num) = sqrt(mean(abs(A(:,num) - a(num)).^2));
    end
end