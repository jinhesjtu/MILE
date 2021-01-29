function [ output_args ] = norml(M, ele )
%NORML  scalar ambiguty removal

[row, col] = size(M);
N = ones(row,1)*M(ele,:);
output_args = M./N;


end

