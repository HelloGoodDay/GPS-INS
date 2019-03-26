function [ data1, nepoch ] = fitting( data )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

data1 = data;
nepoch = size(data);
nepoch = nepoch(1);
% -------------- fitting gy data ---------------------------------
gy = data1(:, 5);
B = ones(nepoch, 2);
L = zeros(nepoch, 1);
for i  = 1: nepoch
    B(i, 2) = i - 1;
    L(i, 1) = gy(i);
end
x = inv(B'*B) * (B'*L);

% ----------- correct --------------------------------------------
for i  = 1: nepoch
    gy_cor(i) = gy(i) - i * x(2);
end
plot(gy);
hold on;
plot(gy_cor, 'r');

% ----------------------------------------------------------------
data1(:, 5) = gy_cor;
end

