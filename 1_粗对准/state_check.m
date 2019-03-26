function [ data1, n ] = state_check( data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% data ax, ay, az, gx, gy, gz
% output: checked data

data1 = data;
% -------------- remove bad data ---------------------------------
for i = 1:6
    if i == 5
        % don't average gy
        continue;
    end
    idata = data1(:, i);
    len = size(idata);
    len = len(1);
    ave = sum(idata) / len;
    res = idata - ave;
    rms = res' * res;
    rms = sqrt(rms / len);
    index_x = find(abs(res) < 3*rms);
    data1 = data1(index_x, :);
end

% -------------- average -----------------------------------------
n = size(data1);
n = n(1);
data1 = sum(data1);
data1 = data1/n;
n = 1;


