function [ Cbn, angle ] = CoarseAlig( data, dt, B, L, H )
%UNTITLED Summary of this function goes here
% input
%   data -- [accX, accY, accZ, gyroX, gyroY, gyroZ]
%     dt -- interval, s
%      B -- latitude, degree
%      L -- longitude, degree
%      H -- height, m
% output
%  angle -- [ psai, fai, theta ]
% '航向角yaw','俯仰角pitch', '翻滚角roll'

% 准备数据
deg2arc = pi/180;
data = data/dt;
B = B*deg2arc;
L = L*deg2arc;
% 重力加速度(m/s^2)，地球自转角速度(rad/h)
g = 9.780325*(1 + 0.00530240*sin(B)*sin(B) - 0.00000582*sin(2*B)*sin(2*B));
g = g - 3.08e-6 * H;
omega_e = 7.292e-5 * 3600;

% 导航系[北东天]，正交规范化矢量
g_n = [0, 0, -g]';
omega_n = [omega_e*cos(B), 0, -omega_e*sin(B)]';
vn1 = orth(g_n);
vn2 = orth(cross(g_n, omega_n));
vn3 = orth(cross(cross(g_n, omega_n), g_n));
% 载体系，正交规范化矢量
g_b = orth([data(4), data(5), data(6)]');
omega_b = orth([data(1), data(2), data(3)]' * 3600);
vb1 = orth(g_b);
try
    vb2 = orth(cross(g_b, omega_b));
catch
    Cbn = zeros(3,3);
    angle = [NaN, NaN, NaN];
    return;
end
vb3 = orth(cross(cross(g_b, omega_b), g_b));

% 求解姿态角
%Anb = [g_b, omega_b, z_b] / [g_n, omega_n, z_n];
Anb = [vb1, vb2, vb3] / [vn1, vn2, vn3];
Cbn = Anb';

% psai, fai, theta
angle(1) = atan(Cbn(2,1)/Cbn(1,1)) / deg2arc;
angle(3) = atan(Cbn(3,2)/Cbn(3,3)) / deg2arc;
angle(2) = asin(-Cbn(3,1)) / deg2arc;