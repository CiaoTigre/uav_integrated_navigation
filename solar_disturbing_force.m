function [ Ng ] = solar_disturbing_force( Lbe )
% 名称：solar disturbing force
% 功能：太阳摄动力
%
% Inputs:
%       Lbe: 载体质心地球质心间距离        m
%       where L denotes distance, b denote body, e denote earth.
% Outputs:
%       Ng: 摄动力加速度量级        g(9.8m/s^2)
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1

%%
% 万有引力常量 (N*m^2/kg^2)
G = 6.67259e-11;
% 地球半径 (m)
Re = 6.371e6;
% 地球公转轨道半径 (m)
Lse = 1.496e11;
% 太阳质量 (kg)
Ms = 2e30;

if nargin == 0
    delta_r = Re;
else
    delta_r = Lbe;
end

% disturbing acceleration
g_dst = 2*G*Ms*delta_r/Lse^3;
Ng = g_dst/9.8;

end

