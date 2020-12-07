function [ dpos ] = deltapos( pos )
% 名称：delta position
% 功能：求位置增量
% 程序@ 捷联惯导系统与组合导航原理 P236

% Inputs:
%       pos: N×3 矩阵，每一行保存一个时刻的位置
% Output:
%       dpos: N×3 矩阵，每一行对应着pos中对应行元素与第一个位置pos(, :)之差

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% cos(latitude)
cl = cos(pos(:, 1));
% 地球半径 m，此处仅考虑地球平均半径
Re = 6378137;

dpos = [(pos(:, 1) - pos(1, 1))*Re, (pos(:, 2) - pos(1, 2)).*cl*Re,...
        pos(:, 3) - pos(1, 3) ];

end
