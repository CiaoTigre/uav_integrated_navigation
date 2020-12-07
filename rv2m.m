function [ m ] = rv2m( rv )
%% **************************************************************
% 名称：Rotation Vector to direction cosine Matrix
% 功能：用等效旋转矢量计算方向余弦矩阵
%       e.g.: rv = PHIibi 表示i系绕转动轴uibi(uibb)转动了phi弧度与b系重合
% Inputs:
%       rv: 3×1 vector, 其模表示转动角度的大小，rv/norm(rv) 为转动轴在两个
%          坐标系下的坐标。 
% Outputs:
%       m: Cbi

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 程序@ 捷联惯导算法与组合导航原理讲义 P233
% 等效旋转矢量的模方：
nm2 = rv'*rv ;

if nm2 < 1.e-8
    a = 1 - nm2*(1/6 - nm2/120);
    b = 0.5 - nm2*(1/24 - nm2/720);
else
    nm = sqrt(nm2);
    a = sin(nm)/nm;
    b = (1 - cos(nm))/nm2;
end

VX = skew(rv);
m = eye(3) + a*VX + b*VX*VX;
end
