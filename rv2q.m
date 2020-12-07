function [ qbn ] = rv2q( rv )
% 名称：Rotation Vector to Quaternion
% 功能：将等效旋转矢量转化为姿态四元数
%
% Inputs:
%        rv: 3×1 vector, 其模表示转动角度的大小，rv/norm(rv) 为转动轴在两个
%            坐标系下的坐标。 
% Outputs:
%       qbn: rn = qbn*rb*qbn'

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 程序@ 捷联惯导算法与组合导航原理讲义 P233
% 等效旋转矢量的模方：
nm2 = rv'*rv ;

if nm2 < 1e-8       % 如果模方很小，可以用泰勒级数展开求前几项三角函数
    q0 = 1 - nm2*(1/8 - nm2/384);
    s = 1/2 - nm2*(1/48 - nm2/3840);
else
    nm = sqrt(nm2);
    q0 = cos(nm/2);
    s = sin(nm/2)/nm;
end

qbn = [q0, s*rv']';

end

