function [ rv ] = q2rv( q )
% 名称：Quaternion to Rotation Vector
% 功能：将姿态四元数转化为等效旋转矢量
%
% Inputs:
%       q: 
% Outputs:
%       rv: 

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 原理@ 捷联惯导算法与组合导航原理讲义 P233
% 首先将四元数转化为标量非负四元数
if q(1) < 0
    q = -q;
end

% 等效旋转矢量模值的一半
nmhalf = acos(q(1));

if nmhalf > 1e-20
    b = 2*nmhalf/sin(nmhalf);
else
    b = 2;
end

rv = b*q(2:4);

end
