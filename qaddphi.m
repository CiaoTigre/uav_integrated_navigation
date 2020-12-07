function [ qbp ] = qaddphi( qbn, phi )
% 名称：Quaternion ADD PHI
% 功能：四元数加失准角误差,用来将真实姿态四元数加上误差来模拟含有误差的姿态
%      四元数。
%
% 因为采用n'来表示计算的导航坐标系会在过渡矩阵的表示过程中产生歧义，因此用p来
% 代表n'系。
% Cbp = Cnp*Cbn, 同理有：qbp = qnp*qbn
% Inputs:
%       qbn: 真实姿态四元数
%       phi: 从n系到p系的等效旋转矢量
% Outputs:
%       pbp: 含有误差的姿态四元数

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 程序@ 捷联惯导算法与组合导航原理讲义 P234
qbp = qmul(rv2q(-phi), qbn);

% 另一种计算方法：
% qpn = rv2q(phi);
% qnp = qconj(qpn);
% 
% qbp = qmul(qnp, qbn);

end

