function [ qbn ] = qdelphi( qbp, phi )
% 名称：Quaternion DELete PHI
% 功能：四元数减失准角误差，将误差从含有误差的姿态四元数种剔除，以获得真值。
%
% 因为采用n'来表示计算导航坐标系会在过渡矩阵的表示过程中产生歧义，因此用p来
% 代表n'系。
% Inputs:
%       qbp: 从b系到计算导航系p的姿态四元数。
%       phi: 由n系到p系这个转动所对应的等效旋转矢量
%       @ 捷联惯导算法与组合导航原理 P234 原文：从真实导航系n到计算导航系p，
%         的失准角为phi，反之，从p系到n系的失准角为-phi，若将-phi视为等效旋
%         转矢量，则与其对应的变换四元数为Qnp。
%         这就是说本书中将失准角(312 欧拉角)视为等效旋转矢量。
% Outputs:
%       qbn: 从b系到真实导航系的姿态四元数

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
qbn = qmul(rv2q(phi), qbp);

end

