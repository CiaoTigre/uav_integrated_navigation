function [ phi ] = qq2phi( qbp, qbn )
% 名称：Quaternion Quaternion to PHI
% 功能：计算含有误差的姿态四元数和真值之间对应的等效旋转矢量
%
% Inputs:
%       qbp: 含有误差的姿态四元数
%       qbn: 真实姿态四元数
% Outputs:
%       phi: n系到p系的等效旋转矢量

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
qerr = qmul(qbn, qconj(qbp));
phi = q2rv(qerr);

end

