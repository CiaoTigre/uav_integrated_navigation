function [ Cbn ] = a2mat( att )
% Function：Attitude to dirction cosine MATrix 
% 完全按照严恭敏老师书上复制的版本
%   
% Inputs:
%       att:	姿态角(pitch-θ,roll-Φ,yaw-ψ), unit:rad
%
% Outputs:
%       Cbn:    从b系到n系的坐标变换矩阵, 也是从n系到b系的坐标系转换矩阵
%
% Reference:
%   - @ 捷联惯导算法与组合导航原理 P231
%
% Notes:
%   - 变量命名时上角标在前
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1

s = sin(att);
c = cos(att);

si = s(1);
sj = s(2);
sk = s(3);

ci = c(1);
cj = c(2);
ck = c(3);

Cbn = [cj*ck - si*sj*sk, -ci*sk, sj*ck + si*cj*sk
       cj*sk + si*sj*ck,  ci*ck, sj*sk - si*cj*ck
                 -ci*sj,     si,            ci*cj];
end
