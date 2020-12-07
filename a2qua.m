function [ qbn ] = a2qua( att )
% Function：Attitude to Quaternion
%   
% Inputs:
%       att:	姿态角(pitch-θ,roll-Φ,yaw-ψ), unit:rad
%
% Outputs:
%       qbn:    从b系到n系的四元数
%
% Reference:
%   - @ 捷联惯导算法与组合导航原理 P232  @原理 P247
%
% Notes:
%   - 变量命名时上角标在前
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1

s = sin(att/2);
c = cos(att/2);

si = s(1);  sj = s(2);  sk = s(3);
ci = c(1);  cj = c(2);  ck = c(3);

qbn = [ci*cj*ck - si*sj*sk
       si*cj*ck - ci*sj*sk
       ci*sj*ck + si*cj*sk
       ci*cj*sk + si*sj*ck];

end

