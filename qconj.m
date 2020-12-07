function [ qout ] = qconj( qin )
% 名称：CONJugate Quaternion
% 功能：输出共轭四元数
%
% Inputs:
%       qin: 输入四元数
% Outputs:
%       qout: 输入四运输的共轭

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%% 
qout = zeros(4, 1);

qout(1) = qin(1);
qout(2:4) = -qin(2:4);

end

