function [ att ] = q2att( qbn )
% 名称：Attitude to QUAternion
% 功能：将姿态四元数转化为姿态角
%
% Inputs:
%       qbn: rn = qbn*rb*qbn'
% Outputs:
%       att: 有三个分量，依次为俯仰theta、横滚gamma和航向psi。 航向角北偏
%       西为正,取值为（-pi, pi]

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
att = m2att(q2mat(qbn));

end

