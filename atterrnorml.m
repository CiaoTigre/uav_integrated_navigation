function [ att_err_normalized ] = atterrnorml( att_err )
% 名称: Attitude Error Normalization (version 1.0)
% Function: 姿态误差标准化 (目前仅修正航向角误差)
%
% Inputs:
%       att_err: 未标准化的姿态角误差
% Outputs:
%       att_err_normalized: 标准化过后的姿态角误差

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1


% 分别提取三个姿态角误差
pitch_err = att_err(1);
roll_err = att_err(2);
yaw_err = att_err(3);

% pitch 角度范围为: [-pi/2, pi/2]
% roll 角度范围为: (-pi, pi]
% yaw 角度范围为: [0, 2*pi)

% 直接计算出的航向角误差∈(-2*pi, 2*pi),但我们希望航向角误差∈(-pi, pi]
% 如果航向角误差绝对值大于pi，就令其标准化
if norm(yaw_err) > pi
    % 如果航向角误差为正
    if sign(yaw_err) == 1
        yaw_err = yaw_err - 2*pi;
    else
        % 航向角误差大于pi，且符号为负
        yaw_err = yaw_err + 2*pi;
    end
end

att_err_normalized = [pitch_err, roll_err, yaw_err]';

end

