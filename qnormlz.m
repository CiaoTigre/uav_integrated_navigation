function [ qbn ] = qnormlz( qbn )
% 名称：Quaternion Normlization
% 功能：将四元数标准化
%
% Inputs:
%       qbn: 未标准化四元数
% Outputs:
%       qbn: 标准化过后的四元数

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 输入四元数的模方
nm2 = qbn'*qbn;

if nm2 < 1e-6
    qbn = [1, 0, 0, 0]';
else
    qbn = qbn/sqrt(nm2);
end

end
