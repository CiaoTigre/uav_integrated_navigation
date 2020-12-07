function [ q ] = qmul( q1, q2 )
% 名称：Quaternion MULtiplication
% 功能：四元数乘法
%
% Inputs:
%       q1, q2: 两个相乘的四元数
% Outputs:
%       q: 相乘之后的结果

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% *** 初始化 ***
% 将输出四元数转化为列向量形式
if size(q1, 1) == 1
    q1 = q1';
end
if size(q2, 1) == 1
    q2 = q2';
end

% 如果输入四元数只有三个元素则将其视为零标量四元数
if length(q1) == 3
    q1 = [0, q1']';
end
if length(q2) == 3
    q2 = [0, q2']';
end
% *** 初始化 结束 ***

MQ1 = [q1(1), -q1(2), -q1(3), -q1(4)
       q1(2),  q1(1), -q1(4),  q1(3)
       q1(3),  q1(4),  q1(1), -q1(2)
       q1(4), -q1(3),  q1(2),  q1(1)];

q = MQ1*q2;

end

