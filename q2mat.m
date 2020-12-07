function [ Cbn ] = q2mat( qbn )
% 名称：Quaternion to directin cosine MATrix
% 功能：用姿态四元数计算方向余弦矩阵
%
% Inputs:
%       qbn: rn = qbn*rb*qbn'
% Outputs:
%       Cbn: 从b系到n系的坐标转换矩阵

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 书上没用采用将四元数分解为标量部分+矢量部分的计算方式，而是采用直接用四元数
% 中各个元素来计算Cbn中每一个元素的计算方式。
% 程序@ 捷联惯导算法与组合导航讲义 P232
q11 = qbn(1)*qbn(1);    q12 = qbn(1)*qbn(2);    q13 = qbn(1)*qbn(3);    q14 = qbn(1)*qbn(4);
q22 = qbn(2)*qbn(2);    q23 = qbn(2)*qbn(3);    q24 = qbn(2)*qbn(4);
q33 = qbn(3)*qbn(3);    q34 = qbn(3)*qbn(4);
q44 = qbn(4)*qbn(4);
    
Cbn = [q11 + q22 - q33 - q44,           2*(q23 - q14),            2*(q24 + q13)
               2*(q23 + q14),   q11 - q22 + q33 - q44,            2*(q34 - q12)
               2*(q24 - q13),           2*(q34 + q12),   q11 - q22 - q33 + q44];

end

