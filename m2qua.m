function [ qbn ] = m2qua( Cbn )
% 名称：direction cosine Matrix to QUAternion
% 功能：将方向余弦矩阵转化为姿态四元数
%
% Inputs:
%       Cbn: 从b系到n系的坐标转换矩阵
% Outputs:
%       qbn: rn = qbn*rb*qbn';

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 程序@ 捷联惯导算法与组合导航讲义 P232
% 原理@  P246
% 
% 转动四元数具有多值性，为了确定四元数中各个元素的正负号，
% 需要首先确定一个元素然后求解其他元素，在
% 实际计算中，不妨先确定绝对值大的元素(目的是确保该元素不为0)，并将
% 该元素符号取为正，在计算其他元素。
C11 = Cbn(1, 1);    C12 = Cbn(1, 2);    C13 = Cbn(1, 3);
C21 = Cbn(2, 1);    C22 = Cbn(2, 2);    C23 = Cbn(2, 3);
C31 = Cbn(3, 1);    C32 = Cbn(3, 2);    C33 = Cbn(3, 3);

if C11 >= C22 + C33
    q1 = 0.5*sqrt(1 + C11 - C22 - C33);
    q0 = (C32 - C23)/(4*q1);
    q2 = (C12 + C21)/(4*q1);
    q3 = (C13 + C31)/(4*q1);
elseif C22 >= C11 + C33
    q2 = 0.5*sqrt(1 - C11 + C22 - C33);
    q0 = (C13 - C31)/(4*q2);
    q1 = (C12 + C21)/(4*q2);
    q3 = (C23 + C32)/(4*q2);
elseif C33 >= C11 + C22
    q3 = 0.5*sqrt(1 - C11 - C22 + C33);
    q0 = (C21 - C12)/(4*q3);
    q1 = (C13 + C31)/(4*q3);
    q2 = (C23 + C32)/(4*q3);
else
    q0 = 0.5*sqrt(1 + C11 + C22 + C33);
    q1 = (C32 - C23)/(4*q0);
    q2 = (C13 - C31)/(4*q0);
    q3 = (C21 - C12)/(4*q0); 
    
end

qbn = [q0, q1, q2, q3]';

end

