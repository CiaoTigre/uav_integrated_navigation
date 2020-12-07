function [ q ] = qmul( q1, q2 )
% ���ƣ�Quaternion MULtiplication
% ���ܣ���Ԫ���˷�
%
% Inputs:
%       q1, q2: ������˵���Ԫ��
% Outputs:
%       q: ���֮��Ľ��

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% *** ��ʼ�� ***
% �������Ԫ��ת��Ϊ��������ʽ
if size(q1, 1) == 1
    q1 = q1';
end
if size(q2, 1) == 1
    q2 = q2';
end

% ���������Ԫ��ֻ������Ԫ��������Ϊ�������Ԫ��
if length(q1) == 3
    q1 = [0, q1']';
end
if length(q2) == 3
    q2 = [0, q2']';
end
% *** ��ʼ�� ���� ***

MQ1 = [q1(1), -q1(2), -q1(3), -q1(4)
       q1(2),  q1(1), -q1(4),  q1(3)
       q1(3),  q1(4),  q1(1), -q1(2)
       q1(4), -q1(3),  q1(2),  q1(1)];

q = MQ1*q2;

end

