function [ Cbn ] = q2mat( qbn )
% ���ƣ�Quaternion to directin cosine MATrix
% ���ܣ�����̬��Ԫ�����㷽�����Ҿ���
%
% Inputs:
%       qbn: rn = qbn*rb*qbn'
% Outputs:
%       Cbn: ��bϵ��nϵ������ת������

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ����û�ò��ý���Ԫ���ֽ�Ϊ��������+ʸ�����ֵļ��㷽ʽ�����ǲ���ֱ������Ԫ��
% �и���Ԫ��������Cbn��ÿһ��Ԫ�صļ��㷽ʽ��
% ����@ �����ߵ��㷨����ϵ������� P232
q11 = qbn(1)*qbn(1);    q12 = qbn(1)*qbn(2);    q13 = qbn(1)*qbn(3);    q14 = qbn(1)*qbn(4);
q22 = qbn(2)*qbn(2);    q23 = qbn(2)*qbn(3);    q24 = qbn(2)*qbn(4);
q33 = qbn(3)*qbn(3);    q34 = qbn(3)*qbn(4);
q44 = qbn(4)*qbn(4);
    
Cbn = [q11 + q22 - q33 - q44,           2*(q23 - q14),            2*(q24 + q13)
               2*(q23 + q14),   q11 - q22 + q33 - q44,            2*(q34 - q12)
               2*(q24 - q13),           2*(q34 + q12),   q11 - q22 - q33 + q44];

end

