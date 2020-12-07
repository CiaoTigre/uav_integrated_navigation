function [ att ] = q2att( qbn )
% ���ƣ�Attitude to QUAternion
% ���ܣ�����̬��Ԫ��ת��Ϊ��̬��
%
% Inputs:
%       qbn: rn = qbn*rb*qbn'
% Outputs:
%       att: ����������������Ϊ����theta�����gamma�ͺ���psi�� ����Ǳ�ƫ
%       ��Ϊ��,ȡֵΪ��-pi, pi]

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
att = m2att(q2mat(qbn));

end

