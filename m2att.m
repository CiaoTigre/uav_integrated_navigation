function [ att ] = m2att( Cbn )
% ���ƣ�direction cosine Matrix to ATTitude
% ���ܣ���̬����ת������̬��
%
% Inputs:
%       Cbn: ��bϵ��nϵ������ת������
% Outputs:
%       att: ����������������Ϊ����theta�����gamma�ͺ���psi�� ����Ǳ�ƫ
%       ��Ϊ��,ȡֵΪ��-pi, pi]

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%% 
% ԭ��@ �����ߵ��㷨����ϵ���ԭ�� P245
% ��theta �� ��pi/2��ʱ��gamma��psi���޷�����ģ���֪����
%   gamma + psi = atan2(Cbn(1, 3), Cbn(1, 1));  (theta �� +pi/2)
%   gamma - psi = atan2(Cbn(1, 3), Cbn(1, 1));  (theta �� -pi/2)
% gamma��psi���ڶ�ֵ�ԣ�ֻ��ָ��������һ����ֵ����ȷ����һ����������psi = 0

% ע��: ��������൱�ڸ����ǵ��ڡ�90����������µ����ϵ����岻��������
% ���˶��켣�����������˶��������������ϴ��ڣ��������ڽ��е�������ʱ�ڼ���
% �в�����

if abs(Cbn(3,2)) <= 0.999999
    att = [asin(Cbn(3, 2)), -atan2(Cbn(3, 1), Cbn(3, 3)), -atan2(Cbn(1, 2), Cbn(2, 2))]';
else
    att = [asin(Cbn(3, 2)), atan2(Cbn(1, 3), Cbn(1, 1)), 0]';
end

end