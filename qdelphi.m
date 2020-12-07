function [ qbn ] = qdelphi( qbp, phi )
% ���ƣ�Quaternion DELete PHI
% ���ܣ���Ԫ����ʧ׼���������Ӻ���������̬��Ԫ�����޳����Ի����ֵ��
%
% ��Ϊ����n'����ʾ���㵼������ϵ���ڹ��ɾ���ı�ʾ�����в������壬�����p��
% ����n'ϵ��
% Inputs:
%       qbp: ��bϵ�����㵼��ϵp����̬��Ԫ����
%       phi: ��nϵ��pϵ���ת������Ӧ�ĵ�Ч��תʸ��
%       @ �����ߵ��㷨����ϵ���ԭ�� P234 ԭ�ģ�����ʵ����ϵn�����㵼��ϵp��
%         ��ʧ׼��Ϊphi����֮����pϵ��nϵ��ʧ׼��Ϊ-phi������-phi��Ϊ��Ч��
%         תʸ�����������Ӧ�ı任��Ԫ��ΪQnp��
%         �����˵�����н�ʧ׼��(312 ŷ����)��Ϊ��Ч��תʸ����
% Outputs:
%       qbn: ��bϵ����ʵ����ϵ����̬��Ԫ��

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
qbn = qmul(rv2q(phi), qbp);

end
