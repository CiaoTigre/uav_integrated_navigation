function [ qbp ] = qaddphi( qbn, phi )
% ���ƣ�Quaternion ADD PHI
% ���ܣ���Ԫ����ʧ׼�����,��������ʵ��̬��Ԫ�����������ģ�⺬��������̬
%      ��Ԫ����
%
% ��Ϊ����n'����ʾ����ĵ�������ϵ���ڹ��ɾ���ı�ʾ�����в������壬�����p��
% ����n'ϵ��
% Cbp = Cnp*Cbn, ͬ���У�qbp = qnp*qbn
% Inputs:
%       qbn: ��ʵ��̬��Ԫ��
%       phi: ��nϵ��pϵ�ĵ�Ч��תʸ��
% Outputs:
%       pbp: ����������̬��Ԫ��

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ����@ �����ߵ��㷨����ϵ���ԭ���� P234
qbp = qmul(rv2q(-phi), qbn);

% ��һ�ּ��㷽����
% qpn = rv2q(phi);
% qnp = qconj(qpn);
% 
% qbp = qmul(qnp, qbn);

end

