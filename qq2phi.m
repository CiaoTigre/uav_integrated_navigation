function [ phi ] = qq2phi( qbp, qbn )
% ���ƣ�Quaternion Quaternion to PHI
% ���ܣ����㺬��������̬��Ԫ������ֵ֮���Ӧ�ĵ�Ч��תʸ��
%
% Inputs:
%       qbp: ����������̬��Ԫ��
%       qbn: ��ʵ��̬��Ԫ��
% Outputs:
%       phi: nϵ��pϵ�ĵ�Ч��תʸ��

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
qerr = qmul(qbn, qconj(qbp));
phi = q2rv(qerr);

end

