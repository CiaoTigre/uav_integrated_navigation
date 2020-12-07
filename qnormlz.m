function [ qbn ] = qnormlz( qbn )
% ���ƣ�Quaternion Normlization
% ���ܣ�����Ԫ����׼��
%
% Inputs:
%       qbn: δ��׼����Ԫ��
% Outputs:
%       qbn: ��׼���������Ԫ��

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ������Ԫ����ģ��
nm2 = qbn'*qbn;

if nm2 < 1e-6
    qbn = [1, 0, 0, 0]';
else
    qbn = qbn/sqrt(nm2);
end

end
