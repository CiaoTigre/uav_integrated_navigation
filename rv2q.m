function [ qbn ] = rv2q( rv )
% ���ƣ�Rotation Vector to Quaternion
% ���ܣ�����Ч��תʸ��ת��Ϊ��̬��Ԫ��
%
% Inputs:
%        rv: 3��1 vector, ��ģ��ʾת���ǶȵĴ�С��rv/norm(rv) Ϊת����������
%            ����ϵ�µ����ꡣ 
% Outputs:
%       qbn: rn = qbn*rb*qbn'

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ����@ �����ߵ��㷨����ϵ���ԭ���� P233
% ��Ч��תʸ����ģ����
nm2 = rv'*rv ;

if nm2 < 1e-8       % ���ģ����С��������̩�ռ���չ����ǰ�������Ǻ���
    q0 = 1 - nm2*(1/8 - nm2/384);
    s = 1/2 - nm2*(1/48 - nm2/3840);
else
    nm = sqrt(nm2);
    q0 = cos(nm/2);
    s = sin(nm/2)/nm;
end

qbn = [q0, s*rv']';

end

