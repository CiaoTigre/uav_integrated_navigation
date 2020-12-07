function [ rv ] = q2rv( q )
% ���ƣ�Quaternion to Rotation Vector
% ���ܣ�����̬��Ԫ��ת��Ϊ��Ч��תʸ��
%
% Inputs:
%       q: 
% Outputs:
%       rv: 

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ԭ��@ �����ߵ��㷨����ϵ���ԭ���� P233
% ���Ƚ���Ԫ��ת��Ϊ�����Ǹ���Ԫ��
if q(1) < 0
    q = -q;
end

% ��Ч��תʸ��ģֵ��һ��
nmhalf = acos(q(1));

if nmhalf > 1e-20
    b = 2*nmhalf/sin(nmhalf);
else
    b = 2;
end

rv = b*q(2:4);

end
