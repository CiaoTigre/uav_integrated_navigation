function [ qout ] = qconj( qin )
% ���ƣ�CONJugate Quaternion
% ���ܣ����������Ԫ��
%
% Inputs:
%       qin: ������Ԫ��
% Outputs:
%       qout: ����������Ĺ���

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%% 
qout = zeros(4, 1);

qout(1) = qin(1);
qout(2:4) = -qin(2:4);

end

