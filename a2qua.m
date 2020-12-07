function [ qbn ] = a2qua( att )
% Function��Attitude to Quaternion
%   
% Inputs:
%       att:	��̬��(pitch-��,roll-��,yaw-��), unit:rad
%
% Outputs:
%       qbn:    ��bϵ��nϵ����Ԫ��
%
% Reference:
%   - @ �����ߵ��㷨����ϵ���ԭ�� P232  @ԭ�� P247
%
% Notes:
%   - ��������ʱ�ϽǱ���ǰ
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1

s = sin(att/2);
c = cos(att/2);

si = s(1);  sj = s(2);  sk = s(3);
ci = c(1);  cj = c(2);  ck = c(3);

qbn = [ci*cj*ck - si*sj*sk
       si*cj*ck - ci*sj*sk
       ci*sj*ck + si*cj*sk
       ci*cj*sk + si*sj*ck];

end

