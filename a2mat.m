function [ Cbn ] = a2mat( att )
% Function��Attitude to dirction cosine MATrix 
% ��ȫ�����Ϲ�����ʦ���ϸ��Ƶİ汾
%   
% Inputs:
%       att:	��̬��(pitch-��,roll-��,yaw-��), unit:rad
%
% Outputs:
%       Cbn:    ��bϵ��nϵ������任����, Ҳ�Ǵ�nϵ��bϵ������ϵת������
%
% Reference:
%   - @ �����ߵ��㷨����ϵ���ԭ�� P231
%
% Notes:
%   - ��������ʱ�ϽǱ���ǰ
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1

s = sin(att);
c = cos(att);

si = s(1);
sj = s(2);
sk = s(3);

ci = c(1);
cj = c(2);
ck = c(3);

Cbn = [cj*ck - si*sj*sk, -ci*sk, sj*ck + si*cj*sk
       cj*sk + si*sj*ck,  ci*ck, sj*sk - si*cj*ck
                 -ci*sj,     si,            ci*cj];
end
