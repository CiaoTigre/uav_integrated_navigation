function [ Ng ] = solar_disturbing_force( Lbe )
% ���ƣ�solar disturbing force
% ���ܣ�̫���㶯��
%
% Inputs:
%       Lbe: �������ĵ������ļ����        m
%       where L denotes distance, b denote body, e denote earth.
% Outputs:
%       Ng: �㶯�����ٶ�����        g(9.8m/s^2)
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1

%%
% ������������ (N*m^2/kg^2)
G = 6.67259e-11;
% ����뾶 (m)
Re = 6.371e6;
% ����ת����뾶 (m)
Lse = 1.496e11;
% ̫������ (kg)
Ms = 2e30;

if nargin == 0
    delta_r = Re;
else
    delta_r = Lbe;
end

% disturbing acceleration
g_dst = 2*G*Ms*delta_r/Lse^3;
Ng = g_dst/9.8;

end

