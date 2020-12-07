function [ dpos ] = deltapos( pos )
% ���ƣ�delta position
% ���ܣ���λ������
% ����@ �����ߵ�ϵͳ����ϵ���ԭ�� P236

% Inputs:
%       pos: N��3 ����ÿһ�б���һ��ʱ�̵�λ��
% Output:
%       dpos: N��3 ����ÿһ�ж�Ӧ��pos�ж�Ӧ��Ԫ�����һ��λ��pos(, :)֮��

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% cos(latitude)
cl = cos(pos(:, 1));
% ����뾶 m���˴������ǵ���ƽ���뾶
Re = 6378137;

dpos = [(pos(:, 1) - pos(1, 1))*Re, (pos(:, 2) - pos(1, 2)).*cl*Re,...
        pos(:, 3) - pos(1, 3) ];

end
