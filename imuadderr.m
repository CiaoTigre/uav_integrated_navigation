function [ wm, vm ] = imuadderr( wm, vm, eb, web, db, wdb, ts )
% ���ƣ�imu add error
% ���ܣ�����ʵimu��������ģ��ʵ���д�������������imu�����
% ����@ �����ߵ��㷨����ϵ���ԭ�� P235
%
% Inputs:
%       wm: n��3 ����ÿһ�ж�Ӧһ����������[theta_x, theta_y, theta_z] rad
%       vm: n��3 ����ÿһ�ж�Ӧһ���ٶ����� m/s
%       ע�⣡wm��vm��������У����һ���Ƕ�Ӧ��ʱ��
%       eb: �����ǳ�ֵƯ�ƣ���������ƫ
%       web:�Ƕ�����������
%       db: ���ٶȼƳ�ֵƫ��
%       wdb: �ٶ�����������
%       ts: ���� s
% Outputs:
%       wm: ��������Ľ�����
%       vm: ����������ٶ�����

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ��������
[m, n] = size(wm);
% ������ƽ��
sts = sqrt(ts);

% �Ƿ񱣴���ʱ���
switch n
    case 3
        % ����Ϊ3��û�б���ʱ��
        wm = wm + [ts*eb(1) + sts*web(1)*randn(m, 1), ...
                   ts*eb(2) + sts*web(2)*randn(m, 1), ...
                   ts*eb(3) + sts*web(3)*randn(m, 1)];
        
        vm = vm + [ts*db(1) + sts*wdb(1)*randn(m, 1), ...
                   ts*db(2) + sts*wdb(2)*randn(m, 1), ...
                   ts*db(3) + sts*wdb(3)*randn(m, 1)];
        
    case 4
        % ����Ϊ4�����һ�б���ʱ��
        wm(:, 1:3) = wm(:, 1:3) + [ts*eb(1) + sts*web(1)*randn(m, 1), ...
                                   ts*eb(2) + sts*web(2)*randn(m, 1), ...
                                   ts*eb(3) + sts*web(3)*randn(m, 1)];
        
        vm(:, 1:3) = vm(:, 1:3) + [ts*db(1) + sts*wdb(1)*randn(m, 1), ...
                                   ts*db(2) + sts*wdb(2)*randn(m, 1), ...
                                   ts*db(3) + sts*wdb(3)*randn(m, 1)];
    otherwise
        disp('imuadderr.m ���������ʽ����');
end

end
