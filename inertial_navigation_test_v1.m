%% **************************************************************
%���� inertial navigation test version 1.0
%���ܣ����ߵ��������
%����@ �����ߵ�ϵͳ����ϵ���ԭ�� P236
%
%

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
close all
clear all
gvar_earth;

% �ߵ�����ʱʹ�õĶ�������Ŀ�Ͳ���ʱ��
nn = 2;
ts = 0.1;
nts = nn*ts;

% ��̬���ٶȡ�λ�õĳ�ʼֵ
att = [1, 1, 30]'*arcdeg;
vn = [0, 0, 0]';
pos = [34*arcdeg, 108*arcdeg, 100]';
qbn = a2qua(att);

eth = earth(pos, vn);

% ����ģ������
wm = qmulv(qconj(qbn), eth.wien)*ts;
vm = qmulv(qconj(qbn), -eth.gn)*ts;
% ���澲̬imu����
wm = repmat(wm', nn, 1);
vm = repmat(vm', nn, 1);

% ʧ׼��
phi = [0.1, 0.2, 3]'*arcmin;
% ����������̬��Ԫ��
qbn = qaddphi(qbn, phi);

% ����ʱ��
len = fix(3600/ts);

% ���ߵ������õ���̬���ٶȡ�λ����Ϣ������Ϊ��Ԥ����洢�ռ�
avp = zeros(len, 10);

kk = 1;
t = 0;
%% pure inertial navigation
for k = 1 : nn : len
    % ��ǰʱ��
    t = t + nts;
    
    % ���ߵ�����
    [qbn, vn, pos] = insupdate(qbn, vn, pos, wm, vm, ts);
    
    % ���������
    avp(kk, :) = [q2att(qbn)', vn', pos', t];
    
    kk = kk + 1;
    
    % �ڳ�������ʱ��ʾ��ǰ����
    if mod(t, 100) < nts
        disp(fix(t));
    end
    
end
%% ��ͼ
avp(kk:end, :) = [];

tt = avp(:, end);

msplot(221, tt, avp(:, 1:2)/arcdeg, 'Att/(\circ)');
legend('\it\theta', '\it\gamma');

msplot(222, tt, avp(:, 3)/arcdeg, '\psi /\circ');

msplot(223, tt ,avp(:, 4:6), 'Vel /m.s^{-1}');
legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U');

msplot(224, tt, deltapos(avp(:, 7:9)), '\DeltaPos /m');
legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith');
