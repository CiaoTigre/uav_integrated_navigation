%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   SINS/GPS Intergration Navigation System test version 2.0
%   Indirect kalman filter�����У������
%
% References:
%   - @�����ߵ��㷨����ϵ���ԭ��,�Ϲ��� P239
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1


close all
clear
clc

gvar_earth;


% loadģ��켣����: avp_SD and imu_SD
load('trajectory_simulator_data.mat');


% ���θ�����ʹ�õ�������
nn = 2;
% ����ʱ��
ts = 0.01;
nts = nn*ts;


% ��ʼ��̬���ٶȡ�λ��
att0 = [0, 0, 90]'*arcdeg;
vn0 = [0, 0, 0]';
pos0 = [34*arcdeg, 108*arcdeg, 100]'; % lattitude, longtitude, height
qbn0 = a2qua(att0);


% ��̬��Ԫ�����ٶȡ�λ��
qbn = qbn0;
vn = vn0;
pos = pos0;

eth = earth(pos, vn);


% *** ������ *** 
% ʧ׼��
phi = [0.1, 0.2, 1]'*arcmin;
qbn = qaddphi(qbn, phi);

% ������ƫ���Ƕ��������
eb_ref = [0.1, 0.15, 0.2]'*dph;
eb = [0.01, 0.015, 0.02]'*dph;
web = [0.001, 0.001, 0.001]'*dpsh;

% �Ӽ���ƫ���ٶ��������
db_ref = [800, 900, 1000]'*ug;
db = [80, 90, 100]'*ug;
wdb = [1, 1, 1]'*ugpsHz;

Qk = diag([web', wdb', zeros(1, 9)]')^2*nts;

rk = [[0.1, 0.1, 0.1], [10/Re, 10/Re, 10]]';
Rk = diag(rk)^2;

% Э�������x = [phi, delta_vn, delta_p, eb, db]
P0 = diag([[0.1, 0.1, 10]*arcdeg, [1, 1, 1], [10/Re, 10/Re, 10]...
           [0.1, 0.1, 0.1]*dph, [80, 90, 100]*ug]')^2;

% �������
Hk = [zeros(6,3), eye(6), zeros(6, 6)];

% Kalman filter initialization
kf = kfinit(Qk, Rk, P0, zeros(15), Hk);

% ��ģ��켣ʱ��һ��
kTime = fix(100/ts);   
err = zeros(kTime, 10);
xkpk = zeros(kTime, 2*kf.n + 1);

kk = 1;
t = 0;
for k = 2 : nn : kTime
    t = t + nts;
    
    % ��ȡģ��켣��Ӧ��imu���: ���������ٶ��������ο�ֵ��
    wm(1:nn,:) = imu_SD.wm(k-nn+1:k,:);
    vm(1:nn,:) = imu_SD.vm(k-nn+1:k,:);
    
    % Ϊimu�ο�����������
    [wm1, vm1] = imuadderr(wm, vm, eb, web, db, wdb, ts);
    
    % �ߵ����£���̬��Ԫ�����ٶȡ�λ�� 
    [qbn, vn, pos, eth] = insupdate(qbn, vn, pos, wm1, vm1, ts);
    
    % ����ģ��Ԥ�⣺�������ϵͳģ�Ϳ������˲�
    kf.Phikk_1 = eye(15) + kfft15(eth, q2mat(qbn), sum(vm1, 1)'/nts)*nts;
    kf = kfupdate(kf);
    
    % ������£�ģ��GPS�������� 5Hz
    if mod(t, 0.2) < nts
        gps = [avp_SD.vn(k,:)'; avp_SD.pos(k,:)'] + rk.*randn(6, 1);
        Zk = [vn', pos']' - gps;
        kf = kfupdate(kf, Zk, 'M');
    end
    
    % ����˲��������������ֻ���µ��������δ�Թߵ�ϵͳ����У����
    qbn_nav = qdelphi(qbn, kf.Xk(1:3));
    vn_nav  = vn - kf.Xk(4:6);
    pos_nav = pos - kf.Xk(7:9);

    % compute the err between the navigation output(estimation) & truth data
    qbn_ref = a2qua(avp_SD.att(k,:));
    vn_ref = avp_SD.vn(k,:)';
    pos_ref = avp_SD.pos(k,:)';
    err(kk, :) = [qq2phi(qbn_nav, qbn_ref)', (vn_nav - vn_ref)', (pos_nav - pos_ref)', t];
    xkpk(kk, :) = [kf.Xk', diag(kf.Pk)', t]';
    
    kk = kk + 1;
    
%     % ��������ʱ��ʾ��ǰ����
%     if mod(t, 50) == 0
% %         disp(fix(t));
%         disp('...');
%     end
end

% Ϊ����err���㹻�Ŀռ䣬�ڳ�ʼ��ʱ���ǽ��䳤������Ϊlen�����ڲ��ö�������
% �����߱��ĳЩԵ�ʣ�errͨ����װ���������ò�������Ϊ�˰Ѷ����0�õ���
err(kk:end, :) = [];
xkpk(kk:end, :) = [];
tt = err(:, end);


%% ��ͼ
msplot(321, tt, err(:, 1:2)/arcdeg, '\it\phi\rm/(\circ)');  % ��
legend('\it\phi\rm_E', '\it\phi\rm_N');

msplot(322, tt, err(:, 3)/arcdeg, '\it\phi\rm_U\rm/(\circ)'); % ��
legend('\it\phi\rm_U');

msplot(323, tt, err(:, 4:6), '\delta\itv^n\rm/(m.s^{-1})');
legend('\delta\itv\rm_E', '\delta\itv\rm_N', '\delta\itv\rm_U');

msplot(324, tt, [err(:, 7)*Re, err(:, 8)*Re*cos(pos(1)), err(:, 9)],...
    '\delta\itp\rm/m');
legend('\delta\itL', '\delta\it\lambda', '\delta\ith');

msplot(325, tt, xkpk(:, 10:12)/dph, '\it\epsilon\rm/(\circ.h^{-1})');
legend('\it\epsilon_x', '\it\epsilon_y', '\it\epsilon_z');

msplot(326, tt, xkpk(:, 13:15)/ug, '\it\nabla\rm/\mu\itg');
legend('\it\nabla_x', '\it\nabla_y', '\it\nabla_z');

% �����������ͼ
spk = sqrt(xkpk(:, 16:end-1 ));

msplot(321, tt, spk(:, 1:2)/arcdeg, '\it\phi\rm/(\circ)');
legend('\it\phi\rm_E', '\it\phi\rm_N');

msplot(322, tt, spk(:, 3)/arcdeg, '\it\phi\rm_U\rm/(\circ)');
legend('\it\phi\rm_U');

msplot(323, tt, spk(:, 4:6), '\delta\itv^n\rm/(m.s^{-1})');
legend('\delta\itv\rm_E', '\delta\itv\rm_N', '\delta\itv\rm_U');

msplot(324, tt, [spk(:, 7)*Re, spk(:, 8)*Re*pos(1), spk(:, 9)],...
       '\delta\itp\rm/m');
legend('\delta\itL', '\delta\it\lambda', '\delta\ith');

msplot(325, tt, spk(:, 10:12)/dph, '\it\epsilon\rm/(\circ.h^{-1})');
legend('\it\epsilon_x', '\it\epsilon_y', '\it\epsilon_z');

msplot(326, tt, spk(:, 13:15)/ug, '\it\nabla\rm/\mu\itg');
legend('\it\nabla_x', '\it\nabla_y', '\it\nabla_z');
