%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   SINS/GPS Intergration Navigation System test version 1.0
%   Indirect kalman filter������У������
%
% Reference:
%   - @�����ߵ��㷨����ϵ���ԭ��,�Ϲ��� P239
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1


close all
clear all
clc

gvar_earth;

% ���θ�����ʹ�õ�������
nn = 2;
% ����ʱ��
ts = 0.01;
nts = nn*ts;

% ��ʼ��̬��λ�á��ٶ�
att0 = [0, 0, 30*arcdeg]';
pos0 = [34*arcdeg, 108*arcdeg, 100]';
vn0 = [0, 0, 0]';
qbn0 = a2qua(att0);

% ��̬��Ԫ�����ٶȡ�λ��
qbn = qbn0;
vn = vn0;
pos = pos0;

eth = earth(pos, vn);

% *** ȫ��̬���棬����û���κ��˶���ת��
wm = qmulv(qconj(qbn), eth.wien)*ts;
vm = qmulv(qconj(qbn), -eth.gn)*ts;

wm = repmat(wm', nn, 1);
vm = repmat(vm', nn, 1);

% *** ������ *** 
% ʧ׼��
phi = [0.1, 0.2, 1]'*arcmin;
qbn = qaddphi(qbn, phi);

% ������ƫ���Ƕ��������
eb = [0.01, 0.015, 0.02]'*dph;
web = [0.001, 0.001, 0.001]'*dpsh;

% �Ӽ���ƫ���ٶ��������
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

len = fix(200/ts);
err = zeros(len, 10);
xkpk = zeros(len, 2*kf.n + 1);

kk = 1;
t = 0;
for k = 1 : nn : len
    t = t + nts;
    
    % ģ��imu���
    [wm1, vm1] = imuadderr(wm, vm, eb, web, db, wdb, ts);
    
    % �ߵ�����
    [qbn, vn, pos, eth] = insupdate(qbn, vn, pos, wm1, vm1, ts);
    
    % ״̬ת�ƾ���
    kf.Phikk_1 = eye(15) + kfft15(eth, q2mat(qbn), sum(vm1, 1)'/nts)*nts;
    
    kf = kfupdate(kf);
    
    % ģ��GPS���� 1Hz
    % mod: ������
    if mod(t, 0.2) < nts
        gps = [vn0; pos0] + rk.*randn(6, 1);
        kf = kfupdate(kf, [vn', pos']' - gps, 'M');
    end
    
    % ���� (����У��)???
    qbn = qdelphi(qbn, kf.Xk(1:3));
    kf.Xk(1:3) = 0;
    
    vn = vn - kf.Xk(4:6);
    kf.Xk(4:6) = 0;
    
    pos = pos - kf.Xk(7:9);
    kf.Xk(7:9) = 0;
    
    err(kk, :) = [qq2phi(qbn, qbn0)', (vn - vn0)', (pos - pos0)', t];    
    xkpk(kk, :) = [kf.Xk', diag(kf.Pk)', t]';
    
    kk = kk + 1;
    
    % ��������ʱ��ʾ��ǰ����
    if mod(t, 200) < nts
        disp(fix(t));
    end
end

% Ϊ����err���㹻�Ŀռ䣬�ڳ�ʼ��ʱ���ǽ��䳤������Ϊlen�����ڲ��ö�������
% �����߱��ĳЩԵ�ʣ�errͨ����װ���������ò�������Ϊ�˰Ѷ����0�õ���
err(kk:end, :) = [];
xkpk(kk:end, :) = [];
tt = err(:, end);
%% ��ͼ
msplot(321, tt, err(:, 1:2)/arcdeg, '\it\phi\rm/(\circ)');
legend('\it\phi\rm_E', '\it\phi\rm_N');

msplot(322, tt, err(:, 3)/arcdeg, '\it\phi\rm_U\rm/(\circ)');
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

msplot(321, tt, spk(:, 1:2)/arcmin, '\it\phi\rm/(\prime)');
legend('\it\phi\rm_E', '\it\phi\rm_N');

msplot(322, tt, spk(:, 3)/arcmin, '\it\phi\rm_U\rm/(\prime)');
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
