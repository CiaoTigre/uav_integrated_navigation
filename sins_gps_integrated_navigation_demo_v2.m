%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   SINS/GPS Intergration Navigation System test version 2.0
%   Indirect kalman filter������У������
%
% Reference:
%   - @�����ߵ��㷨����ϵ���ԭ��,�Ϲ��� P239
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1


close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Uav true trajectory data.
addpath UavTrajectorySim;

disp(' ')
disp('Available UAV Truth Trajectory Data Files:')
dir_mat_files = dir('UavTrajectorySim\*.mat');
for nFile=1:length(dir_mat_files)
    fprintf('   %d: %s\n',nFile,dir_mat_files(nFile).name);
end
% nFileChoice = input('Choose a UAV Truth data file (e.g. 1<Enter>): ');
try
%     load(dir_mat_files(nFileChoice).name)
    load(dir_mat_files(1).name)
catch
    error('Selected UAV Truth Trajectory data file (%d) is invalid.\n',nFileChoice);
end



gvar_earth;

% ���θ�����ʹ�õ�������
nn = 2;
% ����ʱ��
ts = 0.01;
nts = nn*ts;


% ��ʼ��̬���ٶȡ�λ��
att0 = [0, 0, 90]'*arcdeg;
vn0  = [0, 0, 0]';
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

rk = [[0.1, 0.1, 0.1], [5/Re, 5/Re, 5]]';
Rk = diag(rk)^2;

% Э�������x = [phi, delta_vn, delta_p, eb, db]
P0 = diag([[0.1, 0.1, 10]*arcdeg, [1, 1, 1], [10/Re, 10/Re, 10]...
           [0.1, 0.1, 0.1]*dph, [80, 90, 100]*ug]')^2;

% �������
Hk = [zeros(6,3), eye(6), zeros(6, 6)];

% Kalman filter initialization
kf = kfinit(Qk, Rk, P0, zeros(15), Hk);

% ��ģ��켣ʱ��һ��
kTime = fix(t_SD/ts);   
err = zeros(kTime, 10);
xkpk = zeros(kTime, 2*kf.n + 1);

pos_ref = zeros(kTime,3);
pos_est = zeros(kTime,3);
pos_gps = zeros(kTime,3);

kk = 1;
t = 0;
for k = 2 : nn : kTime
    t = t + nts;
    
    % ��ȡģ��켣��Ӧ��imu���: ���������ٶ��������ο�ֵ��
    wm(1:nn,:) = imu_SD.wm(k-nn+1:k,:);
    vm(1:nn,:) = imu_SD.vm(k-nn+1:k,:);
    
    % ΪIMU�ο�����������
    [wm1, vm1] = imuadderr(wm, vm, eb, web, db, wdb, ts);
    
    % �ߵ����£���̬��Ԫ�����ٶȡ�λ�� 
    [qbn, vn, pos, eth] = insupdate(qbn, vn, pos, wm1, vm1, ts);
    
    % ����ģ��Ԥ�⣺�������ϵͳģ�Ϳ������˲�
    kf.Phikk_1 = eye(15) + kfft15(eth, q2mat(qbn), sum(vm1, 1)'/nts)*nts;
    kf = kfupdate(kf);
    
    % ģ��GPS��������
    gps = [avp_SD.vn(k,:)'; avp_SD.pos(k,:)'] + rk.*randn(6, 1);
    pos_gps(kk,:) = gps(4:6)';
    % ������� 5Hz
    if mod(t, 0.2) < nts
        Zk = [vn', pos']' - gps;
        kf = kfupdate(kf, Zk, 'M');
    end
    
    % Indirect Kalman filter��feedback to IMU (����У����)
    qbn = qdelphi(qbn, kf.Xk(1:3));
    vn  = vn - kf.Xk(4:6);
    pos = pos - kf.Xk(7:9);
    pos_est(kk,:) = pos';
    % ����У�������ڷ�����Ĵ��ڵ��¿������˲����������ֵʼ��Ϊ��. Ref: ������
    kf.Xk(1:3) = 0;
    kf.Xk(4:6) = 0;
    kf.Xk(7:9) = 0;
%     kf.Xk(10:12) = 0;
%     kf.Xk(13:15) = 0;
    
        
    % compute the error between estimation & truth data 
    % Note that this 'error' is not the 'state vector' in the Kalman equ. 
    % In indirect kalman filter, the 'state vector' means the error of 
    % the IMU update (respect to True data.)
    qbn_ref = a2qua(avp_SD.att(k,:));
    vn_ref = avp_SD.vn(k,:)';
    pos_ref(kk,:) = avp_SD.pos(k,:);
    err(kk, :) = [qq2phi(qbn, qbn_ref)', (vn - vn_ref)', (pos - pos_ref(kk,:)')', t];
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
pos_ref(kk:end,:) = [];
pos_est(kk:end,:) = [];
pos_gps(kk:end,:) = [];
tt = err(:, end);

%% �����ǻ�ͼ����
figure;
subplot(3,3,[1,4]);  
% �����Ǿ���Lontitu��������γ��Latitude
plot(pos_gps(:,2)/arcdeg,pos_gps(:,1)/arcdeg, 'dg','LineWidth',0.1); hold on;
plot(pos_est(:,2)/arcdeg,pos_est(:,1)/arcdeg, 'r','LineWidth',4); hold on;
plot(pos_ref(:,2)/arcdeg,pos_ref(:,1)/arcdeg, 'b','LineWidth',1); hold on;
plot(pos_ref(1,2)/arcdeg,pos_ref(1,1)/arcdeg, 'oc','LineWidth',4);
% axis equal;
grid on;
xlabel('\it\lambda\rm /(\circ)');
ylabel('\itL\rm /(\circ)');
legend('GPS meas.','Est. pos.','True pos.', 'Start');
title('UAV Position')

subplot(3,3,2);
plot(tt, err(:, 1:2)/arcdeg);
grid on;
axis tight;
xlabel('t/s');
ylabel('\it\phi\rm/(\circ)');
legend('\it\phi\rm_E', '\it\phi\rm_N');
title('Pitch & Roll Est. error')

subplot(3,3,3);
plot(tt, err(:, 3)/arcdeg);
grid on;
axis tight;
% ylim([-10,10])
xlabel('t/s');
ylabel('\it\phi\rm_U\rm/(\circ)');
legend('\it\phi\rm_U');
title('Yaw Est. error')

subplot(3,3,5);
plot(tt, err(:, 4:6));
grid on;
axis tight;
xlabel('t/s');
ylabel('\delta\itv^n\rm/(m.s^{-1})');
legend('\delta\itv\rm_E', '\delta\itv\rm_N', '\delta\itv\rm_U');
title('Velocity Est. error')

subplot(3,3,6);
plot(tt, [err(:, 7)*Re, err(:, 8)*Re*cos(pos(1)), err(:, 9)]);
grid on;
axis tight;
ylim([-10,10]); 
xlabel('t/s');
ylabel('\delta\itp\rm/m');
legend('\delta\itL', '\delta\it\lambda', '\delta\ith');
title('Position Est. error')


subplot(3,3,7);  
plot(tt,pos_gps(:,3), ':g','LineWidth',0.1); hold on;
plot(tt,pos_est(:,3), 'r','LineWidth',2); hold on;
plot(tt,pos_ref(:,3), 'b','LineWidth',1); 
grid on;
axis tight;
xlabel('t/s');
ylabel('\ith\rm /(m)');
legend('GPS meas.','Est. Alt.','True Alt.');
title('UAV Altitude')

subplot(3,3,8);
plot(tt, xkpk(:, 10:12)/dph);
grid on;
axis tight;
xlabel('t/s');
ylabel('\it\epsilon\rm/(\circ.h^{-1})');
legend('\it\epsilon_x', '\it\epsilon_y', '\it\epsilon_z');
title('Gyro biases')

subplot(3,3,9);
plot(tt, xkpk(:, 13:15)/ug);
grid on;
axis tight;
xlabel('t/s');
ylabel('\it\nabla\rm/\mu\itg');
legend('\it\nabla_x', '\it\nabla_y', '\it\nabla_z');
title('Accelerometer biases')

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

% ��ά�켣
figure(3)
plot3(pos_gps(:,2)/arcdeg, pos_gps(:,1)/arcdeg, pos_gps(:,3), ':g','LineWidth',0.1); hold on;
plot3(pos_ref(:,2)/arcdeg, pos_ref(:,1)/arcdeg, pos_ref(:,3),'b','LineWidth',2); hold on;
plot3(pos_est(:,2)/arcdeg, pos_est(:,1)/arcdeg, pos_est(:,3),'r','LineWidth',3); hold on;
plot3(pos_ref(1,2)/arcdeg, pos_ref(1,1)/arcdeg, pos_ref(1,3),'oc','LineWidth',10);
xlabel('\it\lambda\rm /(\circ)');
ylabel('\itL\rm /(\circ)');
zlabel('\ith\rm /(m)')
legend('GPS meas.','True pos.','Est. pos.','start')
grid on;
title('3D Trajectory')
