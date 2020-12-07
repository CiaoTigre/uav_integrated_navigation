%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���ƣ�Opitiaml Based Alignment (GPS velocity aided) version 2.0
% ���ܣ�
%
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
close all
clear all

% ȫ�ֱ���
gvar_earth;

% ��������
% ����avp_SD.(att, pos, vn) ��imu_SD.(acc, gry)
load('trajectory_simulator_data.mat');

% ���ݲ���
ts_imu = 0.01;
% ���ö������㷨ʱʹ�õ�������
num_subsample = 2;
nts = num_subsample*ts_imu;

% ��msr=ref+errģ��imu���
imu_err = imuerrorset('selfdefine');
[imu_msr.gyr, imu_msr.acc] = imuadderr(imu_SD.wb, imu_SD.fb, ...
              imu_err.eb, imu_err.web, imu_err.db, imu_err.wdb, ts_imu);

% ʱ�����
init_time = 0;
k_init = fix(init_time/ts_imu) + 1;
total_time = 100;

% ״̬������ֵ
pos_init_ref = avp_SD.pos(k_init, :)';
vn_init_ref = avp_SD.vn(k_init, :)';
att_init_ref = avp_SD.att(k_init, :)';
Cbn0_ref = a2mat(att_init_ref);

% ������ǰһʱ��ֵ,�ڶ�����δ�����ܻ����õĳ�ֵ���
pos_prv = pos_init_ref + [0, 0, 0]';
vn_prv = vn_init_ref + [0, 0, 0]';
vn_init = vn_init_ref + [0, 0, 0]';
att_prv = att_init_ref + [5, 5, 100]'*arcdeg;
eth_prv = earth(pos_prv, vn_prv);
K_prv = zeros(4, 4);

% ʱ����̬�����ֵ
Cntn0_prv = eye(3);
Cbtb0_prv = eye(3);
alpha_sigma = [0, 0, 0]';
beta_prv = [0, 0, 0]';
beta = [0, 0, 0]';
qbtb0 = [1, 0, 0, 0]';
% Ϊ��̬���������ڴ�
att_cmpt_sv = zeros(length(avp_SD.att), 3);
att0_cmpt_sv = zeros(length(avp_SD.att), 3);
att_err_sv = zeros(length(avp_SD.att), 3);

% ����ѭ����һЩ��Ҫ�ļ�������
% ���imu������Ϣ�Ĵ���
Lp_msr = 0;
% ���и��µĴ���
Lp_update = 0;
% ��ǰʱ��(�Դӵ�����ʼʱ�̶���������ʼʱ��)
crt = ts_imu;
%%
% ���˼·�������ģ������У�ÿ��ts_imuʱ��imu�ͻ����һ����������ٶ�����
% ��Ϣ���ڲ��ö������㷨(�Զ�����Ϊ��)ʱ��ÿ���2��imu����ֵ���������������
% һ�ιߵ����¡�Lp�ĺ������imu�����Ϣ�Ĵ�����
while (crt + init_time) < total_time
    % ����ѭ����������+1
    Lp_msr = Lp_msr+1;
    
    % ��ǰʱ��������avp_SD��imu_SD�ж�Ӧ�ı�� 
    k = round(crt/ts_imu) + k_init;
    
    % ÿ����㹻������������ͽ���һ�ιߵ�����,����ʲôҲ����
    if mod(Lp_msr, num_subsample) == 0
        % ���йߵ����µĴ���
        Lp_update = Lp_update+1;
        
        % ����imu������Ϣ
        % ʵ������µ����������ÿ��⵽һ��imu����ͻᱣ��һ�Ρ�����Ϊ
        % ��̷��㣬�����Ϊ��ÿ��nts��������ƾ�һ���Ӷ���num_subsample
        % �β�����������һ�ζ��벢����һ������
        wm = imu_msr.gyr(k-num_subsample+1 : k, :);
        vm = imu_msr.acc(k-num_subsample+1 : k, :);
        
        % ����Բ׶/�������
        [phim, dvbm] = cnscl(wm, vm);
        
        % Ŀǰ���˼�����̬֮��ó����ݲ����йߵ����£�pos��vn����GPS����
        % ģ��GPS����(��ʱ�������)
        gps.pos = avp_SD.pos(k, :)';
        gps.vn = avp_SD.vn(k, :)';
        pos = gps.pos;
        vn = gps.vn; 
        
        % ����Cbtb0��Cntn0
        Cbtb0 = Cbtb0_prv*rv2m(phim);
        Cntn0 = Cntn0_prv*rv2m(eth_prv.winn*nts);
        
        % *** ����alpha(n0)��beta(b0) ***
        % alpha(n0)
        % 1. Ŀǰ��vn0_ref����alpha, ��ʱ�����ǻ������ڡ�
        % 2. ��gps�����λ�ü���wien��gn
        alpha_sigma = alpha_sigma + ...
            Cntn0_prv*(cross(eth_prv.wien, vn_prv) - eth_prv.gn)*nts;
        alpha = Cntn0*vn - vn_init + alpha_sigma;
        % beta(b0)
        beta = beta_prv + Cbtb0_prv*dvbm;
        
        % QUEST����qbn0
        [ qbn0, K_prv ] = QUEST( beta, alpha, K_prv );
        
        % 
        Cbn0 = q2mat(qbn0); 
        Cbn = Cntn0'*Cbn0*Cbtb0;
        att_cmpt_sv(Lp_update, :) = m2att(Cbn)';
        att_err_sv(Lp_update, :) = atterrnorml(m2att(Cbn)' ...
                                        - avp_SD.att(k, :));
        att0_cmpt_sv(Lp_update, :) = q2att(qbn0)';
        
        % ����һЩ��Ҫ���ٶ�
        eth = earth(pos, vn);
        
        % �����´�ѭ���ĳ�ֵ
        Cbtb0_prv = Cbtb0;
        Cntn0_prv = Cntn0;
        beta_prv = beta;
        eth_prv = eth;
        pos_prv = pos;
        vn_prv = vn;
    else
        % do noting
    end
        
    % �´ν���ѭ����ʱ��
    crt = crt + ts_imu;
    
end

% ɾ������ռ�
att_cmpt_sv(Lp_update+1 : end, :) = [];
att0_cmpt_sv(Lp_update+1 : end, :) = [];
att_err_sv(Lp_update+1 : end, :) = [];
%% ��ͼ
Time_axis = (1:1:Lp_update)*nts;
length_imu = (1:1:length(imu_SD.fb))*ts_imu;

% ������̬(Cbn0)
msplot(311, Time_axis, att0_cmpt_sv(:, 1)*deg, 'ʱ�� /s', 'pitch');
title('Cbn0');
msplot(312, Time_axis, att0_cmpt_sv(:, 2)*deg, 'ʱ�� /s', 'roll');
msplot(313, Time_axis, att0_cmpt_sv(:, 3)*deg, 'ʱ�� /s', 'yaw');

% ������̬(Cbn)
msplot(311, Time_axis, att_cmpt_sv(:, 1)*deg, 'ʱ�� /s', 'pitch');
title('change of Cbn');
msplot(312, Time_axis, att_cmpt_sv(:, 2)*deg, 'ʱ�� /s', 'roll');
msplot(313, Time_axis, att_cmpt_sv(:, 3)*deg, 'ʱ�� /s', 'yaw');

% ��̬���
msplot(311, Time_axis, att_err_sv(:, 1)*deg, 'ʱ�� /s', 'pitch');
title('attitude error');
msplot(312, Time_axis, att_err_sv(:, 2)*deg, 'ʱ�� /s', 'roll');
msplot(313, Time_axis, att_err_sv(:, 3)*deg, 'ʱ�� /s', 'yaw');

% imu���ģ��ֵ
msplot(311, length_imu, imu_SD.fb(:, 1), 'ʱ�� /s', 'fb(1)');
title('accelerometer inpute');
msplot(312, length_imu, imu_SD.fb(:, 2), 'ʱ�� /s', 'fb(2)');
msplot(313, length_imu, imu_SD.fb(:, 3), 'ʱ�� /s', 'fb(3)');

% imu���ģ��ֵ
msplot(311, length_imu, imu_SD.wb(:, 1), 'ʱ�� /s', 'wb(1)');
title('gyroscope inpute');
msplot(312, length_imu, imu_SD.wb(:, 2), 'ʱ�� /s', 'wb(2)');
msplot(313, length_imu, imu_SD.wb(:, 3), 'ʱ�� /s', 'wb(3)');
