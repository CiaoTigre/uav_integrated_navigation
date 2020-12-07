%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���ƣ�attitude update test
% Function: ��̬�����㷨���Գ���
% 
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
clear all
close all
% ����ͱ���һЩȫ�ֱ���
gvar_earth;
pos = [46*arcdeg, 126*arcdeg, 100]';
eth = earth(pos, [0, 0, 0]');

% ��������
load('angle_motion_data');
% ���ݳ���
[data_roll, data_column] = size(att_SD);
% ��ʱ��
total_time_of_data = att_SD(end, 4);
% �������㷨���õ�������
% ע�⣡����ʱ���֣�����"����"һ����Ϊ���Σ���õ���round(��������)������
% fix(�������)���������ʱ��������ݶ�Ӧ�쳣�����⡣
num_sample = round(2);
% imu�����Ӧ����
ts_imu = 0.01;

% ��ʼʱ��,e.g.����init_time=20s,���ö������㷨������׸���̬����t=20.02s
% ����̬
init_time = 0;
% t=init_time��������att_SD�е�index
k_init = round(init_time/ts_imu) + 1;

% ��ʼ��̬
Cbn0_ref = a2mat(att_SD(k_init, 1:3));

Cbn0 = Cbn0_ref;
Cntn0 = eye(3);
Cbtb0 = eye(3);

% Ϊ��̬�����ֱ�洢�ռ�
att_cmpt_sv1 = zeros(data_roll, 3);
att_ref_sv = zeros(data_roll, 3);
att_err_sv = zeros(data_roll, 3);
idx_sv = zeros(data_roll, 3);
% ѭ������
Lp = 0;
Lp_subsample = 0;

% ��ǰʱ��
crt = ts_imu;
i = 0;
%% 
while round((crt + init_time)/ts_imu) <= round(120/ts_imu)
    % ��Lp��ѭ��
    Lp = Lp+1;
    
    % crtʱ����������Ӧ�ı��
    i = round(crt/ts_imu) + k_init;
    idx_sv(Lp, :) = [crt + init_time, i, round((crt + init_time)/ts_imu) - i];
    
    % �ò�ͬ����������̬ 1.��������������Ԫ�� 2.�Ľ����������
    % *** ������ ***
    if mod(round(Lp), num_sample) == 0  
        % i��num_sample����������������㹻�Ĳ��������и���
        Lp_subsample = Lp_subsample+1;
        % imu������ʽ�������������ʱ������ȡƽ��
        vm = 0.5*ts_imu.*(imu_SD.fb_ref(i-num_sample+1 : i, 1:3)...
            + imu_SD.fb_ref(i-num_sample : i-1, 1:3));
        
        wm = 0.5*ts_imu.*(imu_SD.wb_ref(i-num_sample+1 : i, 1:3)...
            + imu_SD.wb_ref(i-num_sample : i-1, 1:3));
        
        % ����Բ׶/�������
        [ phim, dvbm ] = cnscl(wm, vm);
        
        % ������̬�仯��
        Cntn0 = Cntn0*rv2m(eth.winn*num_sample*ts_imu);
        Cbtb0 = Cbtb0*rv2m(phim);
        Cbn = Cntn0'*Cbn0*Cbtb0;
        
        att_cmpt_sv1(Lp_subsample, :) = m2att(Cbn)*deg;
        att_ref_sv(Lp_subsample, :) = att_SD(i, 1:3)*deg;
        att_err_sv(Lp_subsample, :) = att_cmpt_sv1(Lp_subsample, :) -...
            att_ref_sv(Lp_subsample, :);
    else
        % do nothing
    end
    
    % *** �Ľ���������� ***
    
    
    % �ϵ�λ��
    stophere = 1;
    
    % �´θ��µ�ʱ��
    crt = crt + ts_imu;
end
crt = crt - ts_imu;

% ���δʹ�ÿռ�
idx_sv(Lp+1:end, :) = [];
att_cmpt_sv1(Lp_subsample+1:end, :) = [];
att_ref_sv(Lp_subsample+1:end, :) = [];
att_err_sv(Lp_subsample+1:end, :) = [];

%% ��ͼ
Time_axis_subsample = (1:1:Lp_subsample)*num_sample*ts_imu;
num_fig = 1;
% ��̬����ֵ
figure(num_fig)
subplot(311)
plot(Time_axis_subsample, att_cmpt_sv1(:, 1));
xlabel('ʱ�� /s'); ylabel('pitch /deg'); hold on; 
plot(Time_axis_subsample, att_ref_sv(:, 1));
legend('����ֵ','�ο�ֵ');
subplot(312)
plot(Time_axis_subsample, att_cmpt_sv1(:, 2));
xlabel('ʱ�� /s'); ylabel('pitch /deg'); hold on; 
plot(Time_axis_subsample, att_ref_sv(:, 2));
subplot(313)
plot(Time_axis_subsample, att_cmpt_sv1(:, 3));
xlabel('ʱ�� /s'); ylabel('pitch /deg'); hold on; 
plot(Time_axis_subsample, att_ref_sv(:, 3));


% ��̬���
msplot(311, Time_axis_subsample, att_err_sv(:, 1), 'ʱ�� /s', 'pitch err');
msplot(312, Time_axis_subsample, att_err_sv(:, 2), 'ʱ�� /s', 'roll err');
msplot(313, Time_axis_subsample, att_err_sv(:, 3), 'ʱ�� /s', 'yaw err');
