%% **************************************************************
% ���ƣ�angle motion and imu output simulation (version 1.0)
% Function���û�Ԥ����̬���ó���ģ�����̬��������Ӧimu���
%   
% Inputs:
%       ang_motion: �û����õĽ��˶�
%       ang_motion.p(r,y)��ÿһ�б���һ�������˶���Ӧ�Ĳ���
%       theta = A*sin( w*t + phi) + k;
%       ang_motion.p = [A1(��ֵ), w1(���ٶ�), phi1(��ʼ��λ), k1(����ֵ)
%                 A2(��ֵ), w2(���ٶ�), phi2(��ʼ��λ), k2(����ֵ)]
%       Time: [total_time, ts, nn] 
% Outputs:
%       att : �����õĽ��˶�����Ӧ����̬
%       imu : ���� imu.fb_ref, imu.wb_ref, imu.fb_msr, imu.wb_msr

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1



%%
close all;
% ȫ�ֱ�������
gvar_earth;
pos = [46*arcdeg, 126*arcdeg, 100]';
eth = earth(pos, [0, 0, 0]');

total_time = 300;
ts = 0.01;
nn = 2;

% ���ý��˶�1
ang_motion.p = [2.5*arcdeg, 10*arcdeg, 0*arcdeg, 0*arcdeg];
ang_motion.r = [1.5*arcdeg, 10*arcdeg, 0*arcdeg, 0*arcdeg];
ang_motion.y = [10*arcdeg, 12*arcdeg, 0*arcdeg, 0*arcdeg
                6*arcdeg, 10*arcdeg, 0*arcdeg, 0*arcdeg
                8*arcdeg, 5*arcdeg, 0*arcdeg, 0*arcdeg];

% % ���ý��˶�2
% ang_motion.p = [0*arcdeg, 0*arcdeg, 0*arcdeg, 0*arcdeg];
% ang_motion.r = [0*arcdeg, 0*arcdeg, 0*arcdeg, 0*arcdeg];
% ang_motion.y = [10*arcdeg, 5*arcdeg, 0*arcdeg, 0*arcdeg];

% ���ô��������
% imu_err.eb = [0.1, 0.1, 0.1]'*dph;
% imu_err.web = [0.01, 0.01, 0.01]'*dpsh;
% imu_err.db = [500, 500, 500]'*ug;
% imu_err.wdb = [10, 10, 10]'*ugpsHz;
imu_err = imuerrorset('yanP239');

% Ϊ��̬��������洢�ռ�
att_SD = zeros(fix(total_time/ts) + 10, 4);
imu_SD.fb_ref = zeros(fix(total_time/ts) + 10, 4);
imu_SD.wb_ref = zeros(fix(total_time/ts) + 10, 4);
%% simulation
% ��ǰʱ��,ע��crt��0��ʼ��ʱ,att_SD����Ԫ��Ϊt=0ʱ�̵�������̬
crt = 0;
% ѭ������
i = 0;
while crt < total_time
    
    % ѭ������+1
    i = i+1;
    
    % �ֱ���ȡ������̬���˶���Ϣ
    % ����pitch(crt), roll(crt), yaw(crt)
    [pitch, wnb_p] = caculate_theta( ang_motion.p, crt );
    [roll, wnb_r] = caculate_theta( ang_motion.r, crt );
    [yaw, wnb_y] = caculate_theta( ang_motion.y, crt );
    
    % ��̬
    Cbn = a2mat([pitch, roll, yaw]');
    
    % wnbnx, ע��ý��ٶȲ�����ȫͶӰ��nϵ������ֻ��Ϊ�˷����������ʾ
    wnbnx = [wnb_p, wnb_r, wnb_y]';
    
    % ����wibb, Cnxb�е�nx��ʾ(n1, n2, b)����������ϵ
    Cnxb = [ cos(roll), 0, -cos(pitch)*sin(roll)
                    0,  1,            sin(pitch)
            sin(roll),  0,  cos(pitch)*cos(roll)];
    wibb_ref = Cnxb*wnbnx + Cbn'*eth.winn;
    
    % ���������˶����Ӽ����ӦΪbϵ���������ٶ�
    fb_ref = Cbn'*eth.gn;
    
    % ������
    att_SD(i, :) = [pitch, roll, yaw, crt];
    imu_SD.fb_ref(i, :) = [fb_ref', crt];
    imu_SD.wb_ref(i, :) = [wibb_ref', crt];
    
    % ����һ��������ʱ�䣬Ϊ�´μ�����׼��
    crt = crt + ts;
    
end
% ɾ���洢�����еĿ�λ
att_SD(i+1:end, :) = [];
imu_SD.fb_ref(i+1:end, :) = [];
imu_SD.wb_ref(i+1:end, :) = [];

% imu_ref��������imu_msr
[ imu_SD.fb_msr, imu_SD.wb_msr ] = imuadderr( imu_SD.fb_ref, imu_SD.wb_ref, ...
     imu_err.eb, imu_err.web, imu_err.db, imu_err.wdb, ts );

% ���ݱ���
save('angle_motion_data', 'att_SD', 'imu_SD');
%% ��ͼ
% ʱ����
Time_axis = (1 : 1 : i)*ts;

% ��̬��
msplot(311, Time_axis, att_SD(:, 1)*deg, 'ʱ�� /s', 'pitch /\circ');
title('��̬�Ǳ仯');
msplot(312, Time_axis, att_SD(:, 2)*deg, 'ʱ�� /s', 'roll /\circ');
msplot(313, Time_axis, att_SD(:, 3)*deg, 'ʱ�� /s', 'yaw /\circ');

% imu.fb_ref
msplot(321, Time_axis, imu_SD.fb_ref(:, 1), 'ʱ�� /s', 'fbx m/s^2');
title('�ӼƱ�׼���');
msplot(323, Time_axis, imu_SD.fb_ref(:, 2), 'ʱ�� /s', 'fby m/s^2');
msplot(325, Time_axis, imu_SD.fb_ref(:, 3), 'ʱ�� /s', 'fbz m/s^2');
% imu.wb_ref
msplot(322, Time_axis, imu_SD.wb_ref(:, 1)*deg, 'ʱ�� /s', 'wibb x rad/s');
title('���ݱ�׼���');
msplot(324, Time_axis, imu_SD.wb_ref(:, 2)*deg, 'ʱ�� /s', 'wibb y rad/s');
msplot(326, Time_axis, imu_SD.wb_ref(:, 3)*deg, 'ʱ�� /s', 'wibb z rad/s');

% imu fb_msr
msplot(321, Time_axis, imu_SD.fb_msr(:, 1), 'ʱ�� /s', 'fbx m/s^2');
title('�ӼƲ���ֵ');
msplot(323, Time_axis, imu_SD.fb_msr(:, 2), 'ʱ�� /s', 'fby m/s^2');
msplot(325, Time_axis, imu_SD.fb_msr(:, 3), 'ʱ�� /s', 'fbz m/s^2');
% imu.wb_msr
msplot(322, Time_axis, imu_SD.wb_msr(:, 1)*deg, 'ʱ�� /s', 'wibb x rad/s');
title('���ݲ���ֵ');
msplot(324, Time_axis, imu_SD.wb_msr(:, 2)*deg, 'ʱ�� /s', 'wibb y rad/s');
msplot(326, Time_axis, imu_SD.wb_msr(:, 3)*deg, 'ʱ�� /s', 'wibb z rad/s');