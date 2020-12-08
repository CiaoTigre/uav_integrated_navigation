%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   Trajectory simulator
%
% Outputs:
%   - avp_SD.(pos,vn,att):  Reference(true) data of pos, vel & att.
%   - imu_SD.(wm,vm):       Increment of angular & vel.
%
% References:
%   - @�����ߵ��㷨����ϵ���ԭ��,�Ϲ��� P229
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1



close all
clear
clc



% ����ȫ�ֱ���
gvar_earth;
% ����
ts = 0.01;

% ��ʼ��̬���ٶȡ�λ��
att0 = [0, 0, 90]'*arcdeg; % 
vn0 = [0, 0, 0]';
pos0 = [34*arcdeg, 108*arcdeg, 100]'; % lattitude, longtitude, height

% Ԥ��켣
%   w_pitch   w_roll   w_yaw   acc_y    time
% wat = [0,       0,      0,      0,      10      %��ֹ
%        0,       0,      0,      1,      10      %����
%        0,       0,      0,      0,      10      %����
%        5,       0,      0,      0,      4       %̧ͷ
%        0,       0,      0,      0,      10      %����
%        -5,      0,      0,      0,      4       %��ͷ
%        0,       0,      0,      0,      10      %����
%        0,       10,     0,      0,      1       %���
%        0,       0,      9,      0,      10      %ת��
%        0,       -10,    0,      0,      1       %���
%        0,       0,      0,      0,      10      %����
%        0,       0,      0,      -1,     10      %����
%        0,       0,      0,      0,      10    ];%��ֹ

   
% % Ԥ��켣 Loiter
% %   w_pitch   w_roll   w_yaw   acc_y    time
% wat = [ 0,       0,      0,       1,      5     %����
%         0,       0,      0,       0,      10    %����
%         3,       0,      0,       0,      3     %̧ͷ
%         0,       0,      6,       0,      60  	%����ת��
%         0,       0,      7,       0,      50  	%����ת��
%         0,       0,      8,       0,      40  	%����ת��
%         0,       0,      9,       0,      30  	%����ת��
%         0,       0,      10,       0,     20  	%����ת��
%        -3,       0,      11,       0,      3 	%��ͷ
%         0,       0,      11,       -1,     5];	%��

% % Ԥ��켣 Square
% %   w_pitch   w_roll   w_yaw   acc_y    time
wat = [ 0,       0,      0,       1,      4     %����
        0,       0,      0,       0,      30    %����
        0,       0,      0,       -1,      2    %����
        0,       0,      10,       0,     9     %ת��
        0,       0,      0,       1,      2     %����
        
        0,       0,      0,       0,      10    %����
        2,       0,      0,       0,      20    %����
        -2,       0,      0,       0,      20    %����
        0,       0,      0,       -1,      2    %����
        0,       0,      10,       0,     10    %ת��
        0,       0,      0,       1,      2     %����
        
                
        0,       0,      0,       0,      10    %����
        0,       0,      0,       0,      5    %����
        0,       0,      0,       0,      5    %����
        0,       0,      0,       0,      20    %����
        
        0,       0,      0,       -1,      2    %����
        0,       0,      10,       0,     10    %ת��
        0,       0,      0,       1,      2     %����
        
        0,       0,      0,       0,      30    %����
        0,       0,      0,       -1,      2    %����
        0,       0,      10,       0,     9     %ת��
        0,       0,      0,       1,      2     %����
        
        0,       0,      0,       0,      50    %����
        0,       0,      0,       -1,      5];  %����
   

% ����deg/s��ʾ�Ľ��ٶ�ת����rad/s
wat(:, 1:3) = wat(:, 1:3)*pi/180;

[avp_SD.att, avp_SD.vn, avp_SD.pos] = trjprofile(att0, vn0, pos0, wat, ts);
% �����õ�avp2imu��Ϊ���ϵ�av2imu���� @P229
[imu_SD.wm, imu_SD.vm] = avp2imu(avp_SD.att, avp_SD.vn, avp_SD.pos, ts);

% ��imu������"����" t=0ʱ��imu���,��(-ts,0]���ʱ���ڵĽ��������ٶ�������
% ������ʹ��imu_SD��avp_SDӵ����ͬ�ĳ��ȣ��������ݵĵ��á�
imu_SD.wm = [[0, 0, 0]; imu_SD.wm];
imu_SD.vm = [[0, 0, 0]; imu_SD.vm];

t_SD = sum(wat(:,5));

% ��������
% save('./UavTrajectorySim/UavTrajectory_Loiter.mat', 'avp_SD', 'imu_SD', 't_SD');
save('./UavTrajectorySim/UavTrajectory_Square.mat', 'avp_SD', 'imu_SD', 't_SD');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ͼ����
tt = (0 : length(avp_SD.att) - 1 )'*ts;
% ��̬
msplot(221, tt, avp_SD.att/arcdeg, 'Att/(\circ)');
legend('\it\theta', '\it\gamma', '\it\psi')
title('Attitutde');

% �ٶ�
msplot(222, tt, avp_SD.vn, 'Vel /m.s^{-1}');
legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
title('Velocity');

% ���λ��
msplot(223, tt, deltapos(avp_SD.pos), '\DeltaPos /m');
legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')
title('Relative position');

% ����λ�ã������Ǿ��ȣ�������γ�ȣ�ƽ��ͼ
subplot(2,2,4)
plot(avp_SD.pos(:, 2)/arcdeg, avp_SD.pos(:, 1)/arcdeg,'b','LineWidth',1); hold on;
grid on;
axis equal;
xlabel('\it\lambda\rm /(\circ)');
ylabel('\itL\rm /(\circ)');
plot(avp_SD.pos(1,2)/arcdeg, avp_SD.pos(1, 1)/arcdeg, 'co','LineWidth',3);
legend('True pos.','start')
title('Absolute position');

% ��ά�켣
figure;
plot3(avp_SD.pos(:,2)/arcdeg,avp_SD.pos(:,1)/arcdeg,avp_SD.pos(:,3), 'b', 'LineWidth',2);
xlabel('\it\lambda\rm /(\circ)');
ylabel('\itL\rm /(\circ)');
zlabel('\ith\rm /(m)')
hold on;
plot3(avp_SD.pos(1,2)/arcdeg, avp_SD.pos(1,1)/arcdeg, avp_SD.pos(1,3),'oc','LineWidth',3);
legend('True pos.','start')
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ģ���imu�����������Ϣ
msplot(121, tt, imu_SD.wm/ts/arcdeg, ...
    '\it\omega^b_{ib}\rm /(\circ.s^{-1})');
legend('\it\omega^b_{ibx}', '\it\omega^b_{iby}', '\it\omega^b_{ibz}');
title('Angular rate');

msplot(122, tt, imu_SD.vm/ts, '\itf^b\rm_{sf}/(m.s^{-2})');
legend('\itf^b\rm_{sf\itx}', '\itf^b\rm_{sf\ity}', '\itf^b\rm_{sf\itz}');
title('Acceleration');
