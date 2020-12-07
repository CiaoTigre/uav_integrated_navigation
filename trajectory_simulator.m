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
att0 = [0, 0, 90]'*arcdeg;
vn0 = [0, 0, 0]';
pos0 = [34*arcdeg, 108*arcdeg, 100]'; % lattitude, longtitude, height

% Ԥ��켣
%   w_pitch   w_roll   w_yaw   vb_y    time
wat = [0,       0,      0,      0,      10      %��ֹ
       0,       0,      0,      1,      10      %����
       0,       0,      0,      0,      10      %����
       5,       0,      0,      0,      4       %̧ͷ
       0,       0,      0,      0,      10      %����
       -5,      0,      0,      0,      4       %��ͷ
       0,       0,      0,      0,      10      %����
       0,       10,     0,      0,      1       %���
       0,       0,      9,      0,      10      %ת��
       0,       -10,    0,      0,      1       %���
       0,       0,      0,      0,      10      %����
       0,       0,      0,      -1,     10      %����
       0,       0,      0,      0,      10    ];%��ֹ

   
% % Ԥ��켣
% %   w_pitch   w_roll   w_yaw   vb_y    time
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
%        0,       0,      9,      0,      10      %����
%        0,       10,      0,      0,      1      %����
%        0,       0,      0,      -1,     10      %����
%        0,       0,      0,      0,      10    ];%��ֹ
   

% ����deg/s��ʾ�Ľ��ٶ�ת����rad/s
wat(:, 1:3) = wat(:, 1:3)*pi/180;

[avp_SD.att, avp_SD.vn, avp_SD.pos] = trjprofile(att0, vn0, pos0, wat, ts);
% �����õ�avp2imu��Ϊ���ϵ�av2imu���� @P229
[imu_SD.wm, imu_SD.vm] = avp2imu(avp_SD.att, avp_SD.vn, avp_SD.pos, ts);

% ��imu������"����" t=0ʱ��imu���,��(-ts,0]���ʱ���ڵĽ��������ٶ�������
% ������ʹ��imu_SD��avp_SDӵ����ͬ�ĳ��ȣ��������ݵĵ��á�
imu_SD.wm = [[0, 0, 0]; imu_SD.wm];
imu_SD.vm = [[0, 0, 0]; imu_SD.vm];

% ��������
save('trajectory_simulator_data.mat', 'avp_SD', 'imu_SD');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ��ͼ
tt = (0 : length(avp_SD.att) - 1 )'*ts;

msplot(221, tt, avp_SD.att/arcdeg, 'Att/(\circ)');
legend('\it\theta', '\it\gamma', '\it\psi')
title('Attitutde');

msplot(222, tt, avp_SD.vn, 'Vel /m.s^{-1}');
legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
title('Velocity');

msplot(223, tt, deltapos(avp_SD.pos), '\DeltaPos /m');
legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')
title('Relative position');

% �����Ǿ��ȣ�������γ��
msplot(224, avp_SD.pos(:, 2)/arcdeg, avp_SD.pos(:, 1)/arcdeg, ...
    '\it\lambda\rm /(\circ)', '\itL\rm /(\circ)');    hold on
plot(avp_SD.pos(1,2)/arcdeg, avp_SD.pos(1, 1)/arcdeg, 'ro');
legend('trajectory','start position')
title('Absolute position');


% imu�����Ϣ��ͼ
msplot(121, tt, imu_SD.wm/ts/arcdeg, ...
    '\it\omega^b_{ib}\rm /(\circ.s^{-1})');
legend('\it\omega^b_{ibx}', '\it\omega^b_{iby}', '\it\omega^b_{ibz}');
title('Angular rate');

msplot(122, tt, imu_SD.vm/ts, '\itf^b\rm_{sf}/(m.s^{-2})');
legend('\itf^b\rm_{sf\itx}', '\itf^b\rm_{sf\ity}', '\itf^b\rm_{sf\itz}');
title('Acceleration');
