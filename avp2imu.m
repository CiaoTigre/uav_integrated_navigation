function [ wm, vm ] = avp2imu( att, vn, pos, ts )
% 名称：using Attitude & Velocity inverse the IMU output  %%%(Practice version)
% Function: 反演生成imu输出数据  完全没有改动的版本
% 
% Inputs:
%       att: pitch, roll, yaw   (rad)
%       vn: vE, vN, vU          (m/s)
%       pos: latitude, longitude, high       (deg, deg, m)
%       ts: step time 
% Outpus:
%       wm: 角增量
%       vm: 速度增量

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
wm0 = zeros(3, 1);
vm0 = zeros(3, 1);
I33 = eye(3); 

wm = att(2:end, :);
vm = att(2:end, :);

for k = 2 : length(att)
    
    % 用vn和pos计算地球导航参数
    eth = earth((pos(k-1,:) + pos(k, :))/2, (vn(k-1, :) + vn(k, :))/2);
    
    qbb = qmul( qmul( qconj( a2qua(att(k-1, :) ) ), rv2q(eth.winn*ts)), a2qua(att(k, :) ) );
    phim = q2rv(qbb);
    
    % Note: matlab执行除法A/B时,实际上会进行矩阵操作：A*inv(B), (A,B为标量
    % 时也不例外)。通常矩阵乘法 A*B ≠ B*A, 因此矩阵除法分为"左除"和
    % "右除"，即:A/B = A*inv(B), A\B = inv(A)*B
    wm1 = (I33 + skew(1/12*wm0))\phim;
    
    dvnsf = vn(k, :)' - vn(k-1, :)' - eth.gcc*ts;
    Cbn0 = a2mat(att(k-1, :)');
    
    vm1 = (I33 + 1/2*skew(1/6*wm0 + wm1))\...
        (Cbn0'*(I33 + skew(eth.winn*ts/2))*dvnsf - 1/12*cross(vm0, wm1));
    
    wm(k-1, :) = wm1';
    vm(k-1, :) = vm1;
    wm0 = wm1;
    vm0 = vm1;
    
end
