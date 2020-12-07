function [ kf ] = kfupdate( kf, Zk, time_measure_both )
% 名称：Kalman filter update
% 功能：
%
% Inputs:
%       kf: k-1时刻的kalman filter参数
%       Zk: k时刻传感器测得的量测信息
%       time_measure_both: 
% Outputs:
%       kf: k时刻的kalman filter参数

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% T: 进行时间更新，M: 进行量测更新，B: 进行时间更新和量测更新
if nargin == 1
    % 如果没有量测输入，则只进行时间更新
    time_measure_both = 'T';
elseif nargin == 2
    % 有量测输入，进行时间更新和量测更新
    time_measure_both = 'B';  
end


if time_measure_both == 'T' || time_measure_both == 'B'
    % *** 时间更新 ***
    % 状态预测
    kf.Xkk_1 = kf.Phikk_1*kf.Xk;
    kf.Pkk_1 = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak'; 
    
else % time_measure_both == 'M'(原文)
    % 没有模型信息，直接把上个时刻的信息当作时间更新结果
    kf.Xkk_1 = kf.Xk;
    kf.Pkk_1 = kf.Pk;
end


if time_measure_both == 'M' || time_measure_both == 'B'
    % *** 量测更新 ***
    %  Pk|k-1*Hk'
    kf.PXZkk_1 = kf.Pkk_1*kf.Hk';
    % (Hk*Pk|k-1*Hk' + Rk)
    kf.PZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;
    % Pk|k-1*Hk'*inv(Hk*Pk|k-1*Hk' + Rk)
    kf.Kk = kf.PXZkk_1/kf.PZkk_1;
    kf.Xk = kf.Xkk_1 + kf.Kk*(Zk - kf.Hk*kf.Xkk_1);
    % 协方差减小算法：
    kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZkk_1*kf.Kk'; 
    
else % time_measure_both == 'T'(原文)
    % 没有量测信息，相当于Kk=zeros(n),直接用时间更新的结果作为量测更新结果
    kf.Xk = kf.Xkk_1;
    kf.Pk = kf.Pkk_1;
end

% P阵对角化
kf.Pk = (kf.Pk + kf.Pk')/2;

end
