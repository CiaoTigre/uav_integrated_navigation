function [ kf ] = kfinit( Qk, Rk, P0, Phikk_1, Hk, Gammak )
% 名称：Kalman filter initialization
% 功能：卡尔曼滤波器初始化
% 程序@ 捷联惯导系统与组合导航原理 P238
%
% Inputs:
%       Qk: 
%       Rk: 
%       P0: 
%       Phikk_1: 
%       Hk: 
%       Gammak: 
% Outputs:
%       kf: 

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% kf.m:量测维数; kf.n:状态维数
[kf.m, kf.n] = size(Hk);

% 初始状态和初始状态误差协方差矩阵
kf.Xk = zeros(kf.n, 1);
kf.Pk = P0;

% 状态噪声和量测噪声
kf.Qk = Qk;
kf.Rk = Rk;

% 状态转移矩阵和量测矩阵
kf.Phikk_1 = Phikk_1;
kf.Hk = Hk;

% 状态噪声驱动矩阵
if nargin < 6
    kf.Gammak = eye(kf.n);
else
    kf.Gammak = Gammak;
end

end
