function [ kf ] = kfinit( Qk, Rk, P0, Phikk_1, Hk, Gammak )
% ���ƣ�Kalman filter initialization
% ���ܣ��������˲�����ʼ��
% ����@ �����ߵ�ϵͳ����ϵ���ԭ�� P238
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
% kf.m:����ά��; kf.n:״̬ά��
[kf.m, kf.n] = size(Hk);

% ��ʼ״̬�ͳ�ʼ״̬���Э�������
kf.Xk = zeros(kf.n, 1);
kf.Pk = P0;

% ״̬��������������
kf.Qk = Qk;
kf.Rk = Rk;

% ״̬ת�ƾ�����������
kf.Phikk_1 = Phikk_1;
kf.Hk = Hk;

% ״̬������������
if nargin < 6
    kf.Gammak = eye(kf.n);
else
    kf.Gammak = Gammak;
end

end
