function [ kf ] = kfupdate( kf, Zk, time_measure_both )
% ���ƣ�Kalman filter update
% ���ܣ�
%
% Inputs:
%       kf: k-1ʱ�̵�kalman filter����
%       Zk: kʱ�̴�������õ�������Ϣ
%       time_measure_both: 
% Outputs:
%       kf: kʱ�̵�kalman filter����

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% T: ����ʱ����£�M: ����������£�B: ����ʱ����º��������
if nargin == 1
    % ���û���������룬��ֻ����ʱ�����
    time_measure_both = 'T';
elseif nargin == 2
    % ���������룬����ʱ����º��������
    time_measure_both = 'B';  
end


if time_measure_both == 'T' || time_measure_both == 'B'
    % *** ʱ����� ***
    % ״̬Ԥ��
    kf.Xkk_1 = kf.Phikk_1*kf.Xk;
    kf.Pkk_1 = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak'; 
    
else % time_measure_both == 'M'(ԭ��)
    % û��ģ����Ϣ��ֱ�Ӱ��ϸ�ʱ�̵���Ϣ����ʱ����½��
    kf.Xkk_1 = kf.Xk;
    kf.Pkk_1 = kf.Pk;
end


if time_measure_both == 'M' || time_measure_both == 'B'
    % *** ������� ***
    %  Pk|k-1*Hk'
    kf.PXZkk_1 = kf.Pkk_1*kf.Hk';
    % (Hk*Pk|k-1*Hk' + Rk)
    kf.PZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;
    % Pk|k-1*Hk'*inv(Hk*Pk|k-1*Hk' + Rk)
    kf.Kk = kf.PXZkk_1/kf.PZkk_1;
    kf.Xk = kf.Xkk_1 + kf.Kk*(Zk - kf.Hk*kf.Xkk_1);
    % Э�����С�㷨��
    kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZkk_1*kf.Kk'; 
    
else % time_measure_both == 'T'(ԭ��)
    % û��������Ϣ���൱��Kk=zeros(n),ֱ����ʱ����µĽ����Ϊ������½��
    kf.Xk = kf.Xkk_1;
    kf.Pk = kf.Pkk_1;
end

% P��Խǻ�
kf.Pk = (kf.Pk + kf.Pk')/2;

end
