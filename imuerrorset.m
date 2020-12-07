function [ imu_err ] = imuerrorset( condition )
% ���ƣ�imu error set
% ���ܣ�����imu���
%
% Inputs:
%       condition: ��������
% Outputs:
%       imu_err: imu����������ɵ�struct,������
%            eb: ������ƫ
%            web: �Ƕ��������
%            db: �Ӽ���ƫ
%            wdb: �ٶ��������

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
gvar_earth;

% �������condition,����condition != ''��'zero'
if exist('condition', 'var') && ~strcmp(condition, '') ...
                             && ~strcmp(condition, 'zero') 
    switch condition
        case 'selfdefine'
            % *** �Զ��巽�� ***
            imu_err.case = 'selfdefine';
            imu_err.eb = [0.01, 0.01, 0.01]'*dph;
            imu_err.web = [0.01, 0.01, 0.01]'*dpsh;
            imu_err.db = [100, 100, 100]'*ug;
            imu_err.wdb = [1, 1, 1]'*ugpsHz;
        
        otherwise
            % *** Ĭ����� *** 
            % @�����ߵ��㷨����ϵ���ԭ�� P239 ��ϵ��������е��������
            imu_err.case = 'default';
            imu_err.eb = [0.01, 0.015, 0.02]'*dph;
            imu_err.web = [0.001, 0.001, 0.001]'*dpsh;
            imu_err.db = [80, 90, 100]'*ug;
            imu_err.wdb = [1, 1, 1]'*ugpsHz;
    end
    
else
    % *** zero error ***
    % ���û�и���condition����Ĭ��imu���Ϊ�� 
    imu_err.case = 'zero';
    imu_err.eb = [0, 0, 0]'*dph;
    imu_err.web = [0, 0, 0]'*dpsh;
    imu_err.db = [0, 0, 0]'*ug;
    imu_err.wdb = [0, 0, 0]'*ugpsHz;
    
end

end
