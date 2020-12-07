function [ imu_err ] = imuerrorset( condition )
% 名称：imu error set
% 功能：设置imu误差
%
% Inputs:
%       condition: 器件类型
% Outputs:
%       imu_err: imu器件参数组成的struct,包含：
%            eb: 陀螺零偏
%            web: 角度随机游走
%            db: 加计零偏
%            wdb: 速度随机游走

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
gvar_earth;

% 如果存在condition,并且condition != ''和'zero'
if exist('condition', 'var') && ~strcmp(condition, '') ...
                             && ~strcmp(condition, 'zero') 
    switch condition
        case 'selfdefine'
            % *** 自定义方案 ***
            imu_err.case = 'selfdefine';
            imu_err.eb = [0.01, 0.01, 0.01]'*dph;
            imu_err.web = [0.01, 0.01, 0.01]'*dpsh;
            imu_err.db = [100, 100, 100]'*ug;
            imu_err.wdb = [1, 1, 1]'*ugpsHz;
        
        otherwise
            % *** 默认情况 *** 
            % @捷联惯导算法与组合导航原理 P239 组合导航程序中的器件误差
            imu_err.case = 'default';
            imu_err.eb = [0.01, 0.015, 0.02]'*dph;
            imu_err.web = [0.001, 0.001, 0.001]'*dpsh;
            imu_err.db = [80, 90, 100]'*ug;
            imu_err.wdb = [1, 1, 1]'*ugpsHz;
    end
    
else
    % *** zero error ***
    % 如果没有给出condition，则默认imu误差为零 
    imu_err.case = 'zero';
    imu_err.eb = [0, 0, 0]'*dph;
    imu_err.web = [0, 0, 0]'*dpsh;
    imu_err.db = [0, 0, 0]'*ug;
    imu_err.wdb = [0, 0, 0]'*ugpsHz;
    
end

end
