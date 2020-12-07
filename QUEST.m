function [ q_optimal, Km ] = QUEST ( alpha, beta, Km_1 )
% 名称：QUternion ESTimate
% 功能：基于四元数估计的姿态最优姿态计算值
%
% Inputs:
%       qbn*alpha*qbn' = beta
%       Km_1: 即K(tm-1)，由获得第m次之前的m-1组（参考-量测）值构建的二次型矩阵
% Outputs:
%       q_optimal: 姿态四元数qbn的最优值
%       Km: 添加了本次量测的二次型矩阵
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 如果之前不存在量测信息，则将Km_1设置为0矩阵
if ~exist('Km_1', 'var')
    Km_1 = zero(4,4);
end

% 将输入矢量转化为列向量的形式
if size(alpha, 1) == 1
    alpha = alpha';
end
if size(beta, 1) == 1
    beta = beta';
end

% alpha和beta在四元数乘法中对应的系数矩阵：
Ma = [    0,      -alpha'
      alpha, -skew(alpha)];

Mb = [   0,     -beta'
      beta, skew(beta)];

% 构建二次型矩阵Km
Km = Km_1 + (Mb - Ma)'*(Mb - Ma);

% 计算Km的特征值和特征向量
[eigen_vector, eigen_value]=eig(Km);    
[~, number]=min(diag(eigen_value));
q=eigen_vector(:, number);   

% 四元数归一化:
q_optimal = q/norm(q, 2);

end
