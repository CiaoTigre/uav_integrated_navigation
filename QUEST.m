function [ q_optimal, Km ] = QUEST ( alpha, beta, Km_1 )
% ���ƣ�QUternion ESTimate
% ���ܣ�������Ԫ�����Ƶ���̬������̬����ֵ
%
% Inputs:
%       qbn*alpha*qbn' = beta
%       Km_1: ��K(tm-1)���ɻ�õ�m��֮ǰ��m-1�飨�ο�-���⣩ֵ�����Ķ����;���
% Outputs:
%       q_optimal: ��̬��Ԫ��qbn������ֵ
%       Km: ����˱�������Ķ����;���
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ���֮ǰ������������Ϣ����Km_1����Ϊ0����
if ~exist('Km_1', 'var')
    Km_1 = zero(4,4);
end

% ������ʸ��ת��Ϊ����������ʽ
if size(alpha, 1) == 1
    alpha = alpha';
end
if size(beta, 1) == 1
    beta = beta';
end

% alpha��beta����Ԫ���˷��ж�Ӧ��ϵ������
Ma = [    0,      -alpha'
      alpha, -skew(alpha)];

Mb = [   0,     -beta'
      beta, skew(beta)];

% ���������;���Km
Km = Km_1 + (Mb - Ma)'*(Mb - Ma);

% ����Km������ֵ����������
[eigen_vector, eigen_value]=eig(Km);    
[~, number]=min(diag(eigen_value));
q=eigen_vector(:, number);   

% ��Ԫ����һ��:
q_optimal = q/norm(q, 2);

end
