function [  ] = msplot( mnp, x, y, xstr, ystr )
%% **************************************************************
% ���ƣ�
% ���ܣ���ͼ
% ����@ �����ߵ�ϵͳ����ϵ���ԭ�� P236
%
% Inputs:
%       mnp: subplot(mnp) = subplot(m, n, p)
%       x, y: x��y������Ӧ������
%       xstr, ystr: x��y���ע��

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ��ȫ��������copy
% ����ǵ�һ��Сͼ�����½�һ��figure
% subplot(mnp) = subplot(m, n, p),mnp����Ϊ��λ��
if mod(mnp, 10) == 1
    figure;
end

subplot(mnp);
plot(x, y);
grid on;
axis tight;

% x���y���ע��
if nargin == 4
    ystr = xstr;
    xstr = 't/s';
end

xlabel(xstr);   ylabel(ystr);

end
