function [ att_err_normalized ] = atterrnorml( att_err )
% ����: Attitude Error Normalization (version 1.0)
% Function: ��̬����׼�� (Ŀǰ��������������)
%
% Inputs:
%       att_err: δ��׼������̬�����
% Outputs:
%       att_err_normalized: ��׼���������̬�����

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1


% �ֱ���ȡ������̬�����
pitch_err = att_err(1);
roll_err = att_err(2);
yaw_err = att_err(3);

% pitch �Ƕȷ�ΧΪ: [-pi/2, pi/2]
% roll �Ƕȷ�ΧΪ: (-pi, pi]
% yaw �Ƕȷ�ΧΪ: [0, 2*pi)

% ֱ�Ӽ�����ĺ��������(-2*pi, 2*pi),������ϣ�����������(-pi, pi]
% ��������������ֵ����pi���������׼��
if norm(yaw_err) > pi
    % �����������Ϊ��
    if sign(yaw_err) == 1
        yaw_err = yaw_err - 2*pi;
    else
        % �����������pi���ҷ���Ϊ��
        yaw_err = yaw_err + 2*pi;
    end
end

att_err_normalized = [pitch_err, roll_err, yaw_err]';

end

