%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���ƣ�Global variable of earth
% @�����ߵ��㷨����ϵ���ԭ�� P231
% ȫ�ֱ���

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
global GM Re ff wie ge gp g0 ug arcdeg arcmin arcsec hur ...
    dph dpsh ugpsHz lsc deg min sec

% WGS-84 model
GM = 3.986004415e14;
Re = 6.378136998405e6;
wie = 7.2921151467e-5;

% ������� ff = (Re - Rp)/Re
ff = 1/298.257223563;
ee = sqrt(2*ff - ff^2);
e2 = ee^2;
Rp = (1 - ff)*Re;

ge = 9.780325333434361;
gp = 9.832184935381024;
g0 = ge; 
ug = g0*1e-6;

% �Ƕȡ��Ƿֺͽ���,��degת��Ϊ��Ӧ��rad:
% 1��= pi/180 rad, 1��= pi/180/60 rad, 1��= pi/180/60/60 rad
arcdeg = pi/180;
arcmin = arcdeg/60;
arcsec = arcmin/60;
% �ѻ���ת���ɽǶȡ��Ƿֺͽ���
deg = 1/arcdeg;
min = 1/arcmin;
sec = 1/arcsec;

% 1hour = 3600s
hur = 3600;
dph = arcdeg/hur;
dpsh = arcdeg/sqrt(hur);
ugpsHz = ug/sqrt(1);

% lsc: line shape and color
% ע��! lsc�е��ַ������ȱ�����ͬ���ո�Ҳ���볤��
lsc = [' -k'; ' -b'; ' -r'; '-.m'; '--g'; ' :c'];
