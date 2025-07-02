%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%    Eileen Evans    7/15/2024 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  transect parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

transectlength = 12; % km total length, centered at sample EA23-9
transectwidth = 1; % km
doffset = transectlength/2;

range = [-116.5 -116 33 33.34]; % figure range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

svec = -5:1:-1; % uplift rate (mm/yr), negative for thrust motion
dipvec = 80; % dip (degrees)
ld = 3; % locking depth (km)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  load decomposed InSAR (Xu et al., 2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[xe,ye,ze] = grdread2('Xu_et_al_2021_vel_decomp/East_defo.grd');
[xn,yn,zn] = grdread2('Xu_et_al_2021_vel_decomp/North_defo.grd');
[xu,yu,zu] = grdread2('Xu_et_al_2021_vel_decomp/Up_defo.grd');
[X,Y] = meshgrid(xe,ye);
xo = X(:);
yo = Y(:);
east = ze(:);
north = zn(:);
up = zu(:);
IN = inpolygon(xo,yo,[range(1) range(2) range(2) range(1)],[range(3) range(3) range(4) range(4)]); % keep only data in range
xo = xo(IN);
yo = yo(IN);
east = east(IN);
north = north(IN);
up = up(IN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  load samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

sfile = 'samplelocations.txt';
fid = fopen(sfile);
T = textscan(fid,'%s\t%f\t%f','headerlines',1);
fclose(fid);
slon = T{3};
slat = T{2};
sid = T{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Hard-code transect orientation perpendicular to San Felipe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SFaz = 120.6539;
theta = 180 - SFaz;
d = transectlength/2;
% dx_km = d*cosd(theta); % distance in the longitude direction (km)
% dy_km = d*sind(theta); % distance in the latitude direction (km)
% dx = dx_km; % m
% dy = dy_km; % m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  convert horizontal to fault parallel/perpendicular
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fault_perp  = east*cosd(SFaz) - north*sind(SFaz);
fault_par   = east*sind(SFaz) + north*cosd(SFaz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Transect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 6371; % radius of the Earth, km

[Blat,Blon] = reckon_sphere(R,slat(9),slon(9),SFaz+90,d);
[Bprimelat,Bprimelon] = reckon_sphere(R,slat(9),slon(9),SFaz-90,d);

B = [Blon Blat];
Bprime = [Bprimelon; Bprimelat];

Transectaz = SFaz-90;
Transectlen = transectlength;

[S1lat,S1lon] = reckon_sphere(R,Blat,Blon,SFaz,transectwidth);
[S2lat,S2lon] = reckon_sphere(R,Blat,Blon,SFaz-180,transectwidth);
[S3lat,S3lon] = reckon_sphere(R,Bprimelat,Bprimelon,SFaz-180,transectwidth);
[S4lat,S4lon] = reckon_sphere(R,Bprimelat,Bprimelon,SFaz,transectwidth);

S1 = [S1lon, S1lat];
S2 = [S2lon, S2lat];
S3 = [S3lon, S3lat];
S4 = [S4lon, S4lat];

%%% Keep only sar data in the transect swath;
IN = inpolygon(xo,yo,[S1(1) S2(1) S3(1) S4(1) S1(1)],[S1(2) S2(2) S3(2) S4(2) S1(2)]);
tslon = xo(IN);
tslat = yo(IN);
tsup = up(IN);
tspar = fault_par(IN);
tsperp = fault_perp(IN);

%%% Data distance along profile
[dmin, aztoseg] = point_to_line_seg_latlon_nomt(R, tslon,tslat,[S1(1); S2(1)],[S1(2); S2(2)]);
[dmin_s, aztoseg_s] = point_to_line_seg_latlon_nomt(R, slon,slat,[S1(1); S2(1)],[S1(2); S2(2)]);

figure('Position',[1 1 700 730]); 
subplot(3,1,1)
hold on;
plot(dmin,tspar,'.'); % plot distance along transect vs. LOS
xline(dmin_s,'--k');

axis tight
xlabel('distance from B (km)');
ylabel('parallel (mm/yr)');
set(gca,'FontSize',16)

subplot(3,1,2)
hold on;
plot(dmin,tsperp,'.'); % plot distance along transect vs. LOS
xline(dmin_s,'--k');

axis tight
xlabel('distance from B (km)');
ylabel('perpendicular (mm/yr)');
set(gca,'FontSize',16)

subplot(3,1,3)
hold on;
plot(dmin,tsup,'.'); % plot distance along transect vs. LOS
xline(dmin_s,'--k');

axis tight
xlabel('distance from B (km)');
ylabel('vertical (mm/yr)');
set(gca,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Add model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

W = 30000; % very deep fault (km)
nx = 100;
x = linspace(-transectlength/2,transectlength/2,nx);

figure('Position',[50    50   635   730]); 
a = flipud(parula(numel(svec)));

% Plot fault geometry
subplot(2,1,2); hold on; 

U1 = zeros(numel(x),numel(dipvec),numel(svec));
U2 = zeros(numel(x),numel(dipvec),numel(svec));
for dd = 1:numel(dipvec)
    dip = dipvec(dd)*pi/180; % degrees
    xf0 = 0;
    w0 = ld./sin(dip); % fault width sin(dip)
    xf1 = w0*cos(dip);
    yf0 = 0;
    yf1 = w0*sin(dip);
    xf2 = xf1 + W*cos(dip);
    yf2 = yf1 + W*sin(dip);
    for ss = 1:numel(svec)

        slip = (svec(ss))/(sin(dipvec(dd)*pi/180));

        for ii = 1:numel(x)
            [u1, u2] = Thrust2DPartials(x(ii),0,slip,xf1,-yf1,xf2,-yf2,0);
            U1(ii,dd,ss) = u1;
            U2(ii,dd,ss) = u2;
        end
    end
    plot([xf0 xf1],-[yf0 yf1],'--','LineWidth',2,'Color','r');
    plot([xf1 xf2],-[yf1 yf2],'-','LineWidth',2,'Color','r');
    
end
ylim([-20 5])
xlim([-transectlength/2 transectlength/2])
yline(0,'k','LineWidth',2)

xlabel('distance from (locked) fault trace (km)');
ylabel('depth (km)');
set(gca,'FontSize',16);
textcells = {'5 mm/yr', '4 mm/yr','3 mm/yr','2 mm/yr','1 mm/yr'};
axis equal; % no vertical exaggeration

% Plot data and model
subplot(2,1,1)
hold on;
plot(dmin-(transectlength/2),tsup+1,'.','Color',0.5*[1 1 1]); % plot distance along transect vs. LOS
xline(dmin_s-(transectlength/2),'--k');

for dd = 1:numel(dipvec)
    for ss = 1:numel(svec)
        profile = squeeze(U2(:,dd,ss));
        profile = profile - profile(1);
        plot(x,profile,'-','Color',a(ss,:),'LineWidth',2);
        text(7,profile(end)-0.2,textcells{ss},'FontSize',16,'Color',a(ss,:));
    end    
end
colormap(a);
hold on;
xline(0,'r')
xlim([-transectlength/2 transectlength/2])
ylabel('vertical uplift (mm/yr)');
xlabel('distance from (locked) fault trace (km)');
set(gca,'FontSize',16);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fs                           = 18;
sw                           = 2;
ms                           = 3;
cbarystart                   = 0.17;
whl                          = 8;
mp                           = 'm_proj(''mercator'', ''long'', range(1:2) ,''lat'', range(3:4));';
mg1a                         = 'm_grid(''linestyle'', ''none'', ''tickdir'', ''out'', ''yaxislocation'', ''left'', ''xaxislocation'', ''bottom'', ''xlabeldir'', ''middle'',''ticklen'', 0.01, ''FontSize'', fs);';

figure('Position',[1 100 800 750]); hold on;
eval(mp); eval(mg1a);
hold on; 
m_scatter(xo,yo,10,up+1,'filled');

clim([-10 10])
colormap(bluewhitered)

hold on;
m_plot(slon,slat,'sk');
m_plot([B(1) Bprime(1)],[B(2) Bprime(2)],'-k','LineWidth',3);
m_plot([S1(1) S2(1) S3(1) S4(1) S1(1)],[S1(2) S2(2) S3(2) S4(2) S1(2)],'-k');
ch = colorbar('horizontal');
set(get(ch,'xlabel'),'String','vertical velocity (mm/yr)','FontSize',fs,'FontName','Times');
set(ch,'FontSize',fs);
m_text([Blon Bprimelon]+0.01,[Blat Bprimelat],{'B','B'''},'FontSize',fs);


% %%%
% figure('Position',[1 100 800 750]); hold on;
% eval(mp); eval(mg1a);
% hold on; 
% m_scatter(xo,yo,10,fault_perp - mean(fault_perp,'all','omitnan'));
% clim([-10 10])
% colormap(bluewhitered)
% 
% 
% for ss = 1:numel(Section)
%     m_plot(Section(ss).X+360,Section(ss).Y,'-','LineWidth',2,'Color',0*[1 1 1])
% end
% hold on;
% m_plot(SanFelipe(:,1),SanFelipe(:,2),'r');
% m_plot(slon,slat,'sk');
% m_plot([B(1) Bprime(1)],[B(2) Bprime(2)],'-k','LineWidth',3);
% m_plot([S1(1) S2(1) S3(1) S4(1) S1(1)],[S1(2) S2(2) S3(2) S4(2) S1(2)],'-k');
% title('fault-perpendicular')
% 
% %%%
% figure('Position',[1 100 800 750]); hold on;
% eval(mp); eval(mg1a);
% hold on; 
% m_scatter(xo,yo,10,fault_par - mean(fault_par,'all','omitnan'));
% clim([-10 10])
% colormap(bluewhitered)
% 
% 
% for ss = 1:numel(Section)
%     m_plot(Section(ss).X+360,Section(ss).Y,'-','LineWidth',2,'Color',0*[1 1 1])
% end
% hold on;
% m_plot(SanFelipe(:,1),SanFelipe(:,2),'r');
% m_plot(slon,slat,'sk');
% m_plot([B(1) Bprime(1)],[B(2) Bprime(2)],'-k','LineWidth',3);
% m_plot([S1(1) S2(1) S3(1) S4(1) S1(1)],[S1(2) S2(2) S3(2) S4(2) S1(2)],'-k');
% title('fault-parallel')


