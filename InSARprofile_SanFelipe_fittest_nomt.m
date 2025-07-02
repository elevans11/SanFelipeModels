%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%    Eileen Evans    7/15/2024 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  transect parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

transectlength = 12; %km
transectwidth = 1; % km
doffset = transectlength/2;
range = [-116.5 -116 33 33.34];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

svec = -(0:0.1:5); % uplift rates (mm/yr)
ldvec = 0:.1:10; % locking depths (km)
dipvec = 80; % dips to test (degrees)

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
IN = inpolygon(xo,yo,[range(1) range(2) range(2) range(1)],[range(3) range(3) range(4) range(4)]);
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
dx_km = d*cosd(theta); % distance in the longitude direction (km)
dy_km = d*sind(theta); % distance in the latitude direction (km)
dx = dx_km; % m
dy = dy_km; % m

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
R = 6371;
[dmin, aztoseg] = point_to_line_seg_latlon_nomt(R, tslon,tslat,[S1(1); S2(1)],[S1(2); S2(2)]);
[dmin_s, aztoseg_s] = point_to_line_seg_latlon_nomt(R, slon,slat,[S1(1); S2(1)],[S1(2); S2(2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Add model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

W = 30000; % very deep fault (km)
nx = 100;
x = dmin-(transectlength/2);
obs = tsup+1;

a = flipud(parula(numel(svec)));
dd = 1;

MISFIT = NaN*ones(numel(ldvec),numel(svec));
tic
for ll = 1:numel(ldvec)
    ld = ldvec(ll);
    for ss = 1:numel(svec)
        dip = dipvec(dd)*pi/180; % degrees
        xf0 = 0;
        w0 = ld./sin(dip); % fault width sin(dip)
        xf1 = w0*cos(dip);
        yf0 = 0;
        yf1 = w0*sin(dip);
        xf2 = xf1 + W*cos(dip);
        yf2 = yf1 + W*sin(dip);
        

        slip = (svec(ss))/(sin(dipvec(dd)*pi/180));

        U1 = zeros(numel(x),1);
        U2 = zeros(numel(x),1);
        for ii = 1:numel(x)
            [u1, u2] = Thrust2DPartials(x(ii),0,slip,xf1,-yf1,xf2,-yf2,0);
            U1(ii) = u1;
            U2(ii) = u2;
        end

        profile = U2-U2(1);
        MISFIT(ll,ss) = rms((obs-profile),'omitnan');
    end
end
toc

[a,b] = find(MISFIT == min(min(MISFIT)));

figure; contourf(abs(svec),ldvec,MISFIT,100,'Lines','none')
hold on; contour(abs(svec),ldvec,MISFIT,0:.5:10,'LineColor','w','LineWidth',1);
clim([0 3])
ch = colorbar('horizontal');
ch.Color = 'k';
ch.Label.Color = 'k';
ch.Label.String = 'rms misfit (mm/yr)';
ch.TickLength = 0.033;

xlabel('vertical uplift rate (mm/yr)');
ylabel('locking depth (km)');
set(gca,'FontSize',16)

ld_bf = ldvec(a); % best fitting locking depth
svec_bf = svec(unique(b)); % best fitting uplift rate

