%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Eileen Evans    8/24/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function calculating nearest distance between points and a line
%   (representing a trench). Distances given in units of R. 
%                ( 8/24/2021 , 12:26:36 am )
%
%   INPUT
%       1. R - radius of sphere (6371 km for Earth)
%       2. ptlon - column vector of longitudes of input points
%       3. ptlat - column vector of latitudes of input points
%       4. linelon - column vector of longitudes of line/trench
%       5. linelat - column vector of latitudes of line/trench
%
%   OUTPUT
%       1. dmin - minimum distance between trench segments and points
%       2. aztotrench - azimuth to closest point on trench
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dmin, aztoseg] = point_to_line_seg_latlon_nomt(R, ptlon,ptlat,linelon,linelat)

dall = NaN*ones(numel(ptlon),numel(linelon)-1);
azall = NaN*ones(numel(ptlon),numel(linelon)-1);
aztoseg = NaN*ones(numel(ptlon),1);
% dall = NaN*ones(numel(ptlon),numel(linelon)-1);

npts = numel(ptlon);

% figure; sphere()

for ii = 1:numel(linelon)-1
    
    v1lon = linelon(ii);
    v2lon = linelon(ii+1);
    v1lat = linelat(ii);
    v2lat = linelat(ii+1);
    
    
    % solution from Eric Nitardy, https://math.stackexchange.com/questions/23054/how-to-find-the-distance-between-a-point-and-line-joining-two-points-on-a-sphere
    Abar = [cosd(v1lat).*cosd(v1lon), cosd(v1lat).*sind(v1lon), sind(v1lat)];
    Bbar = [cosd(v2lat).*cosd(v2lon), cosd(v2lat).*sind(v2lon), sind(v2lat)];
    Xbar = [cosd(ptlat).*cosd(ptlon), cosd(ptlat).*sind(ptlon), sind(ptlat)];
    Xbar_mag = sqrt(Xbar(:,1).^2 + Xbar(:,2).^2 + Xbar(:,3).^2);
    
    n = cross(Abar,Bbar);
    nbar = repmat(n,numel(Xbar(:,1)),1);
    sinphivec = cross(nbar,Xbar);
    sinphimag = sqrt(sinphivec(:,1).^2 + sinphivec(:,2).^2 + sinphivec(:,3).^2);
    Phi = asin(sinphimag); % for some reason this doesn't work
    
    
    % Find projection onto great circle (https://math.stackexchange.com/questions/2981618/closest-point-on-line-segment-of-a-great-circle)
    Q1 = cross(Abar,Bbar);
    q1bar = repmat(Q1,numel(Xbar(:,1)),1);
    Q2 = cross(q1bar,Xbar);
    
    % should be cleverer here eventually to consider both directions about
    % the globe and pick the smallest.
    Q = cross(Q2,q1bar);
    qmag = sqrt(Q(:,1).^2 + Q(:,2).^2 + Q(:,3).^2);
    Qbar = Q./qmag;
    [ABdist,ABaz] = haversine_distance(R,v1lat,v1lon,v2lat,v2lon); % returns distance units
    ABdist = ABdist./R; % convert to rad
    [XAdist, XAaz] = haversine_distance(R,v1lat,v1lon,ptlat,ptlon);
    XAdist = XAdist./R;
    [XBdist, XBaz] = haversine_distance(R,v2lat,v2lon,ptlat,ptlon);
    XBdist = XBdist./R;
    
    %     Xbar_proj = -Xbar + dot(nbar,Xbar).*nbar;
    % %     Xbar_proj = Xbar - dot(nbar,Xbar);
    %
    %     Xbar_proj_mag = sqrt(Xbar_proj(:,1).^2 + Xbar_proj(:,2).^2 + Xbar_proj(:,3).^2);
    %
    %     xbar = Xbar_proj./Xbar_proj_mag;
    %     hold on; plot3(xbar(1,1),xbar(1,2),xbar(1,3),'*')
    xbar = Qbar;
    M = (Abar + Bbar) ./ 2;
    test = zeros(npts,1);
    for tt = 1:npts
        middistA = dot(Abar,M);
        middistX = dot(xbar(tt,:),M);
        if middistX > middistA % then xbar is between Abar and Bbar
            test(tt) = 1;
        end
    end
    test = logical(test);
    
%     if test
%         keyboard
%     end
%     if ii == 20
% %         % Figure for debugging
%         tt = 1; % replace with desired point
%         figure; sphere()
%     
%         hold on; plot3(Abar(1),Abar(2),Abar(3),'.k')
%         hold on; plot3(Bbar(1),Bbar(2),Bbar(3),'.k')
%         hold on; plot3(Xbar(tt,1),Xbar(tt,2),Xbar(tt,3),'*')
%         xlabel('x');ylabel('y');zlabel('z');
%         [latgc,longc] = track2('gc',v1lat,v1lon,v2lat,v2lon);
%         GCbar = [cosd(latgc).*cosd(longc), cosd(latgc).*sind(longc), sind(latgc)];
%         hold on; plot3(GCbar(:,1),GCbar(:,2),GCbar(:,3),'-k')
%     
%         [latgc2,longc2] = track1(v1lat,v1lon,ABaz);
%         GCbar2 = [cosd(latgc2).*cosd(longc2), cosd(latgc2).*sind(longc2), sind(latgc2)];
%         hold on; plot3(GCbar2(:,1),GCbar2(:,2),GCbar2(:,3),'-','Color',[0.5 0.5 0.5]);
%         hold on; plot3(Qbar(tt,1),Qbar(tt,2),Qbar(tt,3),'o')
%         % end debugging figure
%     
%         keyboard
%     end
        
    qlat = 90 - (acos(Qbar(:,3)).*180/pi);
    qlon = zeros(npts,1);
    qlon(Qbar(:,1)>0) = atan(Qbar(Qbar(:,1)>0,2)./Qbar(Qbar(:,1)>0,1));
    qlon(Qbar(:,1)<0) = atan(Qbar(Qbar(:,1)<0,2)./Qbar(Qbar(:,1)<0,1)) + pi;
    qlon(Qbar(:,1)==0) = pi./2;
    qlon = qlon.*180./pi - 360;
    
    [Qdist, Qaz] = haversine_distance(R,ptlat,ptlon,qlat,qlon);
    Qdist = Qdist/R;

    if test == 0
        endpointdists = [XAdist(~test) XBdist(~test)];
        endpointazs = [XAaz(~test) XBaz(~test)];
        [Qdist(~test),argmin] = min(endpointdists,[],2);
        
        Qnoseg = Qaz(~test);

        Qnoseg(argmin==1) = endpointazs(argmin==1,1);
        Qnoseg(argmin==2) = endpointazs(argmin==2,2);
        
        Qaz(~test) = Qnoseg;
    end
    
    dall(:,ii) = Qdist.*R;
    azall(:,ii) = Qaz;
    
end
[dmin, dargmin] = min(dall,[],2);
for jj = 1:npts
    aztoseg(jj) = azall(jj,dargmin(jj));
end
trenchseg = dargmin;
end