%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%    Eileen Evans    11/11/2013 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%   Describe purpose of script/function here. 
%                ( 11/11/2013 , 10:22:02 am ) 
% 
%   INPUT 
%       1. x1 x-locations of observations
%       2. x2 z-locations of observations (if calculating at depth)
%       3. s = slip
%       4. xf1 = x-location of top of dislocation
%       5. yf1 = z-location
%       6. d = depth of top of dislocation (located at [x,d])
%       7. l = length of dislocation
%       8. ftog - if 1, plots a figure of dislocations
% 
%   OUTPUT 
%       1. Output one here 
% 
%   Outline 
%       1.  
%       2.  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [U1, U2] = Thrust2DPartials(x1,x2,s,xf1,yf1,xf2,yf2,ftog) 
% close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Setup  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% x1 positive right
%%% x2 positive up

xi1                                 = [xf1 xf2];
xi2                                 = [yf1 yf2];

dip                                 = -atan2((yf2-yf1),(xf2-xf1));

s1                                  = s*cos(dip);
s2                                  = s*sin(dip);

nu                                  = 0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Calculate Displacements  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

u1 = zeros(numel(x1),numel(xi1));
u2 = zeros(numel(x1),numel(xi1));
theta1keep = zeros(numel(x1),numel(xi1));
theta2keep = zeros(numel(x1),numel(xi1));
theta1R = zeros(numel(x1),numel(xi1));
theta2R = zeros(numel(x1),numel(xi1));


% figure;
for jj = 1:numel(xi1);
    
    RotMat = [cos(dip-pi/2+pi) sin(dip-pi/2+pi); -sin(dip-pi/2+pi) cos(dip-pi/2+pi)];
    ToRotate = [x1-xi1(jj) x2-xi2(jj)];
    ToRotate2 = [xi1'-xi1(jj) xi2'-xi2(jj)];
    
    Rotated = ToRotate*RotMat;
    Rotated2 = ToRotate2*RotMat;
    
    XR = Rotated(:,1);
    YR = Rotated(:,2);
    XR2 = Rotated2(:,1);
    YR2 = Rotated2(:,2);    
    
    theta1R(:,jj)                              = atan2((XR - XR2(jj)),(YR - YR2(jj)));
    
    RotMat3 = [cos(dip-pi/2 + pi) sin(dip-pi/2+ pi); -sin(dip-pi/2+ pi) cos(dip-pi/2+ pi)];

    ToRotate = [x1-xi1(jj) x2+xi2(jj)];
    Rotated = ToRotate*RotMat3;
    
    ToRotate3 = [xi1'-xi1(jj) -xi2'+xi2(jj)];
    Rotated3 = ToRotate3*RotMat3;
    
    XR = Rotated(:,1);
    YR = Rotated(:,2);
    XR3 = Rotated3(:,1);
    YR3 = Rotated3(:,2);
    
    theta2R(:,jj)                              = atan2((XR - XR3(jj)),(YR + YR3(jj)));
                                    
    for ii = 1:numel(x1)
        theta1                              = theta1R(ii,jj);
        theta2                              = theta2R(ii,jj);
        
        
        theta1keep(ii,jj)                              = theta1;
        theta2keep(ii,jj)                              = theta2;
        
        
        r1                                  = sqrt((x1(ii)- xi1(jj))^2 + (x2(ii)-xi2(jj))^2);
        r2                                  = sqrt((x1(ii)- xi1(jj))^2 + (x2(ii)+xi2(jj))^2);
        
        term1 = ((1 - nu)/2)*(theta2 - theta1);
        term2 = ((x1(ii)-xi1(jj))*(x2(ii) - xi2(jj)))/(4*r1.^2);
        term3 = ((x1(ii)-xi1(jj))*(x2(ii) + (3-4*nu)*xi2(jj)))/(4*r2^2);
        term4 = (xi2(jj)*x2(ii)*(x1(ii)-xi1(jj))*(x2(ii)+xi2(jj)))/(r2^4);
        
        term5 = ((1-2*nu)/4)*(log10(r2/r1));
        term6 = (x2(ii)-xi2(jj))^2/(4*r1^2);
        term7 = (x2(ii)^2 + xi2(jj)^2 - 4*(1-nu)*xi2(jj)*(x2(ii)+xi2(jj)))/(4*r2^2);
        term8 = (x2(ii)*xi2(jj)*(x2(ii)+xi2(jj))^2)/r2^4;
        
        u1(ii,jj) = (s1/(pi*(1-nu))).*(term1 + term2 - term3 + term4)...
            + (s2/(pi*(1-nu))).*(term5 - term6 + term7 + term8);
        
        term9 = ((x2(ii)+xi2(jj))^2 - 2*xi2(jj)^2 - 2*(1-2*nu)*xi2(jj)*(x2(ii)+xi2(jj)))/(4*r2^2);
        term10 = ((1 - nu)/2)*(theta1 - theta2);
        
        u2(ii,jj) = (s1/(pi*(1-nu))).*(term5 + term6 - term9 + term8)...
            + (s2/(pi*(1-nu))).*(term10 + term2 - term3 - term4);
    end

end

U1 = -u1(:,1:end-1) + u1(:,2:end);
U2 = -u2(:,1:end-1) + u2(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if ftog
figure; 
quiver(X1(:),X2(:),U1,U2);
axis equal;
axis tight;
hold on; plot(xi1,xi2,'-or')
end
















