% Written by Erik Lovekin, 1/20/2024
% Calculate the dimensions of the rail, drivetrain, and sliders so that the
% drivetrain doesn't derail or get stuck in t-slot.
% For derailing calculations, max wheel displacement calculated by considering
% drivetrain on circular track.
% For slider calcs, calculate longest slider that can turn about both
% inner and outer t-slots. Assumes rectangular slider with 90 degree
% corners. Inner corner is driving constraint.
close all;

%% Define Dimensions of Drivetrain and Beam
intom = 0.0254;% m/in
ls1 = 0.311*intom;%m, distance from back of drivetrain to back slider
la1 = 1.238*intom;%m, distance from back of drivetrain to back axle
la2 = 3.139*intom;%m, distance from back of drivetrain to front axle
ls2 = 4.075*intom;%m, distance from back of drivetrain to front slider
ld  = 4.386*intom;%m, distance from back of drivetrain to front of drivetrain
w = 1.252*intom;%m, wheelbase width
rw = 0.846*intom;%m, radius of wheel
tr = 0.217*intom;%m, thickness of wheel
ha = 0.276*intom;%m, distance between axle and bottom of drivetrain
hw = rw-ha*intom;%m, clearance below drivetrain on flat rail

lw1 = la1-ls1;%m, distance from back slider to back axle
lw2 = la2-ls1;%m, distance from back slider to front axle
ls = ls2-ls1;%m, distance from back slider to front slider

% Dimensions of Rail
rc = 0.1;%m radius passing through center of beam
tr = 2*intom;%m thickness of beam
ts = 0.4*intom;%m, thickness of slot
ri = rc-tr/2;%m, inner diameter of rail
ro = rc+tr/2;%m, outer diameter of rail
ris = rc-ts/2;%m, inner diameter of slot
ros = rc+ts/2;%m, outer diameter of slot

da = 0.4*intom;%m, depth of t-slot
g = 0.2*intom;%m, gap in t-slot
tslider = 0.1*intom;%m, thickness of slider
rso = ri+da;%m, slot outer radius
rsi = rso-g;%m, slot inner diameter

if ls > 2*rc
    disp('Warning: Drivetrain length cannot exceed center rail loop diameter.')
end

psi = asin(ls/2/rc);%rad, angle 

% Body coordinates of drivetrain relative to back bottom edge midpoint
hB = [0;0];
heB = [0;lw1];
hbB = [0;lw2];
haB = [0;ls];
hgB = [-w/2;lw1];
hfB = [w/2;lw1];
hdB = [-w/2;lw2];
hcB = [w/2;lw2];

%% Transform from body into "circle" coordinates
hC = B2C(hB,psi,rc);
heC = B2C(heB,psi,rc);
hbC = B2C(hbB,psi,rc);
haC = B2C(haB,psi,rc);
hgC = B2C(hgB,psi,rc);
hfC = B2C(hfB,psi,rc);
hdC = B2C(hdB,psi,rc);
hcC = B2C(hcB,psi,rc);

%% Calculate Relevant Values
% Calc for gap between wheels and inner edge of rail and outer edge of
% slider slot
gapc = ros - norm(hcC);
gapd = ri - norm(hdC);
gapg = ri - norm(hgC);
gapf = ros - norm(hfC);

% Calc for max length of slider
lslider = 2*sqrt(rso^2-(rsi+tslider)^2);

%% Calculate min radius for driving on inside and outside curves
%For inside curve, calculate radius of circle touching front and back edge of drivetrain and
%one wheel
%Assume: contact with wheel is at point on wheel furthest from drivetrain

%For outside curve, circle touches both wheels and the midpoint of the
%drivetrain
romin = ((la2-la1)^2/4+ha^2-rw^2)/2/(rw-ha);
%% Plot
icircle = circle(0,0,ri);
ccircle = circle(0,0,rc);
ocircle = circle(0,0,ro);
iscircle = circle(0,0,ris);
oscircle = circle(0,0,ros);

disp('Overhang distance for each wheel (safe if negative) [in]: ')
disp([gapg gapd gapc gapf])
disp('Min rail outer diameter for drivetrain dimensions [in]:')
disp(romin)
disp('Max length of slider [in]: ')
disp(lslider)

% figure(1)
% hold on
% plot([hB(1) haB(1)],[hB(2) haB(2)],'r-')
% plot([hgB(1) hfB(1)],[hgB(2) hfB(2)],'r-')
% plot([hdB(1) hcB(1)],[hdB(2) hcB(2)],'r-')
% ylim([-1,ls+1])
% axis('equal')
% title('Drivetrain Skeleton')
% hold off

figure(2)
hold on
plot([hgC(1) hfC(1) hdC(1) hcC(1)],[hgC(2) hfC(2) hdC(2) hcC(2)],'ro')
plot([hC(1) haC(1)],[hC(2) haC(2)],'r^')
plot(icircle(:,1),icircle(:,2),'k-','HandleVisibility','off')
plot(ccircle(:,1),ccircle(:,2),'k-','HandleVisibility','off')
plot(ocircle(:,1),ocircle(:,2),'k-','HandleVisibility','off')
plot(iscircle(:,1),iscircle(:,2),'k-','HandleVisibility','off')
plot(oscircle(:,1),oscircle(:,2),'k-','HandleVisibility','off')
plot([hC(1) haC(1)],[hC(2) haC(2)],'r-','HandleVisibility','off')
plot([hgC(1) hfC(1)],[hgC(2) hfC(2)],'r-','HandleVisibility','off')
plot([hdC(1) hcC(1)],[hdC(2) hcC(2)],'r-','HandleVisibility','off')
plot(0,0,'k.','HandleVisibility','off')
legend('Wheel','Slider','location','southeast')
ylim([-(ro+.01),ro+.01])
axis('equal')
xlabel('[m]')
ylabel('[m]')
hold off

function vC = B2C(vB,psi,rc) %Body to Circle transformation
    Tb2c = [cos(psi) -sin(psi); sin(psi) cos(psi)];
    trans = [rc;0];
    vC = Tb2c*vB + trans;
end

function coord = circle(x,y,r)
    th = 0:pi/50:2*pi;
    coord = [r * cos(th) + x;r * sin(th) + y]';
end

function [xCenter, yCenter, radius, a] = circlefit(x, y)
    % Found here: https://matlab.fandom.com/wiki/FAQ#How_can_I_fit_a_circle_to_a_set_of_XY_data.3F
    % circlefit(): Fits a circle through a set of points in the x - y plane.
    % USAGE :
    % [xCenter, yCenter, radius, a] = circlefit(X, Y)
    % The output is the center point (xCenter, yCenter) and the radius of the fitted circle.
    % "a" is an optional output vector describing the coefficients in the circle's equation:
    %     x ^ 2 + y ^ 2 + a(1) * x + a(2) * y + a(3) = 0
    % by Bucher Izhak 25 - Oct - 1991
    
    numPoints = numel(x);
    xx = x .* x;
    yy = y .* y;
    xy = x .* y;
    A = [sum(x),  sum(y),  numPoints;
         sum(xy), sum(yy), sum(y);
         sum(xx), sum(xy), sum(x)];
    B = [-sum(xx + yy) ;
         -sum(xx .* y + yy .* y);
         -sum(xx .* x + xy .* y)];
    a = A \ B;
    xCenter = -.5 * a(1);
    yCenter = -.5 * a(2);
    radius  =  sqrt((a(1) ^ 2 + a(2) ^ 2) / 4 - a(3));
end
