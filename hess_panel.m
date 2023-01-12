function [Cl,Cd,Cm,Strength] = hess_panel(point,N,V_inf,AoA,c)
n = length(point);  %Length of Given Station Points

AoA = deg2rad(AoA);

% To Find Upper and Lower Surfaces
Upper = point(1:round(n/2),:);

Bottom = point(round(n/2):n,:);

% Half Cosine x-Spacing for Each Surfaces To Obtain More Panels near LE and TE
beta = linspace(0,pi,(N+1)/2+1)';
x_bottom(:,1) = (0.5*(1 + cos(beta))); % Half cosine based spacing
x_upper(:,1) = flip(x_bottom(:,1));

% % Interpolant Curve Fitting for Each Surfaces
% cs_alt=csape(Bottom(:,1),Bottom(:,2));
% cs_ust=csape(Upper(:,1),Upper(:,2));
% 
% z_bottom = fnval(cs_alt,x_bottom);
%z_upper = fnval(cs_ust,x_upper);

%To Find z for Each Surfaces
p = polyfit(Upper(:,1),Upper(:,2),8);
z_upper = polyval(p,x_upper);

p = polyfit(Bottom(:,1),Bottom(:,2),8);
z_bottom = polyval(p,x_bottom);

% [fitresult] = createFit(Upper(:,1), Upper(:,2));
% p = coeffvalues(fitresult);
% z_upper = polyval(p, x_upper);
% 
% [fitresult] = createFit(Bottom(:,1), Bottom(:,2));
% p = coeffvalues(fitresult);
% z_bottom = polyval(p, x_bottom);

x = [x_bottom; x_upper];
z = [z_bottom; z_upper];

Coordinate = [x, z];   %To obtain panel coordinates clockwise

Coordinate(round((N+1)/2),:)=[];
Coordinate = Coordinate;

% Eq. 3.2.16
for i = 1:N
    for j=1:2
        Control(i,j) = 0.5 * ( Coordinate(i,j) + Coordinate(i+1,j) );
        d(i,j) = (Coordinate(i+1,j)) - (Coordinate(i,j));
    end
    theta(i,1) = atan2(d(i,2),d(i,1));
    ds(i,1) = sqrt(d(i,1)^2 + d(i,2)^2);
end

% Eq. 3.2.16
for i = 1:N
    for j = 1:N
        r(i,j+1) = sqrt((Control(i,1) - Coordinate(j+1,1))^2 ...
            + (Control(i,2) - Coordinate(j+1,2))^2 );
        
        r(i,j) = sqrt((Control(i,1) - Coordinate(j,1))^2 ...
            + (Control(i,2) - Coordinate(j,2))^2 );
        
        %To find inclination to the x-axis with trigonometric logic
        if i == j
            beta(i,j) = pi;
            
        elseif j== N-i+1
            
            beta(i,j) = (atan((Control(i,2) - Coordinate(j+1,2))/...
                (Control(i,1) - Coordinate(j+1,1)))...
                - atan((Control(i,2) - Coordinate(j,2))/...
                (Control(i,1) - Coordinate(j,1)))) -pi;
        else
            beta(i,j) = (atan((Control(i,2) - Coordinate(j+1,2))/...
                (Control(i,1) - Coordinate(j+1,1)))...
                - atan((Control(i,2) - Coordinate(j,2))/...
                (Control(i,1) - Coordinate(j,1)))) ;
        end
    end
end

% To Obtain Influence Coeffs
% Eqs. 3.2.13 - 3.2.14
for i = 1:N
    for j = 1:N
        if i == j
            An(i,j) = 0.5;
            At(i,j) = 0;
        else
            An(i,j) = 1/(2*pi) * ( sin(theta(i) - theta(j)) * log(r(i,j+1)/r(i,j)) + ...
                cos(theta(i)-theta(j)) * beta(i,j));
            
            At(i,j) = 1/(2*pi) * ...
                ((sin(theta(i)-theta(j)) * beta(i,j)) - ...
                (cos(theta(i) - theta(j)) * log( r(i,j+1) / r(i,j))));
        end
    end
end

Bn = -At;   Bt = An;        % Eq. 3.2.15

% To Solve Algebraic Linear Equations
A = zeros(N+1,N+1);                         %Eq. 3.2.20
A(1:N, 1:N) = An;                           %Eq. 3.2.21a
A(1:N, N+1) = sum(Bn,2);                    %Eq. 3.2.21b
A(N+1,1:N) = At(1,:) + At(N,:);             %Eq. 3.2.23a
A(N+1, N+1) = sum((Bt(1,:)+Bt(N,:)) ,2);    %Eq. 3.2.23b

for i = 1:N
    b(i,1) = - V_inf * sin(AoA - theta(i)); %Eq. 3.2.24a
end
%Eq 3.2.24b
b(N+1,1) = -V_inf * cos(AoA - theta(1)) - V_inf * cos(AoA - theta(N));

Strength = linsolve(A,b);  %To solve linear system

% Eq. 3.2.12b
At2 = [At, sum(Bt,2)];
for i = 1:N
    sum1 = 0;
    for j = 1:N+1
        sum1 = At2(i,j) * Strength(j) + sum1;
    end
    
    Vt(i,1) = sum1 + V_inf * cos(AoA - theta(i));   %Eq.3.2.12b
    
    cp(i,1) = 1-(Vt(i)/V_inf)^2;  % calculation of pressure coefficient
end

% Calculation of cl, cd, cm
for i = 1:N
    F(i,:) = (1/c) * (1 - cp(i)) * [-sin(theta(i)); cos(theta(i))] * ds(i);
    cl(i,1) = F(i,2) * cos(AoA) - F(i,1) * sin(AoA);
    cd(i) = F(i,1) * cos(AoA) + F(i,2) * sin(AoA);
    cm(i) = (1-cp(i)) * (-cos(theta(i)) * ds(i) * ...
        (Control(i,1) - c/4*cos(AoA)) -sin(theta(i))* ds(i) *...
        (Control(i,2) + c/4*sin(AoA)));
end

Cl = sum(cl);   Cd = sum(cd);   Cm = sum(cm);

    function [fitresult, gof] = createFit(x, y)
        % Fit: '4th Degree Polynomial'.
        [xData, yData] = prepareCurveData( x, y );
        
        % Set up fittype and options.
        ft = fittype( 'poly4' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = 'LAR';
        
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
    end

end

