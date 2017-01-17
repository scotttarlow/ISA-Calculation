%ISA Orientation Calculations
%requires inside_corner.m and outstide_corner.m
%Works only with square matricies, with equal grid spacing 

clear

%velocity field calculations 
%reduce = 10m0;
dx = 200;%m 
dy = 200; %m
x = 0:dx:3000; %m
y = 0:dy:3000; %m
slabdip = 90 * pi / 180; % 90 degree dip
dip = pi - slabdip; 
v = .005;% m/yr
[X,Y] = ndgrid(x,y); % X/Y grid 
UX = zeros(size(X)); 
UY = UX;


% Calculate velocity field inside and outside corner
u = zeros(length(x),length(y),2); % field velocities, 1 for x and 2 for y
u2 = zeros(length(x),length(y),2);
for i = 1:length(x)
    for j = 1: length(y)
        if (x(i)<0) && (y(j) <= abs(x(i))*tan(slabdip))
            [UX(i,j),UY(i,j)] = inside_corner(x(i),y(j),dip,v);
        else
            [UX(i,j),UY(i,j)] = outside_corner(x(i),y(j),dip,v);
        end
    end
end

% Convert velocities from m/yr to m/s
UX2 = UX * (3.17 * 10^-10);
UY2 = UY * (3.17 * 10^-10); 




%Derivatives (centered) for calculating Strain Rate Tensor 
%(dU(x + dx) - dU(x - dx))/(2*dx) 

%Strain Rate Tesnor 
%      |exx exy exz|
% E =  |eyx eyy eyz|
%      |ezx ezy ezz|

%exx = dUX/dX exy = .5*(dUX/dY + dUY/dX) exz = 0
%eyx = .5*(dUY/dX + dUX/dY) eyy = dUY/UY eyz = 0
%ezx = 0 ezy = 0 ezz = 0


%Calculating strain rate tensor 

E = {};
for i = 2:length(x)-1
    for j = 2:length(y) - 1
        
      exx = (UX2(i+1,j) - UX2(i-1,j)) / (2*dx);
        
      exy = 0.5*(((UX2(i,j+1) - UX2(i,j-1)) / (2*dy)) ...
            + ((UY2(i+1,j) - UY2(i-1,j)) / (2*dx)));
        
      exz = 0;
        
      eyx = 0.5*(((UY2(i+1,j) - UY2(i-1,j)) / (2*dx)) ...
            + ((UX2(i,j+1) - UX2(i,j-1)) / (2*dy)));
        
      eyy = (UY2(i,j+1) - UY2(i,j-1)) / (2*dy);
        
      eyz = 0;
      ezx = 0;
      ezy = 0;
      ezz = 0;
        
      E{i,j} =  [exx exy exz ; eyx eyy eyz ; ezx ezy ezz];
    end
end 





%%%Calculating ISA Orientation (Kam/Ribe 2002, Appendix A)%%%%%%%% 

%Step 1, Calculate L, the Velocity Gradient Tensor
% L = E - (epsi_x * omega_x + epsi_y * omega_y + epsi_z * omega_z)
% epsi is the permulation symbol
% omega is the vorticity vector 


%calculate permutation symbol
%http://en.wikipedia.org/wiki/Levi-Civita_symbol#Three_dimensions
epsi = zeros(3,3,3);
epsi(1,2,3) = 1; epsi(2,3,1) = 1; epsi(3,1,2) = 1;
epsi(3,2,1) = -1; epsi(1,3,2) = -1; epsi(2,1,3) = -1; 


%calculate vorticity (centered derivatives)
omega = {};
for i = 2: length(x)-1
    for j = 2:length(y)-1
        Uyx = (UY2(i+1,j) - UY2(i-1,j)) / (2*dx);
        Uxy = (UX2(i,j+1) - UX2(i,j-1)) / (2*dy);
        omega{i,j} = [ 0 ; 0 ; Uyx - Uxy];
    end
end

%calculate L
L = {};
for i = 2: length(x)-1
    for j = 2:length(y)-1
       L{i,j} = E{i,j} - 0.5*(epsi(:,:,1)*omega{i,j}(1) + epsi(:,:,2)*omega{i,j}(2) + epsi(:,:,3)*omega{i,j}(3));
    end
end


%Step 2, Calculate F, defomration gradient tensor
% F = exp(L*t) = I + L*t + (L^2/2!)*t^2 + (L^3/3!)*t^3
%for t = tmax
%tao_max = 75 = tmax * edot
%edot is the absolute value of the largest eigenvalue of the strain
%rate tensor (Kam/ribe 2002)

%calculate edot and tmax
for i = 2:length(x)-1
    for j = 2:length(y)-1
        if max(max(isnan(E{i,j}))) == 0
            eigval = eig(E{i,j});
            if abs(eigval(1)) > abs(eigval(2)) && abs(eigval(1)) > abs(eigval(3))
                edot(i,j) = abs(eigval(1));
                tmax(i,j) = 75/edot(i,j);
            elseif abs(eigval(2)) > abs(eigval(1)) && abs(eigval(2)) > abs(eigval(3))
                edot(i,j) = abs(eigval(2));
                tmax(i,j) = 75/edot(i,j);
            else
                edot(i,j) = abs(eigval(3));
                tmax(i,j) = 75/edot(i,j);
            end
        else
            edot(i,j) = 0;
            tmax(i,j) = 0;
            
        end
    end
end

%Calculate deformation gradient tensor

F = {};
I = zeros(3,3); I(1,1) = 1; I(2,2) = 1; I(3,3) = 1;
for i = 2:length(x)-1
    for j = 2:length(y)-1
        F{i,j} = I;
        for k=1:4
            F{i,j} = F{i,j} + (((L{i,j}^k) * (tmax(i,j)^k)) / factorial(k));
            
            %(L{i,j}*tmax(i,j))^k/factorial(k); 
        end
    end
end

%step 3, Calculate U,  Left stretch Tensor 
%%% Kam/Ribe 2002 Appenix A Version %%%%% 
% U = F * F', as shown in Conrad 2007 ISA code. Kam/Ribe 2002 Appendix A 
%is incorrect

U = {};
for i = 2:length(x)-1
    for j = 2:length(y)-1 
        U{i,j} = F{i,j} * F{i,j}';
    end 
end


%step 4, Calculate ISA Vector
% ISA vector is the eigen vector corresponding to the greatest eigenvalue 
% of the left stretch tensor (U) as t -> infinity. 
ISAa = {};
for i = 2:length(x)-1
    for j = 2:length(y)-1
        if max(max(isnan(U{i,j}))) == 0
            [V,D] = eig(U{i,j});
            if sum(D(:,1)) > sum(D(:,2)) && sum(D(:,1)) > sum(D(:,3))
                ISAa{i,j} = V(:,1);
            elseif sum(D(:,2)) > sum(D(:,1)) && sum(D(:,2)) > sum(D(:,3))
                ISAa{i,j} = V(:,2);
            else
                ISAa{i,j} = V(:,3);
            end
        else
            ISAa{i,j} = [ 0 0 0];
         
        end    
    end 
end

%formatting for ISAa (plotting)
ISAx = zeros(length(x),length(y));
ISAy = zeros(length(x),length(y));

for i = 2:length(x)-1
    for j = 2:length(y)-1
        ISAx(i,j) = ISAa{i,j}(1);
        ISAy(i,j) = ISAa{i,j}(2);
    end
end

%%%%%%% Calculating The Lag Parameter %%%%%%%%%
%Lag = (1/edot)*D(theta)/DT, where D/DT is the material derivative
%Lag = (1/edot)*D(theta)/DT = (1/edot) * (d(theta)/dt + u dot grad(theta))
%d(theta)/dt = 0 because of instantaneous flow. 

%Calculating theta
%theta = acos(u dot ISA) and u is the local flow direction
%dUv is the divergence of the flow direction
Uv = cell(length(x),length(y));
Uv2 = Uv;
dUv = zeros(length(x),length(y));
theta = [];

for i = 2:length(x)-1
    for j = 2:length(y)-1
        Uv{i,j} = [UX2(i,j),UY2(i,j),0]/sqrt(UX2(i,j)^2 + UY2(i,j)^2);
        theta(i,j) = acos(dot(Uv{i,j},ISAa{i,j}));
    end
end


for i = 3:length(x)-2
    for j = 3:length(y)-2
        dtheta{i,j} = [(theta(i+1,j) - theta(i-1,j))/(2*dx),...
                       (theta(i,j+1) - theta(i,j-1))/(2*dy),...
                       0];
    end
end
lag = zeros(length(x),length(y));
for i = 3:length(x)-2
    for j = 3:length(y)-2
        lag(i,j) = abs(dtheta{i,j}(1)*UX2(i,j) + dtheta{i,j}(2)*UY2(i,j))/edot(i,j);
    end
end

        


%plotting 
h = figure();
quiver(X,Y,UX,UY);
set(gca,'Ydir','reverse');
hold on
h = quiver(X,Y,ISAx,ISAy);
set(h,'ShowArrowHead','off');
xlabel('m')
ylabel('m')
saveas(h,'figure3.pdf')





