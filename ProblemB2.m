clear all
clc

%Define our geometry
geometry = @circleg ;
hhmax=[1/5, 1/20, 1/40, 1/60];
alpha = 4;
delta1 = 0.01;
fn = @(x) x - 1./(x+alpha);
EnE=zeros(1,length(hhmax));

for i = 1:length(hhmax)
    hmax=hhmax(i);
    [p, e, t] = initmesh(geometry, 'hmax', hmax);
    %[A,M]=assema(p,t,delta1,1,0);        % Calculating Stiffness and Mass matrix using built-in function
    M = MassAssembler2D(p,t);             % Alternatively, using self-defined function
    A = delta1*StiffnessAssembler2D(p,t);
    v_old = 1 + 20*rand([size(p,2),1]);
    v_new = v_old;
    dt = 0.01;                           % delta t
    
    for it = 1:200
        v_error = 10;
        while (v_error > 0.01)
            v_tmp = v_new;
            B = LoadAssembler2D(p,t,fn,v_tmp);
            v_new = (2*M + dt*A +dt*M)\(2*dt*M*B + (2*M - dt*A - dt*M)*v_old);
            v_error=norm(v_new-v_tmp);
            if isnan(v_error)
                break
            end
        end
        if isnan(v_error)
            break
        end
        v_old = v_new;
        prt(i,it) = PopulationIntegration(p,t,v_old);
    end
    tl(i) = it;
    i
end

figure
pdesurf(p,t,v_old)
view(-75,40)
set(gca,'BoxStyle','full','Box','off')
daspect([1 1 40])

figure
for i = 1:length(hhmax)
    plot(linspace(0,tl(i)*dt,tl(i)-1), prt(i,1:tl(i)-1),'-','LineWidth',1.5); hold on
end
legend({'h=1/4','h=1/20','h=1/40','h=1/60'}, 'location', 'northeast', 'Fontsize', 14)
ax = gca;
ax.FontSize = 13;
xlabel('Time')
ylabel('Population')

function F = PopulationIntegration(p,t,v)
np = size(p,2);
nt = size(t,2);
F = 0;
for K = 1:nt
    loc2glb = t(1:3,K);
    x = p(1,loc2glb);               % node x-coordinates
    y = p(2,loc2glb);               % node y-
    area = polyarea(x,y);
    tmp = area*sum(v(loc2glb))/3;
    F = F + tmp;
end
end

%%%MassAssembler2D
function M = MassAssembler2D(p,t)
np = size(p,2);                     % number of nodes
nt = size(t,2);                     % number of elements
M = sparse(np,np);                  % allocate mass matrix
for K = 1:nt                        % loop over elements
    loc2glb = t(1:3,K);             % local-to-global map
    x = p(1,loc2glb);               % node x-coordinates
    y = p(2,loc2glb);               % y
    area = polyarea(x,y);           % triangle area
    MK = [2 1 1;
        1 2 1;
        1 1 2]/12*area;             % element mass matrix
    M(loc2glb,loc2glb) = ...
        M(loc2glb,loc2glb) + MK;    % add element masses to M
end
end

%%% StiffnessAssembler2D
function A = StiffnessAssembler2D(p,t)
np = size(p,2);
nt = size(t,2);
A = sparse(np,np);                  % allocate stiffness matrix
for K = 1:nt
    loc2glb = t(1:3,K);             % local-to-global map
    x = p(1,loc2glb);               % node x-coordinates
    y = p(2,loc2glb);               % node y-
    [area,b,c] = HatGradients(x,y);
    AK =(b*b'+c*c')*area;           % element stiffness matrix
    A(loc2glb,loc2glb) = ...
        A(loc2glb,loc2glb)+ AK;     % add element stiffnesses to A
end
end

%Computation of aplha and hat-gradients
function [area,b,c] = HatGradients(x,y)
area=polyarea(x,y);
b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end

%%$ Load Assembler2D
function b = LoadAssembler2D(p,t,fn,v)
np = size(p,2);
nt = size(t,2);
b = zeros(np,1);
for K = 1:nt
    loc2glb = t(1:3,K);
    x = p(1,loc2glb);               % node x-coordinates
    y = p(2,loc2glb);               % node y-
    area = polyarea(x,y);
    xx = v(loc2glb);
    ft = fn(xx);
    bK = [ft(1);...
          ft(2);...
          ft(3)]/3*area;            % element load vector 
    b(loc2glb) = b(loc2glb) + bK;   % add element loads to b
end
end
