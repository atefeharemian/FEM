clear
close all
clc;


%Define parameters
a = -1;     % xl
b = 1;      % xr
N = 13;     % Number of Points
h = (b-a)/(N-1);       % Mesh-size;
d=0.01;
x = a:h:b;            % Node coordinates
TOL = 10^(-9);
lambda = 0.9;
eta2=ones(N-1,1);
ER=[];               %Array to collect errors
Store_N=[];          %Array to collect number of nodes

 while sum(eta2) >= TOL && length(x) <= 10^4
     A=my_stiffness_matrix_assembler(x);
     B=my_load_vector_assembler(x);
     M=my_Mass_matrix_assembler(x);
     zi = A\B/d;            
     yi=-M\(A*zi);          
     N = length(x);
     eta2 = zeros(N,1);
     for i=1:length(x)-1
        h=x(i+1)-x(i);
        eta2(i)=1/2*h^3*(f(x(i))^2+f(x(i+1))^2);  
     end
    ER(end+1)=sum(eta2);
    Store_N(end+1)=length(x);
    for i = 1: length(eta2)          % Refine x
       if eta2(i) > lambda*max(eta2)
           x = [x (x(i+1)+x(i))/2];
       end
    end
    x = sort(x);
 end
 
 A=my_stiffness_matrix_assembler(x);
 B=my_load_vector_assembler(x);
 M=my_Mass_matrix_assembler(x);
 zi = A\B/d; 
 yi=-M\(A*zi);
 for i=1:length(x)-1            %Final eta2
     h=x(i+1)-x(i);
     eta2(i)=1/2*h^3*(f(x(i))^2+f(x(i+1))^2);
 end

 inverse_store_N=zeros(length(Store_N),1);
 for i=1:length(Store_N)
   inverse_store_N(i)=1/Store_N(i);
 end


%%
figure
 plot(x,zi, 'g--')
 xlabel("x")
 ylabel("u_h")
 title("uh")
 ax = gca;
 ax.FontSize = 8;
 set(gca,'DefaultTextFontSize',8);



figure
 scatter(x(2:end),eta2, 'green');
 xlabel("x")
 ylabel("\eta2 ")
 title("Error for each element")
 ax = gca;
 ax.FontSize = 7;
 set(gca,'DefaultTextFontSize',7);


 figure
 plot(x(2:end),(1./diff(x)),'g--')
 xlabel("x")
 title("Grid size distribution")
 ax = gca;
 ax.FontSize = 8;
 set(gca,'DefaultTextFontSize',8);

 figure
 fx=arrayfun(@(x) f(x), x);
 R_uh=fx+(d*yi)';
 plot(x,R_uh, 'g--')
 xlabel("x")
 ylabel("R_uh")
 title("Residual")
 ax = gca;
 ax.FontSize = 8;
 set(gca,'DefaultTextFontSize',8);

figure
plot(x(2:end), sqrt(eta2), 'g--')
xlabel("x")
ylabel("\eta ")
title("eta")
ax = gca;
ax.FontSize = 8;
set(gca,'DefaultTextFontSize',8);


figure
loglog (Store_N, ER, 'r--');hold on
loglog (Store_N,  inverse_store_N, 'b--');hold on
legend(["etasum", "N"])
legend("\eta sum", "N",'Location','southwest')
xlabel("Nmber of nodes")
ax = gca;
ax.FontSize = 8;
set(gca,'DefaultTextFontSize',8);

figure
fx=arrayfun(@(x) f(x), x);
plot(x, zi, 'g.-');hold on
plot(x, fx, 'k.-');hold on
legend(["R_uh","f"])
ax = gca;
ax.FontSize = 8;
set(gca,'DefaultTextFontSize',8);

 %%
function A=my_stiffness_matrix_assembler(x)


N = length(x) - 1;                       % number of elements
A = zeros(N+1, N+1);                     % initialize stiffnes matrix to zero
for i = 1:N                              % loop over elements
    h = x(i+1) - x(i);                     % element length
    n = [i i+1];                           % nodes
    A(n,n) = A(n,n) + [1 -1; -1 1]/h;      % assemble element stiffness
end

A(1,1) = 1.e+6;                          % adjust for BC
A(N+1,N+1) = 1.e+6;
end
%%
function gx = g(x)
      gx=10*x*sin(7*pi*x);
  end

  function fx = f(x)
      if g(x)>abs(x)
          fx=abs(x);
      elseif g(x)<-abs(x)
          fx=-abs(x);
      else
          fx=g(x);
      end
  end
%%
function B = my_load_vector_assembler(x)

N = length(x) - 1;
B = zeros(N+1, 1);
for i = 1:N
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [f(x(i)); f(x(i+1))]*h/2;
end
end
%%
function zi = my_first_fem_solver(x,d)
A=my_stiffness_matrix_assembler(x);
B=my_load_vector_assembler(x);
zi = d.*A\B;            % solve system of equations

end
%%
%Mass-Matrix Assembler
function M = my_Mass_matrix_assembler(x)
 N=length(x)-1;
       M=zeros(N+1,N+1);
       for i=1:N
           h=x(i+1)-x(i);
           n = [i i+1];
           M(n,n)= M(n,n)+[1/3 1/6; 1/6 1/3]*h;
       end
end      

%%
function yi = my_first_femlaplacian_solver(x,zi)
A=my_stiffness_matrix_assembler(x);
M = my_Mass_matrix_assembler(x);
yi = -M\(A*zi);
end
