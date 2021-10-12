clear all;close all;clc;
%a=2, b=0.5, q=0.2, gamma=0.2, alpha=0.3, beta=3)
%Plotting CSS points for varying values of gamma, for b=0.5,0.6 and 0.7.

b1 = linspace(0.01,0.5,200);
for i = 1:length(b1)
    a = 5;
    %b = 0.05;
    b = b1(i);
    gamma = 0.5;
    g(i) = 1./b;
    alpha = 4;
    q = 0.2;
    beta = 6;
    r1 = -0.23; %r'(tau*)
    r2 = -0.25; %r''(tau*)

r = @(tau)(1 - (r1.^2./r2).*(1-exp(r2.*(tau-1)./r1))); %r1 is the derivative of r. 
dr = @(tau)r1.*exp(r2.*(tau-1)./r1) ; %r1 is the derivative of r.. 
X = @(tau) (-b-alpha-gamma+tau)./(r(tau)-beta); %equilibrium points in terms of tau.
Y = @(tau) ((b+alpha+gamma-tau).*(b.*q + a.*r(tau)- b.*r(tau)+ q.*alpha- a.*beta + b.*beta+ q.*gamma- q.*tau))./((r(tau)-beta).*(b.*q-b.*r(tau)+ q.*alpha- r(tau).*alpha + b.*beta + alpha.*beta+ q.*gamma - q.*tau+r(tau).*tau-beta.*tau));
%s = @(tau,taum) (a-q.*(X(tau)+Y(tau))-b-(beta-r(taum)).*Y(tau)).*(alpha-taum+b+gamma)+gamma.*Y(tau).*(beta-r(taum)) is the fitness function;

%s1 = derivative of s w.r.t taum and then substitute taum=tau.

s1 = @(tau) dr(tau).*Y(tau).*(alpha+b-tau)-(a-b-q.*(X(tau)+Y(tau))-(beta-r(tau)).*Y(tau));
        if i==1
            guessval = 1;
            guessval2 = 3.5;
        else
            guessval = singular(i-1);
            guessval2 = singular2(i-1);
        end

singular(i) = fzero(s1,guessval);  %fzero(fun,x0) tries to find a zero of fun near x0, if x0 is a scalar and fun is a function handle.  

end
plot(g,singular,'k','Linewidth',2.5,'MarkerSize',4)
hold on
plot(g,r(singular),'r','Linewidth',2.5,'MarkerSize',4)
set(gca,'FontSize',14) 
legend('Tolerance','Resistance','Location','east','FontSize',15);
xlabel('$Life~span(1/b)$','interpreter','Latex','Fontsize',25)
ylabel('$CSS$ $investment$','interpreter','Latex','Fontsize',25)
H = gca;
H.LineWidth = 2;
set(gca,'box','off')
