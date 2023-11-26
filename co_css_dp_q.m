clear all; close all; clc;

a = 2.5; b = 0.05; alpha = 2; 
% q = 0.5; 
r1 = -1.3996; r2 = -1.5;      %r1 = tau'(f*), r2= tau''(f*) (host trade-off slope & curvature)
r11 = 0.9756; r22 = -1.5;     % r11 = beta'(gamma*), r22 = beta''(gamma*) (parasite trade-off slope & curvature)

q1 = linspace(0.1,1.3,50);
 for i = 1:length(q1)
     q = q1(i);
     g(i) = q;

syms f fm gamma gammam

tau = @(f) (1 - (r1.^2./r2).*(1-exp(r2.*(f - 0.5)./r1)));  %Host trade-off(f* = 0.5, tau* = 1)
beta = @(gamma) (2 - (r11.^2./r22).*(1-exp(r22.*(gamma - 1)./r11)));  %parasite trade-off (gamma* = 1, beta* = 2)

X = @(f,gamma)(b+alpha+gamma-tau(f))./beta(gamma); %equilibrium points as functions of f and gamma inlcuding trade-offs.
Y =  @(f,gamma)(-(1+f).*q.*beta(gamma).*(b+alpha+ gamma -tau(f))+ (beta(gamma)).^2.*(-b+a.*f-alpha+tau(f))+ sqrt(((beta(gamma)).^2.*(-4.*f.*q.*(-a.*beta(gamma)+ ...
    b.*(q + beta(gamma)) + q.*(alpha+gamma-tau(f))).*(b+ alpha +gamma -tau(f))+ (b.*(q + f.*q+ beta(gamma))+(1+f).*q.*(alpha+gamma-tau(f))-beta(gamma).*(a.*f...
    -alpha+tau(f))).^2))))./(2.*f.*q.*(beta(gamma)).^2);

T =  @(f,fm,gamma)(a-q.*(X(f,gamma) +Y(f,gamma))-b-beta(gamma).*Y(f,gamma))-(alpha+b+gamma-tau(fm));  % trace of mutant Jacobian matrix for fm=f
Det = @(f,fm,gamma) (a-q.*(X(f,gamma)+Y(f,gamma))-b-beta(gamma).*Y(f,gamma)).*(tau(fm)-alpha-b-gamma) - beta(gamma).*Y(f,gamma).*(a.*fm - q.*fm.*(X(f,gamma)+Y(f,gamma))+ gamma);

s = @(f,fm,gamma)(T(f,fm,gamma) + sqrt(T(f,fm,gamma).^2 - 4.*Det(f,fm,gamma)))./2; %host fitness proxy
r = @(f, gamma,gammam) beta(gammam).*X(f,gamma) - (alpha + b + gammam - tau(f));   %parasite fitness proxy

s1 = diff(s(f,fm,gamma),fm);                %derivative of s wrt fm
p1 = diff(r(f,gamma,gammam),gammam);        %derivative of r wrt gammam

s11 = diff(s1,fm); p11 = diff(p1,gammam);
EH = matlabFunction(subs(s11,fm,f));   EP = matlabFunction(subs(p11,gammam,gamma));

s12 = diff(s1,f);   p12 = diff(p1,gamma);
MH = matlabFunction(subs(s12,fm,f));   MP = matlabFunction(subs(p12,gammam,gamma));

s22 = diff(s1,gamma); p22 = diff(p1,f);
AH = matlabFunction(subs(s22,fm,f));   AP = matlabFunction(subs(p22,gammam,gamma));
condition = @(f,gamma)(EH(f,gamma) + MH(f,gamma)).*(EP(f,gamma) + MP(f,gamma)) - AH(f,gamma).*AP(f,gamma); %need condtion>0 for absolute convergence

gradH = subs(s1,fm,f);                    %substituting fm = f in s1 (host selection gradient)
gradP = subs(p1,gammam,gamma);            %substituting gammam = gamma in s1 (parasite selection gradient)

[xsol, ysol] = vpasolve([gradH, gradP], [f, gamma],[0.5, 1]);
% xSol = sol.f;      ySol = sol.gamma;
fs(i) = xsol;     gs(i) = ysol;
P = Y(fs,gs)./(X(fs,gs)+Y(fs,gs));
% EH(fs,gs);   EP(fs,gs);
% MH(fs,gs);   MP(fs,gs);
% CS(fs,gs);
% R = beta(gamma).*(a-b)./(q.*(alpha+b+gamma-1))
 end
%%
figure(1)
plot(g, fs,'-b','Linewidth',4,'MarkerSize',4)
hold on
plot(g, gs,'--k','Linewidth',5,'MarkerSize',4)
legend('$f$','$\gamma$','interpreter','Latex','Location','north','FontSize',14)
set(gca,'FontSize',14)
xlabel('$q$','interpreter','Latex','Fontsize',25)
ylabel('$co-CSS$','interpreter','Latex','Fontsize',23)
H=gca; H.LineWidth = 2; box off;
xlim([0.01 1.3])
% ylim([0 2])
%%
figure(2)
plot(g, P,'r','Linewidth',5,'MarkerSize',4)
% legend('$f$','$\gamma$','interpreter','Latex','Location','northwest','FontSize',14)
set(gca,'FontSize',14)
xlabel('$q$','interpreter','Latex','Fontsize',25)
ylabel('$Disease~prevalence (P)$','interpreter','Latex','Fontsize',23)
H=gca; H.LineWidth = 2; box off;
xlim([0.01 1.3])
% ylim([0.3 2])
%%
% figure(2)
% plot(g, EH(fs,gs),'r','Linewidth',5,'MarkerSize',4)
% hold on 
% plot(g, EP(fs,gs),'-.b','Linewidth',5,'MarkerSize',4)
% legend('EH','EP')
% set(gca,'FontSize',14)
% xlabel('$q$','interpreter','Latex','Fontsize',25)
% ylabel('$ES$','interpreter','Latex','Fontsize',23)
% H=gca; H.LineWidth = 2; box off;
% 
% figure(3)
% plot(g, condition(fs,gs),'r','Linewidth',5,'MarkerSize',4)
% set(gca,'FontSize',15)
% xlabel('$q$','interpreter','Latex','Fontsize',25)
% ylabel('$condition$','interpreter','Latex','Fontsize',23)
% H=gca; H.LineWidth = 2; box off;
% 
% figure(4)
% plot(g, X(fs,gs),'r','Linewidth',5,'MarkerSize',4)
% hold on 
% plot(g, Y(fs,gs),'-.b','Linewidth',5,'MarkerSize',4)
% legend('X','Y')
% set(gca,'FontSize',14)
% xlabel('$q$','interpreter','Latex','Fontsize',25)
% ylabel('$Population~densities$','interpreter','Latex','Fontsize',23)
% H=gca; H.LineWidth = 2; box off;
