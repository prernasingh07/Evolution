clear; clc; close all;  
% try
% pool = parpool(4);
% % delete(gcp('nocreate'))
% end

a = 2.5; b = 0.2; alpha = 2; q = 0.5; 
r1 = -1.59368; %r2 = -1;        %r1 = tau'(f*), r2= tau''(f*) (host trade-off slope & curvature)
r11 = 0.909091; %r22 = -2.8;     % r11 = beta'(gamma*), r22 = beta''(gamma*) (parasite trade-off slope & curvature) 

syms f gamma fm gammam r2 r22
tau = @(f,gamma,r2) (1 - (r1.^2./r2).*(1-exp(r2.*(f - 0.5)./r1)));  %Host trade-off(f* = 0.5, tau* = 1)
beta = @(f,gamma,r22) (2 - (r11.^2./r22).*(1-exp(r22.*(gamma - 1)./r11)));  %parasite trade-off (gamma* = 1, beta* = 2)

X = @(f,gamma,r2,r22)(b+alpha+gamma-tau(f,gamma,r2))./beta(f,gamma,r22); %equilibrium points as functions of f and gamma inlcuding trade-offs.
Y =  @(f,gamma,r2,r22)(-(1+f).*q.*beta(f,gamma,r22).*(b+alpha+ gamma -tau(f,gamma,r2))+ (beta(f,gamma,r22)).^2.*(-b+a.*f-alpha+tau(f,gamma,r2))+ sqrt(((beta(f,gamma,r22)).^2.*(-4.*f.*q.*(-a.*beta(f,gamma,r22)+ ...
    b.*(q + beta(f,gamma,r22)) + q.*(alpha+gamma-tau(f,gamma,r2))).*(b +alpha+gamma -tau(f,gamma,r2))+ (b.*(q + f.*q+ beta(f,gamma,r22))+(1+f).*q.*(alpha+gamma-tau(f,gamma,r2))-beta(f,gamma,r22).*(a.*f...
    -alpha+tau(f,gamma,r2))).^2))))./(2.*f.*q.*(beta(f,gamma,r22)).^2);

T =  @(f,fm,gamma,r2,r22)(a-q.*(X(f,gamma,r2,r22) +Y(f,gamma,r2,r22))-b-beta(fm,gamma,r22).*Y(f,gamma,r2,r22))-(alpha+b+gamma-tau(fm,gamma,r2));  % trace of mutant Jacobian matrix for fm=f
Det = @(f,fm,gamma,r2,r22)(a-q.*(X(f,gamma,r2,r22)+Y(f,gamma,r2,r22))-b-beta(fm,gamma,r22).*Y(f,gamma,r2,r22)).*(tau(fm,gamma,r2)-alpha-b-gamma) - beta(fm,gamma,r22).*Y(f,gamma,r2,r22).*(a.*fm - q.*fm.*(X(f,gamma,r2,r22)+Y(f,gamma,r2,r22))+ gamma);

s = @(f,fm,gamma,r2,r22)(T(f,fm,gamma,r2,r22) + sqrt(T(f,fm,gamma,r2,r22).^2 - 4.*Det(f,fm,gamma,r2,r22)))./2; %host fitness proxy
r = @(f, gamma,gammam,r2,r22) beta(f,gammam,r22).*X(f,gamma,r2,r22) - (alpha + b + gammam - tau(f,gammam,r2));   %parasite fitness proxy

s1 = diff(s(f,fm,gamma,r2,r22),fm);                %derivative of s wrt fm
p1 = diff(r(f,gamma,gammam,r2,r22),gammam);        %derivative of r wrt gammam

s11 = diff(s1,fm); p11 = diff(p1,gammam);
EH = matlabFunction(subs(s11,fm,f));   EP = matlabFunction(subs(p11,gammam,gamma));

s12 = diff(s1,f);   p12 = diff(p1,gamma);
MH = matlabFunction(subs(s12,fm,f));   MP = matlabFunction(subs(p12,gammam,gamma));

s22 = diff(s1,gamma); p22 = diff(p1,f);
AH = matlabFunction(subs(s22,fm,f));   AP = matlabFunction(subs(p22,gammam,gamma));
Deter = @(f,gamma,r2,r22)((EH(f,gamma,r2,r22) + MH(f,gamma,r2,r22)).*(EP(f,gamma,r2,r22) + MP(f,gamma,r2,r22))) - AH(f,gamma,r2,r22).*AP(f,gamma,r2,r22);
Trace = @(f,gamma,r2,r22) X(f,gamma,r2,r22).*(EH(f,gamma,r2,r22) + MH(f,gamma,r2,r22))+Y(f,gamma,r2,r22).*(EP(f,gamma,r2,r22) + MP(f,gamma,r2,r22));
% need Trace <0 and Det>0 for convergent stability
gradH = matlabFunction(subs(s1,fm,f));                   %substituting fm = f in s1 (host selection gradient)
gradP = matlabFunction(subs(p1,gammam,gamma));           %substituting gammam = gamma in s1 (parasite selection gradient)
%%
hostcurv = linspace(-1.5,1.5,50);
parcurv = linspace(-1.5,1.5,50);
hostcurv = (nonzeros(hostcurv))';
parcurv = (nonzeros(parcurv))';
% tic
 for ih = 1:length(hostcurv)
     for ip = 1:length(parcurv)
     
     r2 = hostcurv(ih);     r22 = parcurv(ip);

        [xsol, ysol] = vpasolve([gradH(f,gamma,r2,r22), gradP(f,gamma,r2,r22)], [f, gamma],[0.5, 1]);
        fs = double(xsol); gs = double(ysol);
    
    if (isempty(fs) || isempty(gs))
            xsol1(ih,ip) = 1; % no singularity occurs(system has no solution)
    else
        ESH = EH(fs,gs,r2,r22); ESP = EP(fs,gs,r2,r22);
        Xh = X(fs,gs,r2,r22); Yp = Y(fs,gs,r2,r22);
        MIH = MH(fs,gs,r2,r22); MIP = MP(fs,gs,r2,r22);
        IH = ESH+ MIH; %isoclinic stability of host
        IP = ESP+ MIP; %isoclinic stability of parasite
        D = Deter(fs,gs,r2,r22);  T = Trace(fs,gs,r2,r22);

        %%
      if Xh>0 && Yp>0 && xsol>0 && ysol>0 && xsol<=1 && ysol<=2
          if ESH>0 && D>0 && T<0 && MIH<0 && ESP+MIP<0 && ESP<0     %host BP  & Parasite CS  
            xsol = -200;
          elseif ESH>0 && ESH+MIH<0 && MIH<0 && ESP>0 && ESP+MIP>0  %host branch but parasite repeller 
            xsol = -60;
          elseif ESP>0 && ESP+MIP<0 && MIP<0 && ESH>0 && ESH+MIH<0 && MIH<0  %host branch and parasite branch 
            xsol = -120;
           elseif ESH<0 && D>0 && T<0 && ESP<0 && ESP+MIP<0 && ESH+MIH<0 % host and parasite CSS
            xsol = 200;
           elseif D<0 && ESH>0 && ESP+MIP<0 && ESP<0 && ESH+MIH>0   % host repeller & Parasite CS
            xsol = 500;
           elseif D<0 && ESP>0 && ESH+MIH<0 && ESH<0 && ESP+MIP>0   %parasite repeller & host CS
            xsol = 700;
           elseif D<0 && ESH+MIH>0 && ESP+MIP>0 && ESH>0 && ESP>0   %both host and parasite max/min
            xsol = 100;
%           elseif D>0 && T>=0 && ESH+MIH>0 && ESP+MIP>0 && ESH>0 && ESP>0 %cycles 
%           xsol = 420;
%          elseif ESH<0 && ESH+MIH>0 && ESP<0 && ESP+MIP>0 
%           xsol = 350;
   
          else
             xsol = 420;
          end        
     else
        xsol = 1;
     end   
        xsol1(ih,ip) = xsol;
    end
    end
%      toc
 end
%%

figure;
set(gcf, 'Position');
imagesc(hostcurv,parcurv,xsol1'); shading flat;

colormap('jet')
% colorbar
set(gca,'YDir','normal'); set(gca,'FontSize',20)
xlabel('$Host~trade-off~curvature$','interpreter','Latex','Fontsize',21)
ylabel('$Parasite~trade-off~curvature$','interpreter','Latex','Fontsize',21)
sgtitle('$b = 0.2$','interpreter','Latex','Fontsize',22)
%%
% de lete(pool)
