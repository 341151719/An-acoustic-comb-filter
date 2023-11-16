close all 
clear all
clc

flag=1;
BL2=4.6^2; 
freq=40*10.^linspace(0,log10(800),100000);

c0=343;
Wd=0.08;
rcA=1.25*c0*Wd^2;dep=0.15;k0=1000;
freq=freq';w=2*pi*freq;s=i*w;k=w/c0;
mass=4e-3;stiff=(2*pi*250)^2*mass;d=rcA*0.3;
Rc=0.1;Lc=0.15e-3;
fe = input('Please enter an array of frequencies to be perfectly matched, enclosed in brackets. The format should be [f1, f2, f3, f4]:');
Nfe = length(fe);
if Nfe <= 0
    error('The number of frequencies entered must be greater than 0');
end

Noct=sqrt(2);
fe0=80;
Lp=100e-3*ones(size([0:Nfe-1]));

[Rp,Cp]=pd_afa_comb(Lp+Lc,fe,mass,dep,stiff,d,BL2,rcA);
Rp=Rp-Rc;if min(Rp)<0, disp(['Rp=',num2str(min(Rp),3),'ohm <0']),return;end

for n=1:Nfe
    Ap(:,n)=1./(Rp(n)+Lp(n).*s+1./(Cp(n).*s));%admitance
end
Zc=Rc+Lc.*s;Zp=1./sum(Ap,2);
Ze=Zc+Zp;
    Zm=(mass*s+d+stiff./s)/rcA;
    dZ=BL2./Ze./rcA;
    Z=Zm+dZ;    
    beta=abs((1-Z)./(1+Z)).^2; betam=abs((1-Zm)./(1+Zm)).^2;
    alfa=1-beta; alfam=1-betam;
% %
set(gcf,'pos',[200 20 1200 800])
axes('pos',[0.12 0.71 0.78 0.25])
  h1= plot(freq,real(Z),'-','linew',1.2); hold 
   h2= plot(freq,real(Zm),'--','linew',1.2);
   plot(freq,ones(size(freq)),'-.','color',[0.5 0.5 0.5],'linew',.5);
    ylabel('${D/\rho_0 c_0}$','Fontsize',9,'interpreter','latex')
    ymin=0.2;ymax=60;
   for n=1:Nfe
       plot(fe(n)*ones(100),linspace(ymin,ymax,100),'-.k','linew',0.5,'color',[0.5 0.5 0.5])
   end

         set(gca,'xscale','log','yscale','log','fontname','times','fontsize',9,...
         'xlim',[60 900],'ylim',[ymin ymax],'ytick',[1 5 20 50 ],'yticklabel',[1 5 20 50 ],...
         'xtick',round(fe),'xticklabel',[]);
%      subfig(0,1,1)
     
axes('pos',[0.12 0.44 0.78 0.25])
    h1=plot(freq,imag(Z),'-','linew',1.2);hold 
    h2=plot(freq,imag(Zm),'--','linew',1.2);
    plot(freq,zeros(size(freq)),'-.','color',[0.5 0.5 0.5],'linew',0.5);
    ylabel('${\chi/\rho_0 c_0}$','Fontsize',9,'interpreter','latex')
    ymin=-36;ymax=25;
   for n=1:Nfe
       plot(fe(n)*ones(100),linspace(ymin,ymax,100),'-.k','linew',0.5,'color',[0.5 0.5 0.5])
   end
         set(gca,'xlim',[60 900],'xscale','log','fontname','times','fontsize',9,...
         'ylim',[ymin ymax],'ytick',[-30:15:30],'yticklabel',[-30:15:30],...
         'xtick',round(fe),'xticklabel',[]);
             temp=legend([h1 h2],'With shunts','Without shunts');
        set(temp,'orientation','hori','box','off','pos',[0.5 0.36 0.1 0.1]);
% subfig(0,1,2)
               
     axes('pos',[0.56 0.48 0.26 0.08])
         h1=plot(freq,imag(Z),'-','linew',1.2);hold 
         h2=plot(freq,imag(Zm),'--','linew',1.2);
         plot(freq,zeros(size(freq)),'-.','color',[0.5 0.5 0.5],'linew',0.5);
          set(gca,'xlim',[200 350],'xscale','log','fontname','times','fontsize',7,...
         'ylim',[-1 2],'ytick',[0 2],'yticklabel',[0 2],...
         'xtick',[200 226 270 320],'xticklabel',[200 226 270 320]);
     set(gca,'color','w')
   for n=1:Nfe
       plot(fe(n)*ones(100),linspace(ymin,ymax,100),'-.k','linew',0.5,'color',[0.5 0.5 0.5])
   end
          
axes('pos',[0.12 0.09 0.78 0.29])
    semilogx(freq,alfa,'-','linew',1.2);hold 
    semilogx(freq,alfam,'--','linew',1.2)
        set(gca,'xlim',[60 900],'xscale','log','xtick',...
         freq(1)*2.^[0:log2(freq(end)/freq(1))],'fontname','times','fontsize',9,...
         'xlim',[60 900],'ylim',[0 1],'ytick',[0:0.25:1],'yticklabel',[0:0.25:1],...
         'xtick',round(fe),'xticklabel',round(fe));
    for n=1:Nfe
       plot(fe(n)*ones(100),linspace(0,1,100),'-.k','linew',0.5,'color',[0.5 0.5 0.5])
   end
         ylabel('\it{\alpha}','fontname','times','fontsize',9)
         xlabel('Frequency(Hz)','fontname','times','fontsize',9)
         
         yyaxis right
         plot(fe,Rc+Rp,'ko','markersize',3,'markerfacecolor','k');
         set(gca,'yscale','log','ylim',[(Rc+min(Rp))/2 (Rc+max(Rp))*2],...
             'ytick',[0.1 0.3 1 3 10])
         ax=gca;ax.YColor = 'k';
         ylabel('${R_n(\Omega)}$','fontname','times','fontsize',9,'interpreter','latex')
%       subfig(0,1,3)
      
set(gcf,'paperposition',[0 0 9 10]);
print -dtiff -r600 Fig5;
fprintf('\nFe: ');
disp(fe);
fprintf('\nCp: ');
disp(Cp);
fprintf('\nRp: ');
disp(Rp);
function [R,C]=pd_afa_comb(L,fr,m,dep,k0,d,BL2,rcA)

 wr=fr*2*pi;
 s=i*wr;
k=k0;
 Xm=imag(m*s+k./s);
 Xe=BL2*Xm./(Xm.^2+(rcA-d)^2);
 C=1./(wr.*L-Xe)./wr;
 R=BL2*(rcA-d)./(Xm.^2+(rcA-d)^2);
end





