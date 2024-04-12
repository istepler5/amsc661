function driver()
fsz = 16; % Fontsize
nt = 20;
N = 20;
tol = 1e-6;
iter_max = 1000;
[LMf,LMg] = LevenbergMarquardt(nt,N,tol,iter_max);

%
figure;
% subplot(2,1,1);
hold on;
plot((1:length(LMf))',LMf,'Linewidth',2,'Marker','.','Markersize',20,'Displayname','LM');
legend();
grid;
set(gca,'YScale','log','Fontsize',fsz);
xlabel('iteration num','Fontsize',fsz);
ylabel('loss function','Fontsize',fsz);
% subplot(2,1,2);
% hold on;
% plot((1:length(LMg))',LMg,'Linewidth',2,'Marker','.','Markersize',20,'Displayname','LM');
% legend();
% grid;
% set(gca,'YScale','log','Fontsize',fsz);
% xlabel('k','Fontsize',fsz);
% ylabel('|| grad f||','Fontsize',fsz);
end