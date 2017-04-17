% 2012 11 08  Plot fish length distribution

addpath /mnt/storage/net_tows
[L0,nfish,ind] = get_fish_length_info(1,2008,361);  % [cm]
L_bin = L0/100;  % convert [cm] to [m]
L_dist = nfish;

figure;
bar(L_bin*100,L_dist,'barwidth',0.8,'facecolor','k');
axis([18 30 0 0.21]);
xlabel('Fish length (cm)','fontsize',14);
ylabel('PDF','fontsize',14);
set(gca,'fontsize',12,'xtick',18:2:30,'ytick',0:0.05:0.2);

saveas(gcf,'/mnt/storage/bb_echopdf_figs/fish_len_dist/fish_len_dist.fig','fig');
saveas(gcf,'/mnt/storage/bb_echopdf_figs/fish_len_dist/fish_len_dist.png','png');