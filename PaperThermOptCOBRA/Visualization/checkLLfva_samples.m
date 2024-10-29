clear

load('TOF_samples_TIC_rxns_LG_TOCS_m2.mat')
ids_tic= find(ismember(m.rxns,TIC_Rxns));
ids_other = setdiff([1:numel(m.rxns)],ids_tic);

mi_samples_chrr = min(chrr_samples,[],2);
ma_samples_chrr = max(chrr_samples,[],2);

mi_samples_tof = min(chrr_samples_no_tic,[],2);
ma_samples_tof = max(chrr_samples_no_tic,[],2);
[mi,ma]= fluxVariability(m,'optPercentage',100);
[mi_ll,ma_ll] = fluxVariability(m,'optPercentage',100,'allowLoops','LLC-NS');

% % for the case of with loop analysis
% figure()
% scatter(mi_samples_chrr(ids_other),mi(ids_other),'red')
% hold on
% scatter(mi_samples_chrr(ids_tic),mi(ids_tic),'blue')
% 
% figure()
% scatter(ma_samples_chrr(ids_other),ma(ids_other),'red')
% hold on
% scatter(ma_samples_chrr(ids_tic),ma(ids_tic),'blue')
% 
% figure()
% scatter(mi_samples_tof(ids_other),mi_ll(ids_other),'red')
% hold on
% scatter(mi_samples_tof(ids_tic),mi_ll(ids_tic),'blue')
% 
% figure()
% scatter(ma_samples_tof(ids_other),ma_ll(ids_other),'red')
% hold on
% scatter(ma_samples_tof(ids_tic),ma_ll(ids_tic),'blue')

% plots for only TIC reactions
txt_angle=30;
subplot(2,2,1)
scatter(mi_samples_chrr(ids_tic),mi(ids_tic),'blue')
% bold xlabel
xlabel('Minimum flux obtained across all the CHRR samples')
ylabel({'Minimum flux obtained from','the Flux variability analysis'})
xlim([-1000,0])
ylim([-1000,0])
hold on
plot([-1000,0],[-1000,0],'red')
text(-500,-450,'Slope = 1','Color','red','FontWeight','bold','HorizontalAlignment','center',...
    'Rotation', txt_angle)
set(gca,'FontWeight','bold')

subplot(2,2,2)
scatter(ma_samples_chrr(ids_tic),ma(ids_tic),'blue')
xlabel('Maximum flux obtained across all the CHRR samples')
ylabel({'Maximum flux obtained from','the Flux variability analysis'})
xlim([0,1000])
ylim([0,1000])
hold on
plot([0,1000],[0,1000],'red')
text(500,550,'Slope = 1','Color','red','FontWeight','bold','HorizontalAlignment','center',...
    'Rotation', txt_angle)
set(gca,'FontWeight','bold')

subplot(2,2,3)
scatter(mi_samples_tof(ids_tic),mi_ll(ids_tic),'blue')
xlabel('Minimum flux obtained across all the TOF samples')
ylabel({'Minimum flux obtained from the','loop less Flux variability analysis'})
xlim([-250,0])
ylim([-800,0])
hold on
plot([-250,0],[-800,0],'red')
text(-125,-350,'Slope = 1','Color','red','FontWeight','bold','HorizontalAlignment','center',...
    'Rotation', txt_angle)
set(gca,'FontWeight','bold')

subplot(2,2,4)
scatter(ma_samples_tof(ids_tic),ma_ll(ids_tic),'blue')
xlabel('Maximum flux obtained across all the TOF samples')
ylabel({'Maximum flux obtained from the','loop less Flux variability analysis'})
xlim([0,250])
ylim([0,800])
hold on
plot([0,250],[0,800],'red')
text(125,450,'Slope = 1','Color','red','FontWeight','bold','HorizontalAlignment','center',...
    'Rotation', txt_angle)
set(gca,'FontWeight','bold')