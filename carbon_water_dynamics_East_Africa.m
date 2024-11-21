clc;
clear all;
%% read in the data. 
% rangeland = 'C:\Users\VOdongo\OneDrive - CGIAR\ILRIFlux\Manuscript\Submission\JGR-Biogeosciences\Data\KE_KAP_gapfilledNEE_LE.txt'; % Ensure the path to the data is defined e.g. 'C:\Users\Data\KE_KAP_gapfilledNEE_LE.txt'
% cropland = 'C:\Users\VOdongo\OneDrive - CGIAR\ILRIFlux\Manuscript\Submission\JGR-Biogeosciences\Data\KE_AUS_gapfilledNEE_LE.txt'; % Ensure the path to the data is defined e.g. 'C:\Users\Data\KE_AUS_gapfilledNEE_LE.txt'
% aus_biomet = 'C:\Users\VOdongo\OneDrive - CGIAR\ILRIFlux\Manuscript\Submission\JGR-Biogeosciences\Data\KE_AUS_cropland_biomet.csv';% Ensure the path to the data is defined e.g. 'C:\Users\Data\KE_AUS_cropland_biomet.csv'
% kap_biomet = 'C:\Users\VOdongo\OneDrive - CGIAR\ILRIFlux\Manuscript\Submission\JGR-Biogeosciences\Data\KE_KAP_rangeland_biomet.xlsx';% Ensure the path to the data is defined e.g. 'C:\Users\Data\KE_KAP_rangeland_biomet.xlsx'

rangeland = '.\Data\KE_KAP_gapfilledNEE_LE.txt'; % ./ denotes the current directory where the script is located, and Data is a folder within that directory. Or Ensure the path to the data is defined e.g. 'C:\Users\Data\KE_KAP_gapfilledNEE_LE.txt'
cropland = '.\Data\KE_AUS_gapfilledNEE_LE.txt'; 
aus_biomet = '.\Data\KE_AUS_cropland_biomet.csv';
kap_biomet = '.\Data\KE_KAP_rangeland_biomet.xlsx';
%%
[DateTime,Year,DoY,Hour,NEE,QFVar,QFVar_modified,Ustar,Tair,rH,Rg,VPD,LE,qc_LE,qc_LEmodified,NEE_orig,NEE_f,NEE_fqc,NEE_fall,NEE_fall_qc,NEE_fnum,NEE_fsd,NEE_fmeth,NEE_fwin,Tair_orig,Tair_f,Tair_fqc,Tair_fall,Tair_fall_qc,Tair_fnum,Tair_fsd,Tair_fmeth,Tair_fwin,VPD_orig,VPD_f,VPD_fqc,VPD_fall,VPD_fall_qc,VPD_fnum,VPD_fsd,VPD_fmeth,VPD_fwin,Rg_orig,Rg_f,Rg_fqc,Rg_fall,Rg_fall_qc,Rg_fnum,Rg_fsd,Rg_fmeth,Rg_fwin,LE_orig,LE_f,LE_fqc,LE_fall,LE_fall_qc,LE_fnum,LE_fsd,LE_fmeth,LE_fwin,PotRad_NEW,FP_VARnight,FP_VARday,NEW_FP_Temp,NEW_FP_VPD,FP_RRef_Night,FP_qc,FP_dRecPar,FP_errorcode,FP_GPP2000,FP_k,FP_beta,FP_alpha,FP_RRef,FP_E0,FP_k_sd,FP_beta_sd,FP_alpha_sd,FP_RRef_sd,FP_E0_sd,Reco_DT,GPP_DT,Reco_DT_SD,GPP_DT_SD]= KE_KAP_flux_rangeland(rangeland,3, 8930);
[DateTime1,Year1,DoY1,Hour1,NEE1,QFVar1,QFVar_modified1,Ustar1,Tair1,rH1,Rg1,VPD1,LE1,LE_QC,LE_QC_modified,NEE_orig1,NEE_f1,NEE_fqc1,NEE_fall1,NEE_fall_qc1,NEE_fnum1,NEE_fsd1,NEE_fmeth1,NEE_fwin1,Tair_orig1,Tair_f1,Tair_fqc1,Tair_fall1,Tair_fall_qc1,Tair_fnum1,Tair_fsd1,Tair_fmeth1,Tair_fwin1,VPD_orig1,VPD_f1,VPD_fqc1,VPD_fall1,VPD_fall_qc1,VPD_fnum1,VPD_fsd1,VPD_fmeth1,VPD_fwin1,Rg_orig1,Rg_f1,Rg_fqc1,Rg_fall1,Rg_fall_qc1,Rg_fnum1,Rg_fsd1,Rg_fmeth1,Rg_fwin1,LE_orig1,LE_f1,LE_fqc1,LE_fall1,LE_fall_qc1,LE_fnum1,LE_fsd1,LE_fmeth1,LE_fwin1,PotRad_NEW1,FP_VARnight1,FP_VARday1,NEW_FP_Temp1,NEW_FP_VPD1,FP_RRef_Night1,FP_qc1,FP_dRecPar1,FP_errorcode1,FP_GPP1,FP_k1,FP_beta1,FP_alpha1,FP_RRef1,FP_E1,FP_k_sd1,FP_beta_sd1,FP_alpha_sd1,FP_RRef_sd1,FP_E0_sd1,Reco_DT1,GPP_DT1,Reco_DT_SD1,GPP_DT_SD1] = KE_AUS_flux_cropland(cropland,3, 8930);

%% Read biomet data for Kapiti (rangeland)
 [DateTime2,DOY_meteodata,Ta_1_1_1,TA_1_2_1,Rn_1_1_1,SHF_1_1_1,SHF_2_1_1,SHF_3_1_1,RH_1_1_1,SWin_1_1_1,SWC_1_1_1,SWC_2_1_1,SWC_3_1_1,Ts_1_1_1,Ts_2_1_1,Ts_3_1_1,PPFD_1_1_1,P_rain_1_1_1,DOY_meteofull,~,H,LE_range,co2_flux,qc_co2_flux] = KE_KAP_rangeland_biomet(kap_biomet,'Sheet1',2,17521);% 
%% Read biomet data for Ausquest (cropland)
[TIMESTAMP,RECORD,FC,FC_mass,FC_QC,FC_samples,LE1_crop,LE_QC_crop,LE_samples,H1,H_QC,H_samples,NETRAD,G,SG,energy_closure,poor_energy_closure_flg,Bowen_ratio,TAU,TAU_QC,USTAR,TSTAR,TKE,Tair2,TA_SIGMA_1_1_1,RH,T_DP_1_1_1,e,e_sat,TA_2_1_1,RH_2_1_1,T_DP_2_1_1,e_probe,e_sat_probe,H2O_probe,PA,PA_SIGMA,VPD2,U,U_SIGMA,V,V_SIGMA,W,W_SIGMA,T_SONIC,T_SONIC_SIGMA,sonic_azimuth,WS,WS_RSLT,WD_SONIC,WD_SIGMA,WD,WS_MAX,CO2,CO2_SIGMA,CO2_density,CO2_density_SIGMA,H2O,H2O_SIGMA,H2O_density,H2O_density_SIGMA,CO2_sig_strgth_Min,H2O_sig_strgth_Min,Rain,NETRAD_meas,SW_IN,PPFD1,PPFD2,sun_azimuth,sun_elevation,hour_angle,sun_declination,air_mass_coeff,daytime,Tsoil,TS_2_1_1,SWC,SWC_2_1_2,cs65x_ec_1_1_1,cs65x_ec_2_1_1,G_PLATE_1_1_1,G_PLATE_2_1_1,shfp_cal_1_1_1,shfp_cal_2_1_1,FETCH_MAX,FETCH_90,FETCH_55,FETCH_40,UPWND_DIST_INTRST,FTPRNT_DIST_INTRST,FTPRNT_EQUATION] = KE_AUS_cropland_biomet(aus_biomet,2, 10849);

%% define the coincidence period of cropping with rangeland partition
start_datetime = datenum('2019-03-13 00:30:00');
end_datetime = datenum('2019-09-14 18:00:00');

start_date = '03-2019';
end_date = '10-2019';

start_date_num = datenum(start_date, 'mm-yyyy');
end_date_num = datenum(end_date, 'mm-yyyy');

%% Create logical index for selecting data within the specified time range
idx = DateTime >= start_datetime & DateTime <= end_datetime;
idx1 = DateTime1 >= start_datetime & DateTime1 <= end_datetime;
idx2 = DateTime2 >= start_datetime & DateTime2 <= end_datetime;
idx3 = TIMESTAMP >= start_datetime & TIMESTAMP <= end_datetime;

DateTime = DateTime(idx);
DateTime1 = DateTime(idx1);
%% 
P_rain_1_1_1(P_rain_1_1_1<0)            = NaN;
SWC_3_1_1(SWC_3_1_1<0)                  = NaN;
Tair(Tair < 0)                          = NaN;
PPFD_1_1_1(PPFD_1_1_1 < 0)              = NaN;
NEE(NEE<-50)                           = NaN;

VPD              = VPD(idx)/100;
VPD(VPD<0)       = NaN;
VPD1(VPD1<0)     = NaN;
vpd_a_d = VPD1(idx1)*0.1; % covert from hPa to kPa

P_rain_r = P_rain_1_1_1(idx2);
SWC_3_1_1_r = SWC_3_1_1(idx2);
PPFD_1_1_1_r =PPFD_1_1_1(idx2);
LE_r = LE_f(idx);
LE1_a = LE_f1(idx1);
ET_r = 1000*LE_r*1800/(2.45*10^6*1000);
ET_a = 1000*LE1_a*1800/(2.45*10^6*1000);
ET_gr = LE(idx);
ET_ga = LE1(idx);
ET_gap_r = 1000*ET_gr*1800/(2.45*10^6*1000);
ET_gap_a = 1000*ET_ga*1800/(2.45*10^6*1000);
%%
nee_r  = NEE_f(idx)*1E-6*12*1800;
gpp_r  = GPP_DT(idx)*1E-6*12*1800;
reco_r = Reco_DT(idx)*1E-6*12*1800;
%%
gpp_r_d =gpp_r;
nee_r_d =nee_r;
reco_r_d =reco_r;
VPD_r_d=VPD*0.1;

%%
nee_a  = NEE_f1(idx1)*1E-6*12*1800;
nee_rr  = NEE(idx)*1E-6*12*1800;
nee_aa  = NEE1(idx1)*1E-6*12*1800;
gpp_a  = GPP_DT1(idx1)*1E-6*12*1800;
reco_a = Reco_DT1(idx1)*1E-6*12*1800;
%%
Rain(Rain<0)                              = NaN;
Tair1(Tair1 < 0)                          = NaN;
PPFD1(PPFD1 < 0)                          = NaN;
NEE1(NEE1<-50)                          = NaN;
nee_a(nee_a<-50)                          = NaN;
nee_aa(nee_aa<-50)                          = NaN;
nee_rr(nee_rr<-50)                          = NaN;

%% Figure 1: cumulative fluxes before adjusting for lateral carbon
figure('Position', [100, 100, 1024, 768]); % Adjust the values as needed
subplot(3,1,1)
plot(DateTime, cumsum(nee_r), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'LineStyle', '-'); % Orange solid line
xlim([start_date_num, end_date_num]);

hold on
ax = gca;
x_ticks = get(ax, 'XTick');
empty_labels = cell(size(x_ticks));
set(ax, 'XTickLabel', empty_labels);
plot(DateTime1, cumsum(nee_a), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'LineStyle', '--'); % Blue dashed line
%% Draw a horizontal line at y=0
line([start_date_num, end_date_num], [0, 0], 'Color', 'k', 'LineWidth', 1);
hold on

xlim([start_date_num, end_date_num]);
leg = findobj(gcf, 'Type', 'Legend');
legend('Rangeland','Cropland','Location','NorthWest');
ax=gca;
set(ax,'FontSize',16)
title('(a) NEE');

subplot(3,1,2)
% plot(DateTime,cumsum(gpp_r),'r','LineWidth',2);
plot(DateTime, cumsum(gpp_r), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'LineStyle', '-'); % Orange solid line
xlim([start_date_num, end_date_num]);
datetick('x','mmm-yyyy','keeplimits');
ylabel('gC m^{-2}', 'FontSize', 14); % Increase font size for y-labelset(ax,'FontSize',14);
hold on

% plot(DateTime1,cumsum(gpp_a),'b','LineWidth',2)
plot(DateTime1, cumsum(gpp_a), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'LineStyle', '--'); % Blue dashed line
xlim([start_date_num, end_date_num]);
datetick('x','mmm-yyyy','keeplimits');
% set(gca, 'XTick', []);
ax = gca;
x_ticks = get(ax, 'XTick');
empty_labels = cell(size(x_ticks));
set(ax, 'XTickLabel', empty_labels);
set(get(ax, 'YLabel'), 'FontSize', 16);
set(ax,'FontSize',16)
title('(b) GPP');

subplot(3,1,3)
% plot(DateTime,cumsum(reco_r),'r','LineWidth',2);
plot(DateTime, cumsum(reco_r), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'LineStyle', '-'); % Orange solid line
xlim([start_date_num, end_date_num]);

% Set the x-axis to date ticks
datetick('x','mmm-yyyy','keeplimits');

hold on
% plot(DateTime1,cumsum(reco_a),'b','LineWidth',2);
plot(DateTime1, cumsum(reco_a), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'LineStyle', '--'); % Blue dashed line
xlim([start_date_num, end_date_num]);
datetick('x','mmm-yyyy','keeplimits');

ax=gca;
set(ax,'FontSize',14)
set(leg, 'FontSize', 16);
title('(c) R{_e_c_o}');

%%Lateral carbon adjustment after harvest

% Calculate cumulative NEE for cropland before the adjustment
cumulative_nee_a = cumsum(nee_a);

% Adjusted NEE for cropland by adding the lateral C flux at the last point
FW_grain = 1200; % Fresh grain weight harvest in Kg per ha
MC=11; % Percent grain moisture content at harvest
lC_fraction = 49; % Percent legume carbon fraction
H_index = 56; % Percent harvest index
%% Computing for lateral carbon
DW_grain = FW_grain*(1-MC/100);
DW_stover = DW_grain*(1-H_index/100);
ABGM = DW_grain + DW_stover; % Above ground total dry biomass
C_harvest = ABGM * lC_fraction/100; % Harvested lateral carbon in KgC/ha
C_harvest = C_harvest*1000/10000; % Harvested lateral carbon in gC/m2

cumulative_nee_a_adjusted = cumulative_nee_a + C_harvest;


% Plot cumulative fluxes
figure('Position', [100, 100, 1024, 768]); % Adjust the values as needed

% Subplot 1: NEE
subplot(3,1,1)
plot(DateTime, cumsum(nee_r), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'LineStyle', '-'); % Orange solid line
hold on
plot(DateTime1, cumulative_nee_a, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'LineStyle', '--'); % Blue dashed line

% Plot the adjustment as a separate line in a different color
plot([DateTime1(end), DateTime1(end)], [cumulative_nee_a(end), cumulative_nee_a_adjusted(end)], ...
    'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 2.5, 'LineStyle', '-.'); % Green dash-dot line

% Draw a horizontal line at y=0
line([start_date_num, end_date_num], [0, 0], 'Color', 'k', 'LineWidth', 1);

xlim([start_date_num, end_date_num]);
ax = gca;
set(ax, 'FontSize', 16);
x_ticks = get(ax, 'XTick');
empty_labels = cell(size(x_ticks));
set(ax, 'XTickLabel', empty_labels);

% Update the legend
legend('Rangeland', 'Cropland (Before Adjustment)', 'Lateral C Flux Adjustment', 'Location', 'NorthWest');

title('(a) NEE');

% Subplot 2: GPP
subplot(3,1,2)
plot(DateTime, cumsum(gpp_r), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'LineStyle', '-'); % Orange solid line
hold on
plot(DateTime1, cumsum(gpp_a), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'LineStyle', '--'); % Blue dashed line
xlim([start_date_num, end_date_num]);
ylabel('gC m^{-2}', 'FontSize', 14); % Increase font size for y-label
datetick('x', 'mmm-yyyy', 'keeplimits');
ax = gca;
x_ticks = get(ax, 'XTick');
empty_labels = cell(size(x_ticks));
set(ax, 'XTickLabel', empty_labels);
set(ax, 'FontSize', 16);
title('(b) GPP');

% Subplot 3: Reco
subplot(3,1,3)
plot(DateTime, cumsum(reco_r), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'LineStyle', '-'); % Orange solid line
hold on
plot(DateTime1, cumsum(reco_a), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'LineStyle', '--'); % Blue dashed line
xlim([start_date_num, end_date_num]);
datetick('x', 'mmm-yyyy', 'keeplimits');
ax = gca;
set(ax, 'FontSize', 16);
title('(c) R_{eco}');

%%
a=[year(DateTime),month(DateTime), day(DateTime)];
a2=[year(DateTime1),month(DateTime1), day(DateTime1)];
b=unique(a,'rows','stable');
b2=unique(a2,'rows','stable');
d_gpp_r=arrayfun(@(x)sum(gpp_r_d(ismember(a,b(x,:),'rows'))),1:length(b))';
d_nee_r=arrayfun(@(x)sum(nee_r_d(ismember(a,b(x,:),'rows'))),1:length(b))';
d_nee_rr=arrayfun(@(x)sum(nee_rr(ismember(a,b(x,:),'rows'))),1:length(b))';
d_reco_r=arrayfun(@(x)sum(reco_r_d(ismember(a,b(x,:),'rows'))),1:length(b))';
d_ppt=arrayfun(@(x)nansum(P_rain_r(ismember(a,b(x,:),'rows'))),1:length(b))';
d_swc=arrayfun(@(x)nanmean(SWC_3_1_1_r(ismember(a,b(x,:),'rows'))),1:length(b))';
d_evap_r=arrayfun(@(x)nansum(ET_r(ismember(a,b(x,:),'rows'))),1:length(b))';
d_Ts_r=arrayfun(@(x)nanmean(Ts_3_1_1(ismember(a,b(x,:),'rows'))),1:length(b))';
d_Ta_r=arrayfun(@(x)nanmean(Tair(ismember(a,b(x,:),'rows'))),1:length(b))';
d_par_r=arrayfun(@(x)nanmean(PPFD_1_1_1_r(ismember(a,b(x,:),'rows'))),1:length(b))';
d_vpd_r=arrayfun(@(x)nanmean(VPD_r_d(ismember(a,b(x,:),'rows'))),1:length(b))';

newMatrix=[b,d_gpp_r, d_nee_r, d_reco_r];
d_date_r =datenum(b);

%% Daily carbon partition for Agriculture
a_a=[year(DateTime1),month(DateTime1), day(DateTime1)];
a_a1=[year(TIMESTAMP(idx3)),month(TIMESTAMP(idx3)), day(TIMESTAMP(idx3))];
b_a=unique(a_a,'rows','stable');
b_a1=unique(a_a1,'rows','stable');
d_gpp_a=arrayfun(@(x)sum(gpp_a(ismember(a_a,b_a(x,:),'rows'))),1:length(b_a))';
d_nee_a=arrayfun(@(x)sum(nee_a(ismember(a_a,b_a(x,:),'rows'))),1:length(b_a))';
d_nee_aa=arrayfun(@(x)sum(nee_aa(ismember(a_a,b_a(x,:),'rows'))),1:length(b_a))';
d_reco_a=arrayfun(@(x)sum(reco_a(ismember(a_a,b_a(x,:),'rows'))),1:length(b_a))';
d_ppt_a=arrayfun(@(x)nansum(Rain(ismember(a_a,b_a(x,:),'rows'))),1:length(b_a))';
d_evap_a=arrayfun(@(x)nansum(ET_a(ismember(a_a,b_a1(x,:),'rows'))),1:length(b_a1))';
d_Ts_a=arrayfun(@(x)nanmean(Ts_1_1_1(ismember(a,b(x,:),'rows'))),1:length(b))';
d_Ta_a=arrayfun(@(x)nanmean(Tair1(ismember(a,b(x,:),'rows'))),1:length(b))';
d_par_a=arrayfun(@(x)nanmean(PPFD1(ismember(a_a,b_a(x,:),'rows'))),1:length(b_a))';
d_vpd_a=arrayfun(@(x)nanmean(vpd_a_d(ismember(a,b(x,:),'rows'))),1:length(b))';

newMatrix_a=[b,d_gpp_a, d_nee_a, d_reco_a];
d_date_a =datenum(b_a);

%% Figure 2: C dynamics (daily GPP, Reco, NEE) and essential meteorological variables (rainfall and soil moisture) for the (a) rangeland and (b) cropland site
figure('Position', [100, 100, 1024, 768]); % Adjust the values as needed

% Define the tick positions for each month
x_tick_positions = datenum({'01-Mar-2019', '01-Apr-2019', '01-May-2019', '01-Jun-2019', '01-Jul-2019', '01-Aug-2019', '01-Sep-2019', '01-Oct-2019'});

% First subplot (Rangeland - NEE, GPP, Reco)
subplot(2,2,1);
hold on;
p1 = plot(d_date_r, d_nee_rr, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); % Orange
p2 = plot(d_date_r, d_gpp_r, 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250], 'LineStyle', '--'); % Yellowish
p3 = plot(d_date_r, d_reco_r, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); % Sky blue
line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineWidth', 1); % Black line at y=0
hold off;

ylabel('gC m^{-2} day^{-1}');
legend([p1, p2, p3], 'NEE', 'GPP', 'R_{eco}', 'Location', 'northwest', 'Interpreter', 'tex');
ax1 = gca;
set(ax1, 'FontSize', 16);
set(ax1, 'Box', 'on'); % Add the black box around the axes
set(get(ax1, 'YLabel'), 'FontSize', 16);
ylim([-10 15]);
xlim([min(d_date_r) max(d_date_r)]);
set(ax1, 'XTick', x_tick_positions); % Set tick positions for each month
datetick('x', 'mmm-yyyy', 'keeplimits', 'keepticks'); % Display month and year on the x-axis
set(ax1, 'XTickLabel', {}); % Hide x-axis labels for the upper subplot
tickValues1 = -10:5:15;
tickLabels1 = cellstr(num2str(tickValues1', '%.0f'));
set(ax1, 'YTick', tickValues1, 'YTickLabel', tickLabels1);

% Position and display the title inside the plot
title('(a) Rangeland', 'FontSize', 16);

% Second subplot (Rangeland - Rain and Soil Moisture)
subplot(2,2,3);
[hAx,hLine1,hLine2] = plotyy(d_date_r, d_ppt*1000, d_date_r, d_swc, 'bar', 'plot');
set(hLine1, 'FaceColor', [0, 0.4470, 0.7410]); % Set bars to blue color
soilMoistureColor = [0.4940, 0.1840, 0.5560]; % Purplish color used for soil moisture line
set(hLine2, 'LineStyle', '--', 'LineWidth', 2, 'Color', soilMoistureColor);
xlim(hAx(1), [min(d_date_r) max(d_date_r)]);
xlim(hAx(2), [min(d_date_r) max(d_date_r)]);
set(hAx(1), 'XTick', x_tick_positions); % Set tick positions for each month
set(hAx(2), 'XTick', x_tick_positions); % Set tick positions for each month
datetick(hAx(1), 'x', 'mmm-yyyy', 'keeplimits', 'keepticks'); % Display month and year on the x-axis
datetick(hAx(2), 'x', 'mmm-yyyy', 'keeplimits', 'keepticks');
ylabel(hAx(1), 'Rain (mm)');
ylabel(hAx(2), 'swc (m^3/m^3)', 'Color', soilMoistureColor, 'FontSize', 16);
set(hAx, 'FontSize', 16);
set(get(hAx(1), 'YLabel'), 'FontSize', 16);
set(hAx(2), 'ycolor', soilMoistureColor, 'FontSize', 16);
lgd = legend([hLine1(1), hLine2(1)], 'Rain', 'Soil moisture');
set(lgd, 'Location', 'northwest');
set(lgd, 'Box', 'off');
ylim(hAx(1), [0, 80]);
ylim(hAx(2), [0, 1]);
set(hAx(1), 'YTick', 0:10:80);
set(hAx(2), 'YTick', 0:0.25:1);
set(hAx(1), 'Box', 'on');
set(hAx(2), 'Box', 'on');

% Third subplot (Cropland - NEE, GPP, Reco)
subplot(2,2,2);
hold on;
p1 = plot(d_date_a, d_nee_aa, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); % Orange
p2 = plot(d_date_a, d_gpp_a, 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250], 'LineStyle', '--'); % Yellowish
p3 = plot(d_date_a, d_reco_a, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); % Sky blue
line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineWidth', 1); % Black line at y=0
hold off;

ylabel('gC m^{-2} day^{-1}');
legend([p1, p2, p3], 'NEE', 'GPP', 'R_{eco}', 'Location', 'northwest', 'Interpreter', 'tex');
ax2 = gca;
set(ax2, 'FontSize', 16);
set(ax2, 'Box', 'on'); % Add the black box around the axes
set(get(ax2, 'YLabel'), 'FontSize', 16);
ylim([-10 15]);
xlim([min(d_date_a) max(d_date_a)]);
set(ax2, 'XTick', x_tick_positions); % Set tick positions for each month
datetick('x', 'mmm-yyyy', 'keeplimits', 'keepticks'); % Display month and year on the x-axis
set(ax2, 'XTickLabel', {}); % Hide x-axis labels for the upper subplot
tickValues2 = -10:5:15;
tickLabels2 = cellstr(num2str(tickValues2', '%.0f'));
set(ax2, 'YTick', tickValues2, 'YTickLabel', tickLabels2);

% Position and display the title inside the plot
title('(b) Cropland', 'FontSize', 16);

% Fourth subplot (Cropland - Rain)
subplot(2,2,4);
hBar = bar(d_date_a, d_ppt_a, 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'b'); % Set bars to sky blue color with a blue edge
hAx = gca; % Current axes handle
xlim(hAx, [min(d_date_a) max(d_date_a)]);
set(hAx, 'XTick', x_tick_positions); % Set tick positions for each month
datetick(hAx, 'x', 'mmm-yyyy', 'keeplimits', 'keepticks'); % Display month and year on the x-axis
ylabel(hAx, 'Rain (mm)');
set(hAx, 'FontSize', 16);
set(get(hAx, 'YLabel'), 'FontSize', 16);
ylim(hAx, [0, 80]);
set(hAx, 'YTick', 0:10:80);
set(hAx, 'Box', 'on');

% Adjust the legend to show only 'Rain'
lgd = legend(hBar, 'Rain');
set(lgd, 'Location', 'northwest');
set(lgd, 'Box', 'off');


%% Figure 3: Correlation between GPP and Reco for the rangeland and cropland. The stronger correlation in the cropland suggests strong coupling between GPP and Reco.1:1 Reco vs GPP plot
% Define colors for readability
colorRangeland = [0.5, 0, 0];    % Dark red/brown
colorCropland = [0, 0, 0.5];     % Dark blue
figure;
h1 = plot(d_gpp_r, d_reco_r, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorRangeland, 'MarkerEdgeColor', colorRangeland);
hold on;
h2 = plot(d_gpp_a, d_reco_a, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorCropland, 'MarkerEdgeColor', colorCropland);
xlim([0 12]);
ylim([0 10]);
xlabel('GPP (gC m^{-2} d^{-1})', 'FontSize', 14);
ylabel('R_{eco} (gC m^{-2} d^{-1})', 'FontSize', 14);
title('Whole Growing Season', 'FontSize', 14);
set(gca, 'FontSize', 16);

% Fit linear models and calculate goodness of fit
[fit_r, gof_r] = fit(d_gpp_r, d_reco_r, 'poly1');
[fit_a, gof_a] = fit(d_gpp_a, d_reco_a, 'poly1');

% Calculate correlation coefficients (r)
r_r = sqrt(gof_r.rsquare);
r_a = sqrt(gof_a.rsquare);

% Define x range for plotting fit lines
x_range_r = linspace(min(d_gpp_r), max(d_gpp_r), 100);
x_range_a = linspace(min(d_gpp_a), max(d_gpp_a), 100);

% Plot fits and match colors with thicker lines
h3 = plot(x_range_r, feval(fit_r, x_range_r), '--', 'Color', colorRangeland, 'LineWidth', 2);
h4 = plot(x_range_a, feval(fit_a, x_range_a), '--', 'Color', colorCropland, 'LineWidth', 2);

% Add legend for the main plot only
legend([h1, h2, h3, h4], {'Rangeland Data', 'Cropland Data', sprintf('Fit Rangeland: \\it{r}\\rm=%.2f', r_r), sprintf('Fit Cropland: \\it{r}\\rm=%.2f', r_a)}, 'Location', 'northeast');
%% Carbon Use efficiency (CUE) comparisons: Observation period vs peak growing period
%Figure 6: Comparison of NEP and CUE between rangeland and cropland sites. (a) and (b) 
% show the relationship between GPP and NEP for the whole observation period and the peak 
%growing season (90th percentile of NEE) (c) CUE calculated for the entire observation period, 
%showing significantly higher CUE in cropland compared to rangeland (*** indicates p<0.001). 
%(d) CUE during the peak growing period, similarly showing significantly higher CUE in cropland 
%compared to rangeland (p<0.001). Total daily observations (N) are indicated on each bar. 
%Error bar for the cropland observation period is not displayed as it was too large (s.e.m = 25.6). 
figure('Position', [100, 100, 1024, 768]);

% Define colors for readability
colorRangeland = [0.5, 0, 0];    % Dark red/brown
colorCropland = [0, 0, 0.5];     % Dark blue

% Calculate the 90th percentile of NEE for peak growing season identification
nee_90th_percentile_r = prctile(-d_nee_r, 90);
nee_90th_percentile_a = prctile(-d_nee_a, 90);

% Identify indices for the peak growing season based on the 90th percentile
idx_peak_r = -d_nee_r >= nee_90th_percentile_r;  % NEP more negative than the 90th percentile
idx_peak_a = -d_nee_a >= nee_90th_percentile_a;  % NEP more negative than the 90th percentile

% First subplot for whole season data
subplot(2,2,1); % Position (a)
h1 = plot(d_gpp_r, -d_nee_r, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorRangeland, 'MarkerEdgeColor', colorRangeland);
hold on;
h2 = plot(d_gpp_a, -d_nee_a, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorCropland, 'MarkerEdgeColor', colorCropland);
xlim([0 15]);
ylim([-9 10]);
xlabel('GPP (gC m^{-2} d^{-1})', 'FontSize', 14);
ylabel('NEP (gC m^{-2} d^{-1})', 'FontSize', 14);
title('Observation period', 'FontSize', 14);
text(0.5, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 14,'HorizontalAlignment', 'right');
set(gca, 'FontSize', 12);

% Fit linear models and calculate goodness of fit
[fit_r, gof_r] = fit(d_gpp_r, -d_nee_r, 'poly1');
[fit_a, gof_a] = fit(d_gpp_a, -d_nee_a, 'poly1');

% Calculate slopes
slope_r = fit_r.p1;
slope_a = fit_a.p1;

% Calculate correlation coefficients (r)
r_r = sqrt(gof_r.rsquare);
r_a = sqrt(gof_a.rsquare);

% Define x range for plotting fit lines
x_range_r = linspace(min(d_gpp_r), max(d_gpp_r), 100);
x_range_a = linspace(min(d_gpp_a), max(d_gpp_a), 100);

% Plot fits and match colors with thicker lines
h3 = plot(x_range_r, feval(fit_r, x_range_r), '--', 'Color', colorRangeland, 'LineWidth', 2);
h4 = plot(x_range_a, feval(fit_a, x_range_a), '--', 'Color', colorCropland, 'LineWidth', 2);

% Add legend for the main plot (subplot (a)) with slopes
legend([h1, h2, h3, h4], {'Rangeland Data', 'Cropland Data', ...
    sprintf('Rangeland: \\it{r}\\rm=%.2f, Slope=%.2f', r_r, slope_r), ...
    sprintf('Cropland: \\it{r}\\rm=%.2f, Slope=%.2f', r_a, slope_a)}, ...
    'Location', 'northeast');

% Second subplot for peak growing season (based on 90th percentile)
subplot(2,2,2); % Position (b)
d_gpp_rg = d_gpp_r(idx_peak_r);
d_nee_rg = -d_nee_r(idx_peak_r);
d_gpp_ag = d_gpp_a(idx_peak_a);
d_nee_ag = -d_nee_a(idx_peak_a);
h5 = plot(d_gpp_rg, d_nee_rg, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorRangeland, 'MarkerEdgeColor', colorRangeland);
hold on;
h6 = plot(d_gpp_ag, d_nee_ag, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorCropland, 'MarkerEdgeColor', colorCropland);
xlim([0 15]);
ylim([-9 10]);
xlabel('GPP (gC m^{-2} d^{-1})', 'FontSize', 14);
ylabel('NEP (gC m^{-2} d^{-1})', 'FontSize', 14);
title('Peak growing period (90th Percentile)', 'FontSize', 14);
text(0.5, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 14,'HorizontalAlignment', 'right');
set(gca, 'FontSize', 12);

% Fit linear models and calculate goodness of fit for the peak growing season
[fit_rg, gof_rg] = fit(d_gpp_rg, d_nee_rg, 'poly1');
[fit_ag, gof_ag] = fit(d_gpp_ag, d_nee_ag, 'poly1');

% Calculate slopes
slope_rg = fit_rg.p1;
slope_ag = fit_ag.p1;

% Calculate correlation coefficients (r)
r_rg = sqrt(gof_rg.rsquare);
r_ag = sqrt(gof_ag.rsquare);

% Define x range for plotting fit lines
x_range_rg = linspace(min(d_gpp_rg), max(d_gpp_rg), 100);
x_range_ag = linspace(min(d_gpp_ag), max(d_gpp_ag), 100);

% Plot fits and match colors with thicker lines
h7 = plot(x_range_rg, feval(fit_rg, x_range_rg), '--', 'Color', colorRangeland, 'LineWidth', 2);
h8 = plot(x_range_ag, feval(fit_ag, x_range_ag), '--', 'Color', colorCropland, 'LineWidth', 2);

% Add legend for subplot (b) with slopes
legend([h5, h6, h7, h8], {'Rangeland Data', 'Cropland Data', ...
    sprintf('Rangeland: \\it{r}\\rm=%.2f, Slope=%.2f', r_rg, slope_rg), ...
    sprintf('Cropland: \\it{r}\\rm=%.2f, Slope=%.2f', r_ag, slope_ag)}, ...
    'Location', 'northeast');

% Third subplot for CUE bar plot for the whole season
subplot(2,2,3); % Position (c)
CUE_mean_rangeland = sum(-d_nee_r) / sum(d_gpp_r);
CUE_mean_cropland = sum(-d_nee_a) / sum(d_gpp_a);
SD_rangeland = std(-d_nee_r ./ d_gpp_r);
SD_cropland = std(-d_nee_a ./ d_gpp_a);
n_rangeland = numel(d_nee_r);
n_cropland = numel(d_nee_a);
SEM_rangeland = SD_rangeland / sqrt(n_rangeland);
SEM_cropland = SD_cropland / sqrt(n_cropland);
CUE_means = [CUE_mean_rangeland, CUE_mean_cropland];
CUE_errors = [SEM_rangeland, SEM_cropland];

bar(1, CUE_means(1), 'BarWidth', 0.5, 'FaceColor', colorRangeland);
hold on;
bar(2, CUE_means(2), 'BarWidth', 0.5, 'FaceColor', colorCropland);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Rangeland', 'Cropland'}, 'FontSize', 12);
ylabel('CUE', 'FontSize', 14);
ylim([-0.5 1.0]);
text(0.5, 0.95, '(c)', 'Units', 'normalized', 'FontSize', 14,'HorizontalAlignment', 'right');
errorbar(1, CUE_means(1), CUE_errors(1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
% errorbar(2, CUE_means(2), CUE_errors(2), 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Display number of observations
text(1, CUE_means(1) + 0.05, sprintf('N=%d', n_rangeland), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(2, CUE_means(2) + 0.05, sprintf('N=%d', n_cropland), 'HorizontalAlignment', 'center', 'FontSize', 12);

% Perform t-test between rangeland and cropland
[h, p] = ttest2(-d_nee_r ./ d_gpp_r, -d_nee_a ./ d_gpp_a);

% Add significance stars
if h == 1
    if p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end
    text(1.5, max(CUE_means) + 0.1, sig, 'HorizontalAlignment', 'center', 'FontSize', 16);
end
box on;

% Fourth subplot for CUE bar plot for the peak growing season
subplot(2,2,4); % Position (d)
CUE_mean_rangeland = sum(d_nee_rg) / sum(d_gpp_rg);
CUE_mean_cropland = sum(d_nee_ag) / sum(d_gpp_ag);
SD_rangeland = std(d_nee_rg ./ d_gpp_rg);
SD_cropland = std(d_nee_ag ./ d_gpp_ag);
n_rangeland = numel(d_nee_rg);
n_cropland = numel(d_nee_ag);
SEM_rangeland = SD_rangeland / sqrt(n_rangeland);
SEM_cropland = SD_cropland / sqrt(n_cropland);
CUE_means = [CUE_mean_rangeland, CUE_mean_cropland];
CUE_errors = [SEM_rangeland, SEM_cropland];

bar(1, CUE_means(1), 'BarWidth', 0.5, 'FaceColor', colorRangeland);
hold on;
bar(2, CUE_means(2), 'BarWidth', 0.5, 'FaceColor', colorCropland);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Rangeland', 'Cropland'}, 'FontSize', 12);
ylabel('CUE', 'FontSize', 14);
ylim([-0.5 1.0]);
text(0.5, 0.95, '(d)', 'Units', 'normalized', 'FontSize', 14,'HorizontalAlignment', 'right');
errorbar(1, CUE_means(1), CUE_errors(1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
errorbar(2, CUE_means(2), CUE_errors(2), 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Display number of observations
text(1, CUE_means(1) + 0.05, sprintf('N=%d', n_rangeland), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(2, CUE_means(2) + 0.05, sprintf('N=%d', n_cropland), 'HorizontalAlignment', 'center', 'FontSize', 12);

% Perform t-test between rangeland and cropland
[h, p] = ttest2(d_nee_rg ./ d_gpp_rg, d_nee_ag ./ d_gpp_ag);

% Add significance stars
if h == 1
    if p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end
    text(1.5, max(CUE_means) + 0.1, sig, 'HorizontalAlignment', 'center', 'FontSize', 16);
end
box on;

% Force MATLAB to update the figure immediately
drawnow;
%%
%% Water Use Efficiency(WUE)comparisons: Observation period vs peak growing period
%Figure 7: Comparison of GPP and WUE between rangeland and cropland systems during the 
%observation period and peak growing season. (a) Relationship between ET and GPP for rangeland
% and cropland during the entire observation period. (b) Relationship between ET and GPP during 
%the peak growing period (90th percentile NEE threshold). (c) WUE calculated for the entire observation 
%period, showing no significant difference (ns) between rangeland and cropland. (d) WUE during the peak 
%growing period, showing significantly higher WUE in cropland compared to rangeland (*** indicates p<0.001). 
%The number of observations (N) is indicated on each bar.
figure('Position', [100, 100, 1024, 768]);

% Define colors for readability
colorRangeland = [0.5, 0, 0];    % Dark red/brown
colorCropland = [0, 0, 0.5];     % Dark blue

% Calculate the 90th percentile of NEE for peak growing season identification
nee_90th_percentile_r = prctile(-d_nee_r, 90);
nee_90th_percentile_a = prctile(-d_nee_a, 90);

% Identify indices for the peak growing season based on the 90th percentile
idx_peak_r = -d_nee_r >= nee_90th_percentile_r;  % NEP more negative than the 90th percentile
idx_peak_a = -d_nee_a >= nee_90th_percentile_a;  % NEP more negative than the 90th percentile

% First subplot for whole season data
subplot(2,2,1); % Position (a)
h1 = plot(d_evap_r, d_gpp_r, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorRangeland, 'MarkerEdgeColor', colorRangeland);
hold on;
h2 = plot(d_evap_a, d_gpp_a, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorCropland, 'MarkerEdgeColor', colorCropland);
xlim([0 5]);
ylim([0 17]);
xlabel('ET (mm d^{-1})', 'FontSize', 14);
ylabel('GPP (gC m^{-2} d^{-1})', 'FontSize', 14);
title('Observation period', 'FontSize', 14);
text(0.5, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 14,'HorizontalAlignment', 'right');
set(gca, 'FontSize', 12);

% Fit linear models and calculate goodness of fit
[fit_r, gof_r] = fit(d_evap_r, d_gpp_r, 'poly1');
[fit_a, gof_a] = fit(d_evap_a, d_gpp_a, 'poly1');

% Calculate slopes
slope_r = fit_r.p1;
slope_a = fit_a.p1;

% Calculate correlation coefficients (r)
r_r = sqrt(gof_r.rsquare);
r_a = sqrt(gof_a.rsquare);

% Define x range for plotting fit lines
x_range_r = linspace(min(d_evap_r), max(d_evap_r), 100);
x_range_a = linspace(min(d_evap_a), max(d_evap_a), 100);

% Plot fits and match colors with thicker lines
h3 = plot(x_range_r, feval(fit_r, x_range_r), '--', 'Color', colorRangeland, 'LineWidth', 2);
h4 = plot(x_range_a, feval(fit_a, x_range_a), '--', 'Color', colorCropland, 'LineWidth', 2);

% Add legend for the main plot (subplot (a)) with slopes
legend([h1, h2, h3, h4], {'Rangeland Data', 'Cropland Data', ...
    sprintf('Rangeland: \\it{r}\\rm=%.2f, Slope=%.2f', r_r, slope_r), ...
    sprintf('Cropland: \\it{r}\\rm=%.2f, Slope=%.2f', r_a, slope_a)}, ...
    'Location', 'northeast');

% Second subplot for peak growing season (based on 90th percentile)
subplot(2,2,2); % Position (b)
d_gpp_rg = d_gpp_r(idx_peak_r);
d_evap_rg = d_evap_r(idx_peak_r);
d_gpp_ag = d_gpp_a(idx_peak_a);
d_evap_ag = d_evap_a(idx_peak_a);
h5 = plot(d_evap_rg, d_gpp_rg , 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorRangeland, 'MarkerEdgeColor', colorRangeland);
hold on;
h6 = plot(d_evap_ag, d_gpp_ag, 'o', 'MarkerSize', 5, 'MarkerFaceColor', colorCropland, 'MarkerEdgeColor', colorCropland);
xlim([0 5]);
ylim([0 17]);
xlabel('ET (mm d^{-1})', 'FontSize', 14);
ylabel('GPP (gC m^{-2} d^{-1})', 'FontSize', 14);
title('Peak growing period (90th Percentile)', 'FontSize', 14);
text(0.5, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 14,'HorizontalAlignment', 'right');
set(gca, 'FontSize', 12);

% Fit linear models and calculate goodness of fit for the peak growing season
[fit_rg, gof_rg] = fit(d_evap_rg, d_gpp_rg, 'poly1');
[fit_ag, gof_ag] = fit(d_evap_ag, d_gpp_ag, 'poly1');

% Calculate slopes
slope_rg = fit_rg.p1;
slope_ag = fit_ag.p1;

% Calculate correlation coefficients (r)
r_rg = sqrt(gof_rg.rsquare);
r_ag = sqrt(gof_ag.rsquare);

% Define x range for plotting fit lines
x_range_rg = linspace(min(d_evap_rg), max(d_evap_rg), 100);
x_range_ag = linspace(min(d_evap_ag), max(d_evap_ag), 100);

% Plot fits and match colors with thicker lines
h7 = plot(x_range_rg, feval(fit_rg, x_range_rg), '--', 'Color', colorRangeland, 'LineWidth', 2);
h8 = plot(x_range_ag, feval(fit_ag, x_range_ag), '--', 'Color', colorCropland, 'LineWidth', 2);

% Add legend for subplot (b) with slopes
legend([h5, h6, h7, h8], {'Rangeland Data', 'Cropland Data', ...
    sprintf('Rangeland: \\it{r}\\rm=%.2f, Slope=%.2f', r_rg, slope_rg), ...
    sprintf('Cropland: \\it{r}\\rm=%.2f, Slope=%.2f', r_ag, slope_ag)}, ...
    'Location', 'northeast');

% Third subplot for WUE bar plot for the whole season
% Third subplot for WUE bar plot for the whole season
subplot(2,2,3); % Position (c)
WUE_mean_rangeland = sum(d_gpp_r) / sum(d_evap_r);
WUE_mean_cropland = sum(d_gpp_a) / sum(d_evap_a);

wue_data_r=d_gpp_r ./ d_evap_r;
wue_data_r(isinf(wue_data_r))=NaN;
SEM_rangeland = nanstd(wue_data_r) / sqrt(sum(~isnan(wue_data_r)));
% SD_rangeland = std(d_gpp_r ./ d_evap_r);
SD_cropland = std(d_gpp_a ./ d_evap_a);
n_rangeland = numel(d_gpp_r);
n_cropland = numel(d_gpp_a);
% SEM_rangeland = SD_rangeland / sqrt(n_rangeland);
SEM_cropland = SD_cropland / sqrt(n_cropland);
WUE_means = [WUE_mean_rangeland, WUE_mean_cropland];
WUE_errors = [SEM_rangeland, SEM_cropland];

bar(1, WUE_means(1), 'BarWidth', 0.5, 'FaceColor', colorRangeland);
hold on;
bar(2, WUE_means(2), 'BarWidth', 0.5, 'FaceColor', colorCropland);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Rangeland', 'Cropland'}, 'FontSize', 12);
ylabel('WUE (gC m^{-2} mm^{-1})', 'FontSize', 14);
ylim([0 5]);
text(0.5, 0.95, '(c)', 'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right');
set(gca, 'FontSize', 12);

errorbar(1, WUE_means(1), WUE_errors(1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
errorbar(2, WUE_means(2), WUE_errors(2), 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Display number of observations
text(1, WUE_means(1) + 0.2, sprintf('N=%d', n_rangeland), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(2, WUE_means(2) + 0.2, sprintf('N=%d', n_cropland), 'HorizontalAlignment', 'center', 'FontSize', 12);

% Perform t-test between rangeland and cropland
[h, p] = ttest2(d_gpp_r ./ d_evap_r, d_gpp_a ./ d_evap_a);

% Add significance stars or 'ns' if not significant
if h == 1
    if p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end
else
    sig = 'ns';
end

text(1.5, max(WUE_means) + 0.3, sig, 'HorizontalAlignment', 'center', 'FontSize', 16);

box on;


% Fourth subplot for WUE bar plot for the peak growing season
subplot(2,2,4); % Position (d)
WUE_mean_rangeland = sum(d_gpp_rg) / sum(d_evap_rg);
WUE_mean_cropland = sum(d_gpp_ag) / sum(d_evap_ag);
SD_rangeland = std(d_gpp_rg ./ d_evap_rg);
SD_cropland = std(d_gpp_ag ./ d_evap_ag);
n_rangeland = numel(d_evap_rg);
n_cropland = numel(d_evap_ag);
SEM_rangeland = SD_rangeland / sqrt(n_rangeland);
SEM_cropland = SD_cropland / sqrt(n_cropland);
WUE_means = [WUE_mean_rangeland, WUE_mean_cropland];
WUE_errors = [SEM_rangeland, SEM_cropland];

bar(1, WUE_means(1), 'BarWidth', 0.5, 'FaceColor', colorRangeland);
hold on;
bar(2, WUE_means(2), 'BarWidth', 0.5, 'FaceColor', colorCropland);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Rangeland', 'Cropland'}, 'FontSize', 12);
ylabel('WUE (gC m^{-2} mm^{-1})', 'FontSize', 14);
ylim([0 5]);
text(0.5, 0.95, '(d)', 'Units', 'normalized', 'FontSize', 14,'HorizontalAlignment', 'right');

errorbar(1, WUE_means(1), WUE_errors(1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
errorbar(2, WUE_means(2), WUE_errors(2), 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Display number of observations
text(1, WUE_means(1) + 0.2, sprintf('N=%d', n_rangeland), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(2, WUE_means(2) + 0.2, sprintf('N=%d', n_cropland), 'HorizontalAlignment', 'center', 'FontSize', 12);

% Perform t-test between rangeland and cropland
[h, p] = ttest2(d_gpp_rg ./ d_evap_rg, d_gpp_ag ./ d_evap_ag);

% Add significance stars
if h == 1
    if p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end
    text(1.5, max(WUE_means) + 0.3, sig, 'HorizontalAlignment', 'center', 'FontSize', 16);
end
box on;

% Force MATLAB to update the figure immediately
drawnow;


%% Supplimentary figures
WUE_r = d_gpp_r./d_evap_r;
WUE_a = d_gpp_a./d_evap_a;
WUE_r(WUE_r<0) = 0;
%% Daily Rangeland GPP, ET and WUE timeseries
% Define custom RGB values for a darker forest green
darker_forest_green_rgb = [0, 100, 0] / 255;

% figure(11);
figure('Position', [100, 100, 1024, 768]); % Adjust the values as needed

subplot(3,1,1);
% scatter(d_date_r,d_gpp_r,'Marker','*','MarkerFaceColor',darker_forest_green_rgb,'MarkerEdgeColor', darker_forest_green_rgb);
plot(d_date_r,d_gpp_r,'Color',darker_forest_green_rgb, 'LineWidth',2);
datetick('x','mmm-yyyy');
xlim([start_date_num, end_date_num]);
ylim([0 15]);

% Manually set x-axis tick positions for the first subplot
x_tick_positions = linspace(start_date_num, end_date_num, 8);

% Hide x-axis tick labels
set(gca, 'XTick', x_tick_positions, 'XTickLabel', {});

ax=gca;
set(ax,'FontSize',14)
leg = findobj(gcf, 'Type', 'Legend');
set(get(ax, 'YLabel'), 'FontSize', 16);
set(leg, 'FontSize', 16);
ylabel('GPP (gC m^{-2} day^{-1})');

subplot(3,1,2);
plot(d_date_r,d_evap_r,'-b','LineWidth',2);
datetick('x','mmm-yyyy');
xlim([start_date_num, end_date_num]);
ylim([0 6]);

% Manually set x-axis tick positions for the first subplot
x_tick_positions = linspace(start_date_num, end_date_num, 8);

% Hide x-axis tick labels
set(gca, 'XTick', x_tick_positions, 'XTickLabel', {});

ax=gca;
set(ax,'FontSize',14)
leg = findobj(gcf, 'Type', 'Legend');
set(get(ax, 'YLabel'), 'FontSize', 16);
set(leg, 'FontSize', 16);
ylabel('ET (mm)');


subplot(3,1,3);
plot(d_date_r,WUE_r,'-r','LineWidth',2);
datetick('x','mmm-yyyy');
ylim([0 15]);

ax=gca;
set(ax,'FontSize',14)
leg = findobj(gcf, 'Type', 'Legend');
set(get(ax, 'YLabel'), 'FontSize', 16);
set(leg, 'FontSize', 16);
ylabel('WUE (gC m^{-2}mm^{-1} day^{-1})');
%% Daily Cropland GPP, ET and WUE timeseries
figure('Position', [100, 100, 1024, 768]); % Adjust the values as needed

subplot(3,1,1);
% scatter(d_date_r,d_gpp_r,'Marker','*','MarkerFaceColor',darker_forest_green_rgb,'MarkerEdgeColor', darker_forest_green_rgb);
plot(d_date_a,d_gpp_a,'Color',darker_forest_green_rgb, 'LineWidth',2);
datetick('x','mmm-yyyy');
xlim([start_date_num, end_date_num]);
ylim([0 15]);

% Manually set x-axis tick positions for the first subplot
x_tick_positions = linspace(start_date_num, end_date_num, 8);

% Hide x-axis tick labels
set(gca, 'XTick', x_tick_positions, 'XTickLabel', {});

ax=gca;
set(ax,'FontSize',14)
leg = findobj(gcf, 'Type', 'Legend');
set(get(ax, 'YLabel'), 'FontSize', 16);
set(leg, 'FontSize', 16);
ylabel('GPP (gC m^{-2} day^{-1})');

subplot(3,1,2);
plot(d_date_a,d_evap_a,'-b','LineWidth',2);
datetick('x','mmm-yyyy');
xlim([start_date_num, end_date_num]);
ylim([0 6]);

% Manually set x-axis tick positions for the first subplot
x_tick_positions = linspace(start_date_num, end_date_num, 8);

% Hide x-axis tick labels
set(gca, 'XTick', x_tick_positions, 'XTickLabel', {});

ax=gca;
set(ax,'FontSize',14)
leg = findobj(gcf, 'Type', 'Legend');
set(get(ax, 'YLabel'), 'FontSize', 16);
set(leg, 'FontSize', 16);
ylabel('ET (mm)');


subplot(3,1,3);
plot(d_date_a,WUE_a,'-r','LineWidth',2);
xlim([start_date_num, end_date_num]);
datetick('x','mmm-yyyy');
ylim([0 15]);

ax=gca;
set(ax,'FontSize',14)
leg = findobj(gcf, 'Type', 'Legend');
set(get(ax, 'YLabel'), 'FontSize', 16);
set(leg, 'FontSize', 16);
ylabel('WUE (gC m^{-2}mm^{-1} day^{-1})');
%%
