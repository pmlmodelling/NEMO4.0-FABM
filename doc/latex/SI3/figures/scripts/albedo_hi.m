clear all, close all

% plots SI3 albedo dependence versus ice thickness

z1_href_pnd = 1. / 0.05
z1_c1 = 1. / ( log(1.5) - log(0.05) );
z1_c2 = 1. / 0.05;
z1_c3 = 1. / 0.02
z1_c4 = 1. / 0.03

rn_alb_idry = 0.60;
rn_alb_imlt = 0.50;
rn_alb_sdry = 0.85
rn_alb_smlt = 0.75
rn_alb_oce  = 0.06;
rn_alb_dpnd = 0.27;

h_i = 0:0.001:3;

%--- Bare ice albedo (for hi > 150cm)
zalb_ice    = h_i;
zalb_ice(:) = rn_alb_idry;

%--- Bare ice albedo (for hi < 150cm)
zaddr = find( ( 0.05 < h_i ) & ( h_i <= 1.5 ) ); % 5cm < hi < 150cm
zalb_ice(zaddr) = zalb_ice(zaddr) + ( 0.18 - zalb_ice(zaddr) ) * z1_c1 .* ( log(1.5) - log(h_i(zaddr)) );
zaddr = find( h_i <= 0.05 );
zalb_ice(zaddr) = rn_alb_oce  + ( 0.18 - rn_alb_oce ) * z1_c2 * h_i(zaddr);

figure; 


subplot(1,3,1)

plot(h_i,zalb_ice,'k','LineWidth', 3)
xlabel('h_i(m)')
ylabel('\alpha')
ylim([0.1 0.9])
set(gca,'fontsize', 16)
set(gca, 'FontName', 'Helvetica LT Std')



%--- Snow covered albedo

h_s = 0:0.001:0.6;
zalb_ice    = 0.6
zalb_snw    = rn_alb_sdry - ( rn_alb_sdry - zalb_ice ) * exp( -h_s * z1_c3 )

subplot(1,3,2)
plot(h_s,zalb_snw,'k','LineWidth', 3)
xlabel('h_s(m)')
ylabel('\alpha')
ylim([0.1 0.9])
set(gca,'fontsize', 16)
set(gca, 'FontName', 'Helvetica LT Std')

%--- Ponded ice albedo
zalb_ice    = 0.5
h_pnd = 0:0.001:0.6;
zalb_pnd = rn_alb_dpnd - ( rn_alb_dpnd - zalb_ice ) * exp( -h_pnd * z1_href_pnd )

subplot(1,3,3)
plot(h_pnd,zalb_pnd,'k','LineWidth', 3)
xlabel('h_{pnd}(m)')
ylabel('\alpha')
ylim([0.1 0.9])
set(gca,'fontsize', 16)
set(gca, 'FontName', 'Helvetica LT Std')


