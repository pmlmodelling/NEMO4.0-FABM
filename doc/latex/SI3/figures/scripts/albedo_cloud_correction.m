clear all; close all;

% Script to plot the cloud fraction dependency of albedo correction

%------------------
% Sea Ice
%------------------

% overcast albedo
alpha_os = 0:0.001:1.;

% SI3 correction
f = ( -0.1010 * alpha_os.^2 + 0.1933 * alpha_os - 0.0148 );

% clear sky sea ice albedo
alpha_cs = alpha_os - f;

% Grenfell and Perovich, Table 3 data
alpha_GP04 = [ 0.114 0.151 0.141 0.266 0.367 0.309 0.332 0.492 0.678 0.928 ];
f_GP04     = [ -0.010 -0.022 0.012 0.027 0.020 0.036 0.037 0.056 0.071 0.077 ];

%--------------
% Ocean
%--------------
alpha_oce_cs = 0.05 / ( 1.1 * 0.40^1.4 + 0.15 ) % --- Briegleb and Ramanathan 1972
alpha_oce_os = 0.06 %--- Payne 1972
f_oce = alpha_oce_os - alpha_oce_cs

%--------------
% Plot
%--------------
plot(alpha_os,f,'k'); hold on
zaddr = find( ( alpha_os >= alpha_GP04(3)  )   & ...
              ( alpha_os <= alpha_GP04(10) ) );
plot(alpha_os(zaddr),f(zaddr),'k', 'LineWidth', 3)
plot(alpha_GP04(3:10),f_GP04(3:10),'ksq','MarkerFaceColor', 'k')
plot(alpha_GP04(5),f_GP04(5),'rsq','MarkerFaceColor','r')
text(alpha_GP04(5)+0.02,f_GP04(5),'dirty ice')
plot(alpha_GP04(1:2),f_GP04(1:2),'rsq','MarkerFaceColor','r')
text(alpha_GP04(1)+0.02,f_GP04(1),'wet tundra')
text(alpha_GP04(2)+0.02,f_GP04(2),'dry tundra')
plot(alpha_oce_os,f_oce,'bo','MarkerFaceColor', 'b')
text(alpha_oce_os+0.02,f_oce,'ocean')
plot( [ alpha_GP04(3) alpha_GP04(3) ], [ -0.06 0.08 ], 'k--' )
plot( [ alpha_GP04(10) alpha_GP04(10) ], [ -0.06 0.08 ], 'k--' )

xlabel('\alpha^{os}')
ylabel('\alpha^{os}-\alpha^{cs}')

set(gca,'fontsize', 16)
set(gca, 'FontName', 'Helvetica LT Std')
