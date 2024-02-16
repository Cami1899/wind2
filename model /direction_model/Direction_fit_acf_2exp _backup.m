function [fitresult, gof] = fit_acf_2exp_backup(lag, acf)


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( lag, acf );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.148248723149503 -0.119737271211406 0.765744689976673 -0.00121014417781055];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Create a figure for the plots.
figure
h = plot( fitresult, xData, yData );
legend( h, 'acf vs. lag', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'lag', 'Interpreter', 'none' );
ylabel( 'acf', 'Interpreter', 'none' );
title('Fitted ACF','FontSize', 16)
grid on

figure( 'Name', 'untitled fit 1' );

% % Plot fit with data.
% subplot( 2, 1, 1 );
% h = plot( fitresult, xData, yData );
% legend( h, 'acf vs. lag', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'lag', 'Interpreter', 'none' );
% ylabel( 'acf', 'Interpreter', 'none' );
% title('Fitted ACF','FontSize', 16)
% grid on
% 
% % Plot residuals.
% subplot( 2, 1, 2 );
% h = plot( fitresult, xData, yData, 'residuals' );
% legend( h, 'untitled fit 1 - residuals', 'Zero Line', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'lag', 'Interpreter', 'none' );
% ylabel( 'residuals acf', 'Interpreter', 'none' );
% title('Residuals of ACF','FontSize', 16)
% grid on


