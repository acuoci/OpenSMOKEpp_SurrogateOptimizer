function drawPlots(target,database,solution)
% Plots
plots.spx = 2; % # of SubPlots (horizontal / vertical)
if target.doIDTLBV == true
    plots.spx = 2; plots.spy = 2;
else
    plots.spx = 2; plots.spy = 1;
end
plots.idx = 1;
figure;
% Plot DC
subplot(plots.spy,plots.spx,plots.idx), plots.idx = plots.idx+1;
    % Target DC
plot (target.vol, target.Td, 'bo', 'LineWidth', 1.5,'DisplayName','target'); hold on;
    % Surrogate DC
plot (solution.vol,solution.Td, 'r-', 'LineWidth', 1.5 ,'DisplayName','Current x');
    % Boiling points of pure compounds
% Tpures = zeros(1,length(database.name));
% colors = lines(length(database.name));
% for i = 1:length(database.name)
% xp = zeros(1,length(database.name));    xp(i) = 1;
% Tpures(i) = Flash(Td(1), target.P, xp, 0., database.vpType, database.vpCoeffs);
% yline(Tpures(i),'Color',colors(i,:)); % ,'DisplayName',database.name{i}
% end
legend({'target','solution'},'Location', 'northwest');
title('Distillation curve');
xlabel('volume recovery (%)');
ylabel('temperature (°C)');

if target.doIDTLBV == true
% Plot IDT curve
subplot(plots.spy,plots.spx,plots.idx); plots.idx=plots.idx+1;
plot(1000./target.TIdt, target.Idt, 'bo','LineWidth',1.5);
hold on;
plot(1000./target.TIdt, solution.idt, 'r-','LineWidth',1.5)
set(gca, 'YScale', 'log');
legend({'target', 'solution'},'Location', 'northwest');
xlim([min(1000./target.TIdt),max(1000./target.TIdt)]); ylim([1e-5 1e-1]);
xlabel('1000/T [K]'); ylabel('IDT [s]');  title('Ignition delay times');

% Plot LBV curve
subplot(plots.spy,plots.spx,plots.idx); plots.idx=plots.idx+1;
    % Plot target LBV
plot(target.phiLBV, target.lbv, 'bo', 'LineWidth', 1.5);
hold on;
    % Plot surrogate LBV
plot(solution.phiLBV, solution.lbv, 'r-', 'LineWidth', 1.5);
xlim([min(target.phiLBV),max(target.phiLBV)]);
a = axis; a(3) = 0; axis(a); clear a;
title('Laminar flame speeds'); xlabel('Equivalence ratio'); ylabel('s_L [cm/s]'); 
    % Plot LBVs of pure compounds
%for i = 1:6, plot(lbvmodel.phipures, lbvmodel.sLpures(:,i)), hold on, end
l=legend('target', 'solution'); l.Location = 'southeast';
end


% Plot others
subplot(plots.spy,plots.spx,plots.idx); plots.idx=plots.idx+1;
scalarprop = {"","H/C","MW","CN","TSI","\mu","YSI","\rho"};
for i = 1:7
    switch(i)
        case 1
            handleBar=barh(i,solution.HC /target.HC, 'BarWidth', 0.8); hold on;
            lbl = sprintf(" %4.2f (target %4.2f)", solution.HC, target.HC);
            optimized = target.optimizedProperties(2);
        case 2
            handleBar=barh(i,solution.MW /target.MW, 'BarWidth', 0.8);
            lbl = sprintf(" %4.2f (target %4.2f)", solution.MW, target.MW);
            optimized = target.optimizedProperties(1);
        case 3
            handleBar=barh(i,solution.CN /target.CN, 'BarWidth', 0.8);
            lbl = sprintf(" %4.2f (target %4.2f)", solution.CN, target.CN);
            optimized = target.optimizedProperties(3);
        case 4
            handleBar=barh(i,solution.TSI/target.TSI, 'BarWidth', 0.8);
            lbl = sprintf(" %4.2f (target %4.2f)", solution.TSI, target.TSI);
            optimized = target.optimizedProperties(4);
        case 5
            handleBar=barh(i,solution.mu /target.mu, 'BarWidth', 0.8);
            lbl = sprintf(" %4.2f (target %4.2f)", solution.mu, target.mu);
            optimized = target.optimizedProperties(5);
        case 6
            handleBar=barh(i,solution.YSI /target.YSI, 'BarWidth', 0.8);
            lbl = sprintf(" %4.2f (target %4.2f)", solution.YSI, target.YSI);
            optimized = target.optimizedProperties(6);
        case 7
            handleBar=barh(i,solution.rho/target.rho, 'BarWidth', 0.8);
            lbl = sprintf(" %4.2f (target %4.2f)", solution.rho, target.rho);
            optimized = target.optimizedProperties(8);
    end
    if optimized, set(handleBar, 'FaceColor', [0.8,0.8,0.8], 'EdgeColor','none');
    else,         set(handleBar, 'FaceColor', [0.9,0.9,0.9], 'EdgeColor','none');
    end
    t = text(0,i,lbl);
end
set(gca, 'YTick', (0:length(scalarprop)), 'YTickLabel', scalarprop, 'YDir','reverse')
xlabel('property_{surrogate}/property_{target}'); title('Scalar Properties');
xlim([0 1.5]); xline(1,'k');
end