clear all
close all

% tollist = exp(linspace(log(1e-6),log(1e-2),30));
tollist = exp(linspace(log(1e-14),log(1e-6),30));
% logspace(1e-6, 1e-2, 30)


N = 16;

difflist = zeros(size(tollist));
timelist = zeros(size(tollist)); 
iterlist = zeros(size(tollist)); 
for ii = 1:length(tollist)
    disp(tollist(ii))
    [totaldiff, elapsed_time, iter] = funNS_solver_tol(tollist(ii), N);
    difflist(ii) = totaldiff; 
    timelist(ii) = elapsed_time; 
    iterlist(ii) = iter; 
end
% 
% figure
% plot(tollist, difflist, 'LineWidth', 2)


width = 15;
height = 12;

figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto')
semilogx(tollist, difflist, 'LineWidth', 2)

figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto') 
semilogx(tollist32, difflist32, 'LineWidth', 2)




% 
% figure
% plot(tollist, timelist)
% 
% 
% semilogx(tollist, timelist)

figure 
plot(tollist, iterlist, 'LineWidth', 2)

figure
plot(tollist32, iterlist32, 'LineWidth', 2)
% 
% figure 
% semilogx(tollist, iterlist, 'LineWidth', 2)




%%
width = 22;
height = 12;
figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto') 
subplot(2,2,1)
semilogx(tollist, difflist, 'LineWidth', 2)
grid on 
legend('$N = 16$','Location', 'NorthWest', 'Interpreter', 'LaTeX')
xlabel('Stop tolerance $\epsilon$', 'Interpreter', 'LaTeX')
ylabel('Error', 'Interpreter', 'LaTeX');
set(gca,'FontSize',10)

subplot(2,2,3)
semilogx(tollist32, difflist32, 'LineWidth', 2)
grid on 
legend('$N = 32$','Location', 'NorthWest', 'Interpreter', 'LaTeX')
xlabel('Stop tolerance $\epsilon$', 'Interpreter', 'LaTeX')
ylabel('Error', 'Interpreter', 'LaTeX');
set(gca,'FontSize',10)

subplot(2,2,2)
plot(tollist, iterlist, 'LineWidth', 2)
grid on 
legend('$N = 16$','Location', 'NorthEast', 'Interpreter', 'LaTeX')
xlabel('Stop tolerance $\epsilon$', 'Interpreter', 'LaTeX')
ylabel('n iterations', 'Interpreter', 'LaTeX');
set(gca,'FontSize',10)

subplot(2,2,4)
plot(tollist32, iterlist32, 'LineWidth', 2)
grid on 
legend('$N = 32$','Location', 'NorthEast', 'Interpreter', 'LaTeX')
xlabel('Stop tolerance $\epsilon$', 'Interpreter', 'LaTeX')
ylabel('n iterations', 'Interpreter', 'LaTeX');
set(gca,'FontSize',10)

print -depsc2 figures/stoptol

















