clear;
clc;

%% Main program

x = 1:20
dir = 'E:\0127\1MBA\'
prefix = '1mba-3pt8-HO+ER'
id = 'T0766'


y1 = importdata([dir,prefix,'_0.25_maxQ.txt'])
y1 = sort(y1)
y2 = importdata([dir,prefix,'_0.5_maxQ.txt'])
y2 = sort(y2)
y3 = importdata([dir,prefix,'_1.0_maxQ.txt'])
y3 = sort(y3)
y4 = importdata([dir,prefix,'_1.5_maxQ.txt'])
y4 = sort(y4)
y5 = importdata([dir,prefix,'_2.0_maxQ.txt'])
y5 = sort(y5)

hold off

plot(x,y1,'d-','Color','r','MarkerFaceColor','r','MarkerSize',5)
hold on
plot(x,y2,'d-','Color','b','MarkerFaceColor','b','MarkerSize',5)
hold on
plot(x,y3,'d-','Color','k','MarkerFaceColor','k','MarkerSize',5)
hold on
plot(x,y4,'d-','Color','m','MarkerFaceColor','m','MarkerSize',5)
hold on
plot(x,y5,'d-','Color',[0.4 0.8 0.3],'MarkerFaceColor',[0.4 0.8 0.3],'MarkerSize',5) % Yellow color looks too shine in white background, change another mid green by RGB value see https://cn.mathworks.com/help/matlab/ref/plot.html

title([id, ' Template free-annealing profile'])
xlabel('Annealing Index')
ylabel('Qw-Best')
legend('0.25','0.5','1.0','1.5','2.0','Location','southeast')
saveas(gcf,['E:\',prefix,'-annealing-tempfree.png'])

