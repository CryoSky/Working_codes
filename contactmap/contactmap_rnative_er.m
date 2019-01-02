f1 = load('rnative.dat');
[d,p]=size(f1);
mat = zeros(d);
for i = 1:d
    for j = (i+1):d
        dis1 = (f1(i,j)+f1(j,i))/2;
        if dis1 <= 10 && dis1 > 0
            mat(j,i)=1;
        end
    end
end

f2 = load('contactmap.txt');
[d,p]=size(f2);
for i = 1:d
        if f2(i,5)>=0.4
            mat(f2(i,1),f2(i,2))=f2(i,5);
        end
end



figure(1)

imagesc(mat);
%colormap(flipud(jet))
colormap(jet)
% title('Protein contact map');
% set(gca,'xtick',[0:1:seqlen1])
% set(gca,'ytick',[0:1:seqlen3])
grid on;
grid minor
% set(gca,'ytick',[0:10:seqlen1])
% set(gca,'xtick',[0:50:seqlen3*2])
set(gca,'fontsize',20);
