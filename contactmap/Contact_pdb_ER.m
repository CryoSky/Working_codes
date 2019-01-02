function Contact_pdb_ER(pdb1,er_len)

struct1 = pdbread(pdb1);

x1 = [struct1.Model.Atom.X];
y1 = [struct1.Model.Atom.Y];
z1 = [struct1.Model.Atom.Z];
pos1 = find(strcmp({struct1.Model.Atom.AtomName},'CA'));
ln = numel(pos1);

if er_len >= ln
    ln = er_len;
end

figure(1) 
mat = zeros(ln); % returns an n-by-n matrix of zeros.
for i = 1:ln
    index1 = pos1(i);
    for j = (i+1):ln
        index2 = pos1(j);
        dis1 = ((x1(index1)-x1(index2)).^2 + (y1(index1)-y1(index2)).^2 +(z1(index1)-z1(index2)).^2).^0.5;
        if dis1 <= 10
            mat(j,i) = 1;
        else
            mat(j,i) = 0;
        end
                
         
    end
end

f = load('contactmap.txt');
[d,p]=size(f);
for i = 1:d
        if f(i,5)>=0.4
            mat(f(i,1),f(i,2))=f(i,5);
        end
end

imagesc(mat);colormap(jet);colorbar

grid on;
grid minor
% set(gca,'ytick',[0:10:seqlen1])
% set(gca,'xtick',[0:50:seqlen3*2])
set(gca,'fontsize',20);