function Contact_pdb_pdb(pdb1,pdb2,prefix)
% Draw contact map for two pdb files, written by Shikai Jin on 2018-03-12
% Latest modified on 2018-03-12
% Modified from Mingchen's version
% Run in Contact_map('1mba.pdb')

struct1 = pdbread(pdb1);

x1 = [struct1.Model.Atom.X];
y1 = [struct1.Model.Atom.Y];
z1 = [struct1.Model.Atom.Z];

pos1 = find(strcmp({struct1.Model.Atom.AtomName},'CA'));

struct2 = pdbread(pdb2);
 
 x2 = [struct2.Model.Atom.X];
 y2 = [struct2.Model.Atom.Y];
 z2 = [struct2.Model.Atom.Z];

 pos2 = find(strcmp({struct2.Model.Atom.AtomName},'CA'));

 
 if numel(pos1)>= numel(pos2)
    size = numel(pos1);
 else
    size = numel(pos2);
 end

figure(1)

mat = zeros(size); % returns an n-by-n matrix of zeros.

ln = numel(pos1)
for i = 1:ln
    index1 = pos1(i);
    for j = (i+1):ln
        index2 = pos1(j);
        dis1 = (x1(index1)-x1(index2)).^2 + (y1(index1)-y1(index2)).^2 +(z1(index1)-z1(index2)).^2;
            mat(i,j) = dis1.^0.5;
            %mat(j,i) = dis1.^0.5;
    end
end

ln = numel(pos2)
for i = 1:ln
    index1 = pos2(i);
    for j = (i+1):ln
        index2 = pos2(j);
        dis2 = (x2(index1)-x2(index2)).^2 + (y2(index1)-y2(index2)).^2 +(z2(index1)-z2(index2)).^2;
            mat(j,i) = dis2.^0.5;
            %mat(j,i) = dis1.^0.5;
    end
end



imagesc(mat);colormap('jet');colorbar % https://cn.mathworks.com/help/matlab/ref/colormap.html#inputarg_map
%set(gca,'YDir','normal') 
% imagesec will change y axis to reverse
% direction, use above line to return it
fsize = 14;
xlabel('Native Residue Index', 'fontsize', fsize); ylabel('Template Residue Index', 'fontsize', fsize); % default value is 11
title([prefix, ' Protein Distance map']);
saveas(gcf,['E:\',prefix,'.png'])
end
%function Distance_PDB(pdb1, cutoff)
%disp('inputpdb1 sequencelength label_separation');
%end
