function Contact_map(pdb1)
% Draw contact map for one pdb file, written by Shikai Jin on 2018-03-12
% Latest modified on 2018-03-12
% Modified from Mingchen's version
% Run in Contact_map('1mba.pdb')

struct1 = pdbread(pdb1);

x1 = [struct1.Model.Atom.X];
y1 = [struct1.Model.Atom.Y];
z1 = [struct1.Model.Atom.Z];

pos1 = find(strcmp({struct1.Model.Atom.AtomName},'CA'));

ln = numel(pos1);

figure(1)
mat = zeros(ln); % returns an n-by-n matrix of zeros.
for i = 1:ln
    index1 = pos1(i);
    for j = (i+1):ln
        index2 = pos1(j);
        dis1 = (x1(index1)-x1(index2)).^2 + (y1(index1)-y1(index2)).^2 +(z1(index1)-z1(index2)).^2;
            mat(i,j) = dis1.^0.5;
            mat(j,i) = dis1.^0.5;
    end
end

imagesc(mat);colormap('jet');colorbar % https://cn.mathworks.com/help/matlab/ref/colormap.html#inputarg_map
%set(gca,'YDir','normal') 
% imagesec will change y axis to reverse
% direction, use above line to return it
fsize = 14;
xlabel('Residue Index', 'fontsize', fsize); ylabel('Residue Index', 'fontsize', fsize); % default value is 11
title('Protein Distance map');
end

%function Distance_PDB(pdb1, cutoff)
%disp('inputpdb1 sequencelength label_separation');
%end
