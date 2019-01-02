% Written by Shikai Jin on 2018-03-10, latest modified on 2018-03-12
% Python version to generate contact matrix for one pdb file
% Run Contactmap_in_python.py in advance to generate contact matrix

contactfile = '.\contact_matrix.txt';
[A,delimit1] = importdata(contactfile);
[rows,columns] = size(A);
for i = 1 : rows
    for j = 1 : columns
        if A(i,j) <= 12
            rectangle('Position',[i,j,1,1],'Facecolor',[0,0,1]);
        end
    end
end

axis([0,rows,0,columns]);
fsize = 14
xlabel('Residue Index', 'fontsize', fsize); ylabel('Residue Index', 'fontsize', fsize);
title('Protein Distance map');

% for i = 1 : 1596
%    if A(i,3) == 1
%      rectangle('Position',[A(i,1),A(i,2),1,1],'Facecolor',[0,0,A(i,3)])
%    end
% end
% for i = 1 : 2926
%    if B(i,3) == 1 
%     rectangle('Position',[B(i,1),B(i,2),1,1],'Facecolor',[B(i,3),0,0])
% 
%    end
% end
