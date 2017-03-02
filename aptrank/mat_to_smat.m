function mat_to_smat(A,filename)
[i,j,v] = find(A);
M = [i';j';v'];
M(1,:) = M(1,:)-1;
M(2,:) = M(2,:)-1;
fileID = fopen(strcat(filename,'.smat'),'w');
fprintf(fileID,'%3d %3d %3d\n',size(A,1),size(A,2),length(v));
fprintf(fileID,'%3d %3d %3d\n',M);
fclose(fileID);
end