function ExportHomeFEToMedit(filename, mesh)
% Export a mesh (as defined in HomeFE code) in a Medit file, readable by
% Freefem++.
%
%  syntax: ExportHomeFEToMedit(filename, mesh)
%
%  filename: output file name string
%
%  mesh: a mesh cell. Must contain the cells "X" (with the coordinates) and
%        "T" (with the triangles)
%
% Prepare dimensions
dimension = size(mesh.X,2);
nbOfVertices = size(mesh.X,1);
nbOfTriangles = size(mesh.T,1);

% Open file and print header info
outputF = fopen(filename,'w');

fprintf(outputF,'MeshVersionFormatted 2\n');

fprintf(outputF,'Dimension\n');
fprintf(outputF,'%i\n',dimension);

% Print vertices
fprintf(outputF,'Vertices\n');
fprintf(outputF,'%i\n',nbOfVertices);
for iii=1:nbOfVertices
    for jjj=1:dimension
        fprintf(outputF,'   %.16f',mesh.X(iii,jjj));
    end
    fprintf(outputF,'   1\n');
end

% Print triangles
fprintf(outputF,'Triangles\n');
fprintf(outputF,'%i\n',nbOfTriangles);
for iii=1:nbOfTriangles
    fprintf(outputF,'   %i   %i   %i   1\n',mesh.T(iii,:));
end

fclose(outputF);

end

