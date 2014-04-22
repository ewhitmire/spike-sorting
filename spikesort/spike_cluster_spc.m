function clusters=spike_cluster_spc(features)

    fname='param';
    
    f = features';
    save('data.in', 'f', '-ascii', '-tabs');
    
    % DELETE PREVIOUS FILES
    fileexist = exist(['data.out.dg_01.lab'],'file');
    if(fileexist~=0)
        delete(['data.out.dg_01.lab']);
        delete(['data.out.dg_01']);
    end

    fid=fopen(sprintf('%s.run',fname),'wt');
    fprintf(fid,'NumberOfPoints: %s\n',num2str(size(features, 2)));
    fprintf(fid,'DataFile: %s\n','data.in');
    fprintf(fid,'OutFile: %s\n','data.out');
    fprintf(fid,'Dimensions: %s\n',num2str(size(features, 1)));
    fprintf(fid,'MinTemp: %s\n',num2str(0));
    fprintf(fid,'MaxTemp: %s\n',num2str(.201)); %changed, orig. value: 0.201
    fprintf(fid,'TempStep: %s\n',num2str(0.01));
    fprintf(fid,'SWCycles: %s\n',num2str(100));
    fprintf(fid,'KNearestNeighbours: %s\n',num2str(11)); % orig. value: 11
    fprintf(fid,'MSTree|\n');
    fprintf(fid,'DirectedGrowth|\n');
    fprintf(fid,'SaveSuscept|\n');
    fprintf(fid,'WriteLables|\n');
    fprintf(fid,'WriteCorFile~\n');
    if 0
        fprintf(fid,'ForceRandomSeed: %s\n',num2str(handles.par.randomseed));
    end    
    fclose(fid);

    dos(sprintf('Cluster.exe %s.run',fname));

    clu=load(['data.out.dg_01.lab']);
    tree=load(['data.out.dg_01']); 
    delete(sprintf('%s.run',fname));
    delete *.mag
    delete *.edges
    delete *.param
    delete('data.in'); 
    
    temp_index = 2;
    data = clu(temp_index, 3:end);
    non_clustered = find(data >= 3);
    data(non_clustered) = 3;
    clusters = data;
end
