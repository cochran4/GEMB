function GenerateHocFiles

% Loop through files
for i=1:320
    fileID = fopen(['ParameterInstance',num2str(i),'.hoc'],'w');
    
    fprintf(fileID,'// Load parameters\n');
    fprintf(fileID,'objref para_file\n');
    fprintf(fileID,'para_file = new File()\n');
    fprintf(fileID,['xopen("para',num2str(i),'.dat")\n']);
    fprintf(fileID,'para_file.close()\n\n');
    
    fprintf(fileID,['// Mark File to Write (to signify, time to move to next simulation)\n', ...
    'objref sav_cai\n',...
    'strdef file_name\n',...
    'sprint(file_name,"ave_cai%%d.dat",par_nit)\n',... 
    'sav_cai = new File()\n',...
    'sav_cai.wopen(file_name)\n',...
    'sav_cai.close()\n\n']);


    fprintf(fileID,'load_file("Fig6C-F_Edited.hoc")\n\n');

    fprintf(fileID,['// Record variables\n',...
    'objref rec_t, rec_cai\n',...
    'rec_t   = new Vector()\n',...
    'rec_cai = new Vector()\n',...
    'tLower  = 0*1000\n',...
    'tUpper  = 3*1000-1\n',...
    'rec_t   = new Vector(tUpper-tLower+1)\n',...
    'rec_t.indgen(tLower,tUpper,1)\n',...
    'rec_cai.record(&apical[112].cai(0.5),rec_t)\n\n']); 

    fprintf(fileID,['// Run Simulation\n',...
    'run()\n\n']); 

    fprintf(fileID,['// Post-process data (save average cai between time values)\n',...
    'objref sav_cai\n',...
    'strdef file_name\n',...
    'sprint(file_name,"ave_cai%%d.dat",par_nit)\n',... 
    'sav_cai = new File()\n',...
    'sav_cai.wopen(file_name)\n',...
    'sav_cai.printf("%%-20.10e",rec_cai.mean())\n',...
    'sav_cai.close()\n']);
    
    fclose(fileID);
end