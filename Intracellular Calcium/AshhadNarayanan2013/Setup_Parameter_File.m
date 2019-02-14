function Setup_Parameter_File(k)

load('X')
fileID = fopen(['para',num2str(k),'.dat'],'w');
fprintf(fileID,'par_gna         = %20.10e\n',0.09*X(k,1)); 
fprintf(fileID,'par_gkdr        = %20.10e\n',0.003*X(k,2));
fprintf(fileID,'par_gka         = %20.10e\n',0.022*X(k,3));
fprintf(fileID,'par_gcal        = %20.10e\n',0.316e-3*X(k,4));
fprintf(fileID,'par_gcat        = %20.10e\n',0.1e-3*X(k,5)); 
fprintf(fileID,'par_gip3        = %20.10e\n',1.85*X(k,6));
fprintf(fileID,'par_vmax        = %20.10e\n',1e-4*X(k,7));
fprintf(fileID,'par_gamma       = %20.10e\n',8*X(k,8));
fprintf(fileID,'par_nmda0       = %20.10e\n',1*X(k,9));
fprintf(fileID,'par_ampa0       = %20.10e\n',1*X(k,10));
fprintf(fileID,'par_PLC         = %20.10e\n',0.83*X(k,11));
fprintf(fileID,'par_G           = %20.10e\n',100*X(k,12));
fprintf(fileID,'par_initmGluR   = %20.10e\n',0.3e-3*X(k,13));
fprintf(fileID,'par_nit         = %d\n',k);
fclose(fileID);
