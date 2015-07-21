%============================= plot_eke.m ================================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in the KE and EKE from .mat files & =%
%=  plot maps of mean KE and mean EKE. For the KE calculation, the       =%
%=  velocities are broken down into u=u_ + u', where u_ is the time mean =%
%=  velocity, and u' is the time-varying perturbation velocity. The KE   =%
%=  is then KE=0.5*(u_^2 + v_^2) and EKE=0.5*(u'^2 + v'^2).              =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/he_data/he/actodd/PLOTS/USeast/validate/KE/';
datadir='/he_data/he/zhigang/Project/USeast/Out_exp';
matdir ='/he_data/he/actodd/DATA/eke/';

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS output');

ncid=netcdf.open([datadir,'22/his_0001.nc'],'nowrite');
 Vstretching=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vstretching'));
 Vtransform =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vtransform' ));
 theta_s    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_s'    ));
 theta_b    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_b'    ));
 hc         =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hc'         ));
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ));
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ));
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ));
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ));
netcdf.close(ncid);

nanmask=mask./mask;
zw=set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,36,5,h,h.*0,0);
h_layer=zw(:,:,2:end)-zw(:,:,1:end-1);

h100mask=mask; h50mask=mask;
h100mask(find(h<100))=0;
h50mask( find(h<50 ))=0;

%=========================================================================%
%=                  Find deepest layers to 50m and 100m                  =%
%=========================================================================%
for i=1:size(h,1)
  for j=1:size(h,2)
    z100ind(i,j)=min(find(zw(i,j,:)>=-100));
     z50ind(i,j)=min(find(zw(i,j,:)>=-50 ));
  end
end


%=========================================================================%
%=                    Load EKE and find upper layer EKE                  =%
%=========================================================================%
eval(['load ',matdir,'mean_JUNE_vels.mat']);

projs=[1,2,3,4,5,6,11,12,14,16,17,18,19,20,21,22];
for ex=length(projs):length(projs)

  for i=1:36
    [uu,vv]=regridromsvels2d(squeeze(squeeze(meanu(ex,i,:,:)))',...
                             squeeze(squeeze(meanv(ex,i,:,:)))');
     mke(i,:,:)=0.5*((uu.*100)'.^2 + (vv.*100)'.^2);			     
  end
 
  %--- Load and find mean of EKE ---%
  eval(['load ',matdir,'eke_JUNE_exp',sprintf('%02i',projs(ex)),'.mat']);
  
  for i=1:size(h,1)
    for j=1:size(h,2)
      %--- Upper 100m mean KE ---%
      if (z100ind(i,j)==37)
        eke100(i,j)=squeeze(squeeze( ke(end,i,j)));
	 ke100(i,j)=squeeze(squeeze(mke(end,i,j)));
      else
        eke100(i,j)=sum(squeeze(squeeze( ke(z100ind(i,j)+1:end,i,j))).*...
	                squeeze(squeeze(abs(h_layer(i,j,z100ind(i,j)+1:end)))))./...
	            squeeze(squeeze(abs(h_layer(i,j,end)-h_layer(i,j,z100ind(i,j)))));

        mke100(i,j)=sum(squeeze(squeeze(mke(z100ind(i,j)+1:end,i,j))).*...
	                squeeze(squeeze(abs(h_layer(i,j,z100ind(i,j)+1:end)))))./...
	            squeeze(squeeze(abs(h_layer(i,j,end)-h_layer(i,j,z100ind(i,j)))));

      end
      
      %--- Now for upper 50m ---%
      if (z50ind(i,j)==37)
        eke50(i,j)=squeeze(squeeze(ke(end,i,j)));
      else
        eke50(i,j)=sum(squeeze(squeeze(ke(z50ind(i,j)+1:end,i,j))).*...
	               squeeze(squeeze(abs(h_layer(i,j,z50ind(i,j)+1:end)))))./...
	           squeeze(squeeze(abs(h_layer(i,j,end)-h_layer(i,j,z50ind(i,j)))));
        mke50(i,j)=sum(squeeze(squeeze(mke(z50ind(i,j)+1:end,i,j))).*...
                       squeeze(squeeze(abs(h_layer(i,j,z50ind(i,j)+1:end)))))./...
	               squeeze(squeeze(abs(h_layer(i,j,end)-h_layer(i,j,z50ind(i,j)))));

      end   
    end
  end
  plot_procedure(lon,lat,mask,log(eke50).*h50mask./h50mask,[5.5 10],...
                 ['Exp ',sprintf('%02i',projs(ex))],['Upper 50m EKE']);
  eval(['export_fig ',plotdir,'EKE50_JUNE_exp',sprintf('%02i',projs(ex)),'.png ',...
       '-png -r150 -painters']);

  plot_procedure(lon,lat,mask,log(mke50).*h50mask./h50mask,[5.5 10],...
                 ['Exp ',sprintf('%02i',projs(ex))],['Upper 50m KE']);
  eval(['export_fig ',plotdir,'KE50_JUNE_exp',sprintf('%02i',projs(ex)),'.png ',...
       '-png -r150 -painters']);

  plot_procedure(lon,lat,mask,log(eke100).*h100mask./h100mask,[5.5 10],...
                ['Exp ',sprintf('%02i',projs(ex))],['Upper 100m EKE']);
  eval(['export_fig ',plotdir,'EKE100_JUNE_exp',sprintf('%02i',projs(ex)),'.png ',...
       '-png -r150 -painters']);

  plot_procedure(lon,lat,mask,log(mke100).*h100mask./h100mask,[5.5 10],...
                ['Exp ',sprintf('%02i',projs(ex))],['Upper 100m KE']);
  eval(['export_fig ',plotdir,'KE100_JUNE_exp',sprintf('%02i',projs(ex)),'.png ',...
       '-png -r150 -painters']);
end
return;

%=========================================================================%
%=                      Load mean vels calculate KE                      =%
%=========================================================================%
eval(['load ',matdir,'mean_JUNE_vels.mat']);

for i=1:length(projs)
    
  %--- First, mean Kinetic Energy ---%
  [uu,vv]=regridromsvels2d(squeeze(meanu(i,:,:))',squeeze(meanv(i,:,:))');
  uu=uu'.*100; vv=vv'.*100;
  ke = 0.5*(uu.^2 + vv.^2);
  
  plot_procedure(lon,lat,mask,log(ke),[4.5 7.5],['Exp ',sprintf('%02i',projs(i))],...
                                          ['Upper 50m KE']);
  eval(['export_fig ',plotdir,'meanKE_JUNE_exp',sprintf('%01i',projs(i)),'.png ',...
       '-png -r150 -painters']);

  %--- Load and find mean of EKE ---%
  eval(['load ',matdir,'eke_exp',sprintf('%1i',i),'.mat']);
  eke=squeeze(nanmean(ke,1));
  plot_procedure(lon,lat,mask,log(eke),[4.5 7.5],['TNU=40'],...%['Exp ',sprintf('%02i',i)],...
                                           [labs{i},' Eddy KE']);
  eval(['export_fig ',plotdir,'EKE_JUNE_exp',sprintf('%01i',i),'.png ',...
       '-png -r150 -painters']);
       
end
  

