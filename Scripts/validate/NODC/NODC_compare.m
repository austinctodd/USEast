%=========================== NODC_compare.m ==============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program compare results from a the US East Coast water age ROMS =%
%=  simulation to a whole set of data from the NODC, including CTD, XBT, =%
%=  and moored buoy observations.                                        =%
%=                                                                       =%
%=========================================================================%

addpath('/home/actodd/MYMATLAB/');

%=========================================================================%
%=                    Load in each individual dataset                    =%
%=========================================================================%
display('Loading CTD data');
load /home/actodd/ROMS-utils/USeast-age/validate/CTD_ROMS_data.mat


%=========================================================================%
%=                  Calculate statistics for each cruise                 =%
%=========================================================================%
for i=1:length(CTD.cruisetemp)
  xx=[]; yy=[];
  for j=1:length(CTD.cruisetemp{i})
    
    %--- Get cruise Temperature Profile  
    t=roms.cruisetemp{i}{j}; z1=-roms.cruisez{ i}{j};
    y= CTD.cruisetemp{i}{j}; z2=  CTD.cruisetz{i}{j};
    
    %--- Interpolate ROMS values to OBS depths ---%
    if (length(z1)==length(find(isnan(z1)==1)))
        x=z2{:}.*nan;
    else
      x=interp1(z1,t,z2{:});
    end
    
    %--- Only use values without nan for calculations ---%
    a=find(isnan(y)==0); 
    if (length(a)>0) 
      x=x(a); y=y(a); 
    else
      x=[]; y=[];    
    end
    a=find(isnan(x)==0);
    if (length(a)>0) 
      x=x(a); y=y(a); 
    else
      x=[]; y=[];    
    end
    xx=[xx;x]; yy=[yy;y];  
  end
  p=polyfit(xx,yy,1);
  yfit=polyval(p,xx);
  r=cov(xx,yy)./(std(xx).*std(yy));
  
  CTD.cruise_r(  i)=r(1,2);
  CTD.cruise_std(i)=std(yy);
 roms.cruise_std(i)=std(xx); 
  CTD.cruise_reg(i)=std(xx)./std(yy);
  
end



