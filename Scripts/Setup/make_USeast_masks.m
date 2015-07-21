%==============================================================================%
%=                                                                            =%
%=   created by Austin C. Todd, NCSU (2013)                                   =%
%=                                                                            =%
%=  Program reads in contents from an existing netCDF grid file for the US    =%
%=  East Coast domain, and outputs the information into a new grid file that  =%
%=  will be used for experiments testing the age and residence time of water  =%
%=  in different basins.  After the data is all read in, the program will     =%
%=  make the masks in the same way they are created using
%=  create_useast_tracermask.m, then masks are output to the new grid file.   =%
%=                                                                            =%
%==============================================================================%

addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/m_map/'));
addpath(genpath('/he_data/he/zhigang/Matlab_Codes/mexcdf2009'));
format long; format compact;

%==============================================================================%
%                  Load in all and write out the grid variables               =%
%==============================================================================%
infile ='/share/actodd/ROMS-exp/frc/USeast-newbathy.nc';
%outfile='/home/actodd/WORK/USeast/Data/USeast-grid-age.nc';
h    =nc_varget(infile,'h'       );
lon  =nc_varget(infile,'lon_rho' );
lat  =nc_varget(infile,'lat_rho' );
hmask=nc_varget(infile,'mask_rho');

%varnames={'depthmin';'depthmax';'f';'pm';'pn';'dndx';'dmde';'x_rho';...
%          'y_rho';'x_psi';'y_psi';'x_u';'y_u';'x_v';'y_v';'lat_psi';'lon_psi';...
%          'lat_u';'lon_u';'lat_v';'lon_v';'mask_u';'mask_v';'mask_psi';...
%          'angle';'xl';'el';'spherical';'hformer'};
%for i=1:length(varnames)
%  var=nc_varget( infile,char(varnames(i)));
%      nc_varput(outfile,char(varnames(i)),var);
%end

%var=nc_varget( infile,'hraw');
%    nc_varput(outfile,'hraw',var,[0 0 0],[2 482 402])

%nc_varput(outfile,'lon_rho',lon); nc_varput(outfile,'h'       ,h    );
%nc_varput(outfile,'lat_rho',lat); nc_varput(outfile,'mask_rho',hmask);

%clear var

%==============================================================================%
%                             First guess at masks                            =%
%==============================================================================%
mask200 =hmask; mask200( find(h>200 ))=0;
mask2500=hmask; mask2500(find(h<2500))=0;

%maskMISS=hmask;
%maskMISS(  1:146,101:end)=0; maskMISS(   :   ,186:end)=0;
%maskMISS(    147,129:end)=0; maskMISS(    148,132:end)=0;
%maskMISS(    149,134:end)=0; maskMISS(    150,136:end)=0;
%maskMISS(    151,138:end)=0; maskMISS(    152,141:end)=0;
%maskMISS(    153,143:end)=0; maskMISS(    154,144:end)=0;
%maskMISS(154:162,152:end)=0; maskMISS(234:end,176:end)=0;

maskMISS=hmask;

%--- Area southeast of Yucatan and East of FL Straits ---%
maskMISS(  1:146,101:end)=0; maskMISS(   :   ,189:end)=0;
x=127:145; y=round((1/3).*(x-127)+147);
for i=1:length(x)
  maskMISS(y(1):y(i),x(i):end)=0;
end

%--- FL Bay ---%
maskMISS(190:191,189)=1; maskMISS(191,190)=1;

%--- Cuba Bay ---%
maskMISS(154:162,156:188)=0;

%--- Jacksonville Area --%
maskMISS(235:276,182:188)=0;

%--- For scope, make Key West 1 ---%
scopeMISS=maskMISS;
scopeMISS(185:187,174:188)=1;


%==============================================================================%
%             Fine tune shelf mask to cut off areas of disinterest            =%
%==============================================================================%

%--- Carribean and S America ---%
mask200(  1:250,200:402)=0;
mask200(172:179,192:201)=0;

mask200(100:170,140:201)=0;
mask200( 88:97 ,188:193)=0;
mask200(  5:30 ,160:200)=0;
mask200( 85:143,105:144)=0;

%--- Add some shelf waters in Western Gulf of Mexico ---%
mask200(118:119,39:41)=1; 
mask200(132:133,27)=1; 
mask200(131,28)=1;
mask200(151,19)=1;

%-- Extend shelf waters in W GOM ---%
mask200(118,42:43)=1;mask200(123,33)=1; mask200(123:124,32)=1;
mask200(127:128,30)=1;   mask200(129:131,29)=1;
mask200(132,28)=1; mask200(133,28)=1;
mask200(134,27)=1; mask200(135,26)=1;
mask200(148:150,20)=1; mask200(152,19)=1;
mask200(153,18)=1; 
mask200(120,37:39)=1; mask200(122,34)=1;
mask200(119,42)=1; mask200(118,44)=1;
mask200(151,20)=1;

%--- Split MAB & SAB ---%
mask200(314:316,243)=0;
mask200(315,242:244)=0;

%--- Gulf of Maine ---%
mask200(408:432,300:337)=1;
mask200(433:437,321:326)=1;
mask200(409:414,338:342)=1;
mask200(431:440,362:374)=1;
mask200(410:411,343)=1;
mask200(411,344)=1;

%--- Bermuda ---%
mask200(278:279,352:353)=1;

return;

%--- Cut off shelf waters south of Yucatan ---%
mask200(85:144,100:135)=0;

%--- Add some shelf waters in Western Gulf of Mexico ---%
mask200(160:167,16)=1; mask200(151,19:20)=1;
mask200(148:150,20)=1; mask200(135,26)=1;
mask200(131:134,27)=1; mask200(131,28)=1;
mask200(    129,29)=1; mask200(118,41)=1;

%--- Cut off Carribbean sea and island chains ---%
mask200(54:155,187:396)=0;

%--- Cut off South American shelf ---%
mask200(1:31,158:end)=0;
mask200(30:55,240:end)=0;

%--- Cut off Cubaan shelf ---%
mask200(158:162,208:217)=0;
mask200(162:163,208:210)=0;
mask200(150:170,145:207)=0;
mask200(173:178,193:199)=0;

%--- Cut off Bahemian shelf ---%
mask200(154:220,202:264)=0;

%--- Add basins from the Gulf of Maine ---%
mask200(408:433,300:337)=1;
mask200(434:437,321:326)=1;
mask200(431:440,362:374)=1;

%--- Cut off Keys and FL Straights ---%
mask200(196:209,193:196)=0; mask200(183    ,175:177)=0;
mask200(183:184,178:186)=0; mask200(193:195,191:194)=0;
mask200(185    ,181:189)=0; mask200(186    ,185:191)=0;
mask200(187    ,187:191)=0; mask200(188    ,189:192)=0;
mask200(189:192,192:194)=0;
mask200(189    ,191    )=0; mask200(210    ,196    )=0;

%--- Split SAB and MAB ---%
mask200(320,240:241)=0; mask200(319,238:242)=0;
mask200(318,239:243)=0; mask200(317,240:244)=0;
mask200(316,241:245)=0; mask200(315,242:246)=0;
mask200(314,243:245)=0; mask200(321,240)    =0;

%==============================================================================%
%             Fine tune deep sea mask to cut off areas of disinterest         =%
%==============================================================================%
%--- Cut off deep pocket off Cuba ---%
mask2500(140:160,235:259)=0;

%--- Cut off Jamaica ---%
mask2500(119:121,255:262)=0;
mask2500(108:112,345:354)=0;
mask2500(45:81,363:390)=0;

%--- Random spot south of Dominican Rep. ---%
mask2500(75,264)=1;

%==============================================================================%
%                            Output new mask values                           =%
%==============================================================================%
%nc_varput(outfile,'mask_shelf',mask200 );
%nc_varput(outfile,'mask_Miss' ,maskMISS);
%nc_varput(outfile,'mask_deep' ,mask2500);

