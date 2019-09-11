function MagB_inv
clc;clear all;clear global;
%Code by 
%Oksum. E. : eroksum@gmail.com / erdincoksum@sdu.edu.tr
%Pham L.T. : luanpt@hus.edu.vn
%%%%%%%%%%%%%%%%%%%% Descriptions %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INVERSION 
     % 1- Import gridded data [2-D grid (*.grd) or columwise data X,Y,Z ]
     %    by Import Data menu 
     % 2- Set field parameters, filter parameters, termination criterion and maximum number
     %    of iteration to related parts of table 
     % 3- Start Iteration..
     % Note: If Rms criterion is set to 0 then the divergence mode is 
     % used for the inversion stopping, otherwise the convergence mode. 
     % The user can set SH and WH parameters also interactively on raps
     % plot (by clicking on Filter pushbutton) by moving the positions of 
     % vertical lines assigning the limits of the roll-off frequencies..
     %%%
     % View Outputs : Basement depth, Inverted anomaly, Anomaly Difference, 
     %                RMS plot 
     %                (interactively selective by a listbox activated by 
     %                mouse click on outputs GUI)
     % View Options : 2D, 3D (only for depth model), 
     %                Cross-Section (depth,observed and calculated anomalies)
     % Export Outputs: Image file (*.png, 500 dpi resolution),
     %                 Data file (maps-> *.grd, graphics-> *.dat)
     %                 Export as project file (*mgb.mat). Project files
     %                 can be loaded any time for re-view of interpretation
     %
%%%%%%%%% FORWARD CALCULATION
     %    1- Activate Forward GUI by Forward-Calc menu at main GUI 
     %    2- Import gridded depth data [2-D grid (*.grd)]
     %       by Load Depth Grid  menu
     %    3- Set field parameters by Settings & Calc menu
     %       (after confirmation of inputs the code performs forward and
     %        plots the calculated magnetic anomaly) 
     %    4- Export as image file/data file by Save menu
     %       [switch between depth/mag maps by Plot menu]

%%%%%%%%%%%%%%% Tasks of some main Functions
%%% impgrd > NEW 2-D GRID and initialize/memorize parameters
         %uses grd2loader,lodgrd6txt,lodgrd7bin
%%% impgrdat> NEW XYZ GRID and initialize/memorize parameters
         %uses lodgrdatfile_Callback
%%% startiter > retrievs inputs from table, calls function maininv,
%%%             memorize outputs
%%% maininv > performs inversion sheme
%%% getfreqs > calculates wavenumbers k, kx,ky
%%% paddData > zero padding of data
%%% filt_LP > lowpass 2-D filter design
%%% raps_data > radially averaged spectrum
%%% calcmag > retriev inputs, initialize parameters, calls forwardmag function
%%%           memorize forward model
%%% forwardmag > performs forward calculation of magnetic anomalies

%%% please see Manual pdf file for the tasks of all functions used
%%% in MagBinv code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating the Structure of main figure window 
createMainwindow 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORT DATA MENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1- NEW 2-D GRID 
function impgrd(~,~)
%imports a new 2-D grid, memorizes to temporary file, updates inputs table content.
closwindows
set(findobj(gcf,'Tag','mainfig'),'Pointer','watch')
[sourcfil,T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy,errm]=grid2loader;
if errm>0;
if errm<5; errormess (errm);end
set(findobj(gcf,'Tag','mainfig'),'Pointer','arrow');
return;
end
%%%%%%%set an initial value for sh and wh 
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
sh_initial=0.1/dx;
wh_initial=0.5*sh_initial;
setto(7)=sh_initial;
setto(8)=wh_initial;
set(findobj(gcf,'Tag','tabl1'),'Data',setto);
%%%%%%%
[~,filnam] = fileparts(sourcfil);
list1w(T,nx,ny,xmin,xmax,ymin,ymax,dx,dy,sourcfil); % print statistics of data to listbox
save('magbinv.mat','T','x','y','nx','ny','xmin','xmax','ymin','ymax','dx','dy',...
    'filnam','-append'); % save/update variables of new input
drawnow; mapper(x,y,T,'nT','Observed Anomaly',1); %% map the input grid
%%% configure gui items 
set(findobj(gcf,'Tag','startbut'),'Enable','on')
set(findobj(gcf,'Tag','rpsbut'),'Enable','on','String','Filter')
set(findobj(gcf,'Tag','mainfig'),'Pointer','arrow')
end

function [sourcfil,matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy,errm]=grid2loader
%Checks format of 2-D grid input, searches for format errors
errm=0;
[filename, pathname] = uigetfile('*.grd', 'Import Golden Software Binary/Text grid (*.grd)');
sourcfil=[pathname filename]; %%source to file
if ischar(sourcfil)
fidc=fopen(sourcfil);
header= fread(fidc,4,'*char' )';
fclose(fidc);
c1=strcmp(header,'DSAA');
c2=strcmp(header,'DSRB');
sumc=sum([c1 c2]);
if sumc>0
switch c1
    case 1
[matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(sourcfil);
    case 0
[matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd7bin(sourcfil);        
end
if dx~=dy; errm=1;return;end
if any(isnan(matrix(:))); errm=3;return;end
else
errm=2;matrix=0;x=0;y=0;nx=0;ny=0;xmin=0;xmax=0;ymin=0;ymax=0;dx=0;dy=0;
end
else
errm=5;matrix=0;x=0;y=0;nx=0;ny=0;xmin=0;xmax=0;ymin=0;ymax=0;dx=0;dy=0;
end
end

function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(sourcfil)
% read Surfer 6 text grid(*.grd) format
surfergrd=fopen(sourcfil,'r'); % Open grid file
dsaa=fgetl(surfergrd);  % Header
% Get the map dimension [NX: East NY: North];
datasize=str2num(fgetl(surfergrd)); nx=datasize(1); ny=datasize(2);
% Map limits: xmin, xmax, ymin ymax
xcoor=str2num(fgetl(surfergrd)); xmin=xcoor(1); xmax=xcoor(2);
ycoor=str2num(fgetl(surfergrd)); ymin=ycoor(1); ymax=ycoor(2);
% check intervals in x and y direction 
dx=(xmax-xmin)/(nx-1);dx=abs(dx);
dy=(ymax-ymin)/(ny-1);dy=abs(dy);
% data limits
anom=str2num(fgetl(surfergrd)); t0min=anom(1); t0max=anom(2);
% data matrix 
[T,numb] = fscanf(surfergrd, '%f', [nx,ny]);
T=T'; % Traspose matrix
fclose(surfergrd);
% map coordinate matrix
[x,y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
end

function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy] = lodgrd7bin(sourcfil)
%read Surfer 7 Binary grid format
fid= fopen(sourcfil);
fread(fid,4,'*char' )';
fread(fid,1,'uint32');fread(fid,1,'uint32');
fread(fid,4,'*char' )';fread(fid,1,'uint32');
ny= fread(fid,1,'uint32'); nx= fread(fid,1,'uint32');
xmin= fread(fid,1,'double'); ymin= fread(fid,1,'double');
dx= fread(fid,1,'double'); dy= fread(fid,1,'double');
fread(fid,1,'double');fread(fid,1,'double');
fread(fid,1,'double');
parm= fread(fid,1,'double');
fread(fid,4,'*char' )';
nn= fread(fid,1,'uint32');
if ny*nx ~= nn/8 ; error('error') ;end
T= nan(nx,ny);
T(1:end) = fread(fid,numel(T),'double');
T=T';
fclose(fid);
T(T==parm) = nan;
xv = xmin + (0:nx-1)*dx;
yv = ymin + (0:ny-1)*dy;
[x,y]=meshgrid(xv,yv);
xmax=xv(end);
ymax=yv(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2- NEW XYZ GRID
function impgrdat(~,~)
%imports new XYZ column data, memorizes to temporary file, updates inputs table content.
[filename, pathname] = uigetfile({'*.dat'; '*.txt'; '*.csv'}, 'Import XYZ grid data');
sourcfil=[pathname filename];
if ischar(sourcfil)
closwindows 
set(findobj(gcf,'Tag','mainfig'),'Pointer','watch')
[T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy,errm]=lodgrdatfile_Callback(sourcfil);
if dx~=dy; errormess (1);return;end
if errm==4; errormess (4);return;end
%%%%%%%set an initial value for sh and wh 
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
sh_initial=0.1/dx;    
wh_initial=0.5*sh_initial;
setto(7)=sh_initial;
setto(8)=wh_initial;
set(findobj(gcf,'Tag','tabl1'),'Data',setto);
%%%%%%%
[~,filnam] = fileparts(sourcfil);
list1w(T,nx,ny,xmin,xmax,ymin,ymax,dx,dy,sourcfil)% print statistics of data to listbox
save('magbinv.mat','T','x','y','nx','ny','xmin','xmax','ymin','ymax','dx','dy',...
    'filnam','-append'); %save/update variables of new input
drawnow
mapper(x,y,T,'nT','Observed Anomaly',1)% map the input data 
%%% configure menu items 
set(findobj(gcf,'Tag','startbut'),'Enable','on')
set(findobj(gcf,'Tag','rpsbut'),'Enable','on','string','Filter')
set(findobj(gcf,'Tag','mainfig'),'Pointer','arrow')
end
end

function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy,errm]=lodgrdatfile_Callback(sourcfil)
%reads XYZ ascci format file (*.dat,*.txt,*.csv), searches for format errors
data=load(sourcfil);
x=data(:,1);y=data(:,2);T=data(:,3); 
unix=unique(x);xmin=unix(1);xmax=unix(end);dx=unix(2)-unix(1);
uniy=unique(y);ymin=uniy(1);ymax=uniy(end);dy=uniy(2)-uniy(1);
nx=numel(unix);
ny=numel(uniy);
dxcount=numel(unique(diff(unix)));
dycount=numel(unique(diff(uniy)));
if sum([dxcount dycount])>2;errm=4;return;end;
errm=0;
if x(2)==x(1);
T=reshape(T,ny,nx);
x=reshape(x,ny,nx);
y=reshape(y,ny,nx);
else
T=reshape(T,nx,ny);
x=reshape(x,nx,ny);
y=reshape(y,nx,ny);
T=T';x=x';y=y';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3- OPEN PROJECT
function openproj(~,~)  
%Re-loads a complete interpretation saved as a project file (matlab binary file). 
[filename, pathname] = uigetfile('*_mgb.mat', 'Import MAGBinv Project file (*_mgb.)');
k=[pathname filename];
if ischar(k)
%%%%%%%%%%%%%%%%%%%%%Get previous interpretation info    
if exist('magbinv.mat','file')==2;delete('magbinv.mat');end
copyfile(k,'magbinv.mat');
load('magbinv.mat');
list1w(T,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k);% retrieve data info
set(findobj(gcf,'Tag','tabl1'),'Data',setto)%%%%retrieve table 1 info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Set Windows and mappings
%%%%%%%%%%%%%%close existing open figures
closwindows
drawnow
mapper(x,y,T,'nT','Observed Anomaly',1)
set(findobj(gcf,'Tag','startbut'),'Enable','on')
set(findobj(gcf,'Tag','rpsbut'),'Enable','on','string','Filter')
%%%%%%%%%%%%%%%%%%%%
%Re-open second window
fig2 = createwindow2;
drawnow
createmenu2(fig2); %%% create menu items into window2
mapper(x,y,-Zcalc,'Km','Calculated Depth',1)%% map the basement
set(findobj(gcf,'Tag','tabl2'),'Data',[max(Zcalc(:)) min(Zcalc(:)) ...
    mean(Zcalc(:)) numel(rmstor)])
end
end

function impordep(~,~) 
% IMPORT DEPTH MODEL GRID FOR FORWARD MAG 
drawnow
set(findobj(gcf,'Tag','Action'),'Pointer','watch')
[~,z,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy,errm]=grid2loader;
if errm>0;%pop up error messages 
if errm<5; errormess (errm);end
set(findobj(gcf,'Tag','Action'),'Pointer','arrow')
return;
end
z=abs(z);
defp = {'90';'0';'90';'0';'2';'10'};
save('magbinvForw.mat','z','x','y','nx','ny','dx','dy','xmin','xmax',...
    'ymin','ymax','defp');
mapper(x,y,-z,'Km','Depth Model',1)
set(findobj(gcf,'Tag','swiftm'),'Enable','off');
set(findobj(gcf,'Tag','savf'),'Enable','off');
set(findobj(gcf,'Tag','Action'),'Pointer','arrow')
end

function [If,Df,Im,Dm,M,z0,SH,WH,criterio,mxiter,alpha,itrmod]=get_parameters
% reads inputs table content, initialize parameters.  
setto=get(findobj(gcf,'Tag','tabl1'),'Data'); 
If=setto(1)*pi/180;
Df=setto(2)*pi/180;
Im=setto(3)*pi/180;
Dm=setto(4)*pi/180;
M=setto(5);
z0=setto(6);
SH=setto(7);
WH=setto(8);
criterio=setto(9);
mxiter=setto(10);
alpha=0;
if criterio==0; itrmod=1;else itrmod=0;end 
save('magbinv.mat','setto','-append')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% START BUTTON AND INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startiter(~,~) 
% retrievs inputs from table content and temp. file, starts inversion procedure, 
% memorizes outputs (basement, calculated magnetic anomaly, 
% anomaly difference, rms variation)
load('magbinv.mat','T','x','y','nx','ny','dx','dy');
[If,Df,Im,Dm,M,z0,SH,WH,criterio,mxiter,alpha,itrmod]=get_parameters;
if SH==0 || WH==0;msgbox({'please set filter configuration';...
        'SH and WH >0'});return;end
closwindows
drawnow
fig2=createwindow2;%%% Open a new figure window
drawnow
createmenu2(fig2); %%% create menu items into window2
freqorder=sort([SH WH]);WH=freqorder(1);SH=freqorder(2);
T0=T; % memorize orginal T 
T=T./(10^9);
%%%%%%%%%%Data padding
[Te,nxe,nye]=paddData(nx,ny,T);
nxm=2*nxe; nym=2*nye;
%%%%%%%%%%%%%%%%%%%%%% get frequencies   
[kx,ky,k]=getfreqs(nxm,nym,dx,dy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Call inversion sheme 
[Zcalc,rmstor]=maininv(Te,kx,ky,k,nx,ny,nxe,nye,Im,Dm,If,Df,M,alpha,...
                       z0,WH,SH,mxiter,itrmod,criterio);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate forward mag and mag difference between actual and calculated..
z=Zcalc-z0; % extract mean depth
[ze,nxe,nye]=paddData(nx,ny,z);
Tcalc=forwardmag(ze,nx,ny,nxe,nye,kx,ky,k,M,Im,If,Dm,Df,alpha,z0);
Tdiff=T0-Tcalc;
%%%%%%%%%%%%%%%% Store outputs  
save('magbinv.mat','Zcalc','Tcalc','Tdiff','rmstor','-append')
set(findobj(gcf,'Tag','tabl2'),'Data',[max(Zcalc(:)) min(Zcalc(:)) ...
    mean(Zcalc(:)) numel(rmstor)])
drawnow
mapper(x,y,-Zcalc,'Km','Calculated Depth',1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVERSION PROCEDURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Zcalc,rmstor]=maininv(T,kx,ky,k,nx0,ny0,nxe,nye,Im,Dm,If,Df,M,alpha,...
                       z0,WH,SH,mxiter,itrmod,criterio) 
% performs the inversion procedure 
cm=10^-7;
[mx,my,mz,fx,fy,fz]=dircosin(Im,Dm,If,Df,alpha);% unit vectors
thetam = mz+(1i).*((mx.*kx+my.*ky)./abs(k)); 
thetaf=fz+(1i).*((fx.*kx+fy.*ky)./abs(k)); 
hs=M*2.*pi.*cm.*thetam.*thetaf.*exp(-k.*z0).*k;
%%%%%%%%%%%%%%%%%%%%Filter design
dimk=size(k);
filt_d=filt_LP(k,dimk(2),dimk(1),WH,SH);% filter cofficients
%%%%%%%%%%%%%%%%%%%% Filter data and get the first approximation             
FT=fft2(T);
Fh=-FT./hs;
Fh=Fh.*filt_d;
Fh(1,1)=0;
h=real(ifft2(Fh));
h_old=h;%%%first approximation of depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Start of iterative procedure  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rms_old=1000;rmstor=1000;
for istep=1:mxiter
     m=istep+1;
     Fh=Fh-((-k).^(m-1)).*fft2(h.^m)./factorial(m);
     Fh=Fh.*filt_d;
     h=real(ifft2(Fh)); %new h in istep 
     dh=h-h_old; 
     dh2=dh.^2;
     rms=sqrt(sum(sum(dh2))./(numel(dh2(:)))); % new rms from new h
 %%%case convergence mode
 if itrmod==0; 
     if rms<criterio;rmstor(istep)=rms;h_old=h;
     instat(h_old,z0,nxe,nye,nx0,ny0,istep);% instant print to table 2 
         break;
     end;
 end
 %%%    
 if rms>rms_old; break;end % new h is not suitable..use h_old of best approximation
 rmstor(istep)=rms; %new rms is stored for next iteration
 h_old=h; % new h is stored
 rms_old=rms;
 instat(h_old,z0,nxe,nye,nx0,ny0,istep)% instant print to table 2  
 end
%%%%%% End of Inversion sheme
%%%%%%%%%%% Initialize output z
z=h_old+z0;
z=z(nye/2+1:nye/2+ny0,nxe/2+1:nxe/2+nx0);
Zcalc=z; %calculated depth 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CALC FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,nx,ny]=paddData(nx,ny,T) 
% PADDING DATA
T(1,nx+floor(nx/2))=0;
T(ny+floor(ny/2),1)=0;
T=rot90(rot90(T));
T(1,nx+2*floor(nx/2))=0;
T(ny+2*floor(ny/2),1)=0;
T=rot90(rot90(T));
if (mod(nx,2)~=0) nx=nx-1; T(:,end)=[]; end
if (mod(ny,2)~=0) ny=ny-1; T(end,:)=[]; end
end

function [kx,ky,k]=getfreqs(nxm,nym,dx,dy) 
% WAVENUMBERS
dkx= 2.*pi./((nxm-1).*dx);
dky= 2.*pi./((nym-1).*dy);
nyqx= (nxm/2)+1;
nyqy= (nym/2)+1; 
[kx,ky]=meshgrid([(0:nxm/2) (nyqx+1:nxm)-(nxm+1)].*dkx,...
    [((0:nym/2)) (nyqy+1:nym)-(nym+1)].*dky);
k= sqrt(bsxfun(@plus,kx.^2,ky.^2));
k(1,1)=0.00000001;
end

function [mx,my,mz,fx,fy,fz]=dircosin(Im,Dm,If,Df,alpha) 
% calculates unit vectors of the magnetization and the ambient field
mx=cos(Im)*sin(Dm-alpha);
my=cos(Im)*cos(Dm-alpha);
mz=sin(Im);
fx=cos(If).*sin(Df-alpha);
fy=cos(If).*cos(Df-alpha);
fz=sin(If);
end

function filt_d=filt_LP(k,nxm,nym,WH,SH)
% 2D filter design
filt_d=zeros(nym,nxm); %% preallocating
ktotal=k./(2*pi); 
for j=1:nym;
  for i=1:nxm;
if ktotal(j,i)<WH
      filt_d(j,i)=1;  
elseif ktotal(j,i)<SH
      filt_d(j,i)=0.5.*(1+cos((((2*pi)*ktotal(j,i))-(2*pi*WH))/(2*(SH-WH))));
else
filt_d(j,i)=0;
end
  end
end
end

function Tcalc=forwardmag(z,nx0,ny0,nxe,nye,kx,ky,k,M,Im,If,Dm,Df,alpha,z0)
% performs forward calculation of magnetic anomalies from gridded depth model
order=10; 
cm=10^-7; 
[mx,my,mz,fx,fy,fz]=dircosin(Im,Dm,If,Df,alpha);
thetam = mz+(1i).*((mx.*kx+my.*ky)./abs(k)); 
thetaf=fz+(1i).*((fx.*kx+fy.*ky)./abs(k)); 
hs=M*2.*pi.*cm.*thetam.*thetaf.*exp(-k.*z0);
SumF=0;
for m=1:order;
SumF=SumF+((-k).^(m))./(factorial(m)).*fft2(z.^m);
end;
FT=hs.*SumF;
T=(ifft2(FT));
T=real(T);
T=T.*10^9;
Tcalc=T(nye/2+1:nye/2+ny0,nxe/2+1:nxe/2+nx0);
end

function calcmag(~,~)
% retrievs inputs, initialize parameters, 
% calls forwardmag function, memorize forward model
prompt = {'Inc. angle of ambient field (Incf)?'; 'Dec. angle of ambient field (Df)?';...
    'Inc. angle of magnetization(Im) ?'; 'Dec. angle of magnetization(Dm)?';...
'Magnetization Contrast (M)?';
'Average interface depth Z0 (unit in Km)?'};
dlg_title = 'parameters';
num_lines = 1;
load('magbinvForw.mat');
answer= inputdlg(prompt,dlg_title,num_lines,defp);
if numel(answer)>0
If=str2num(cell2mat(answer(1)))*pi/180;
Df=str2num(cell2mat(answer(2)))*pi/180;
Im=str2num(cell2mat(answer(3)))*pi/180;
Dm=str2num(cell2mat(answer(4)))*pi/180;
M=str2num(cell2mat(answer(5)));
z0=abs(str2num(cell2mat(answer(6))));
alpha=0;
set(findobj(gcf,'Tag','Action'),'Pointer','watch')
%%%%%%%%%%%%%%%%%%%%%padding depth data
z=z-z0; % extract mean depth
[ze,nxe,nye]=paddData(nx,ny,z);
nxm=2*nxe; nym=2*nye;
%%%%%%%%%%%%%%%%%%%%%% get frequencies   
[kx,ky,k]=getfreqs(nxm,nym,dx,dy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tcalc=forwardmag(ze,nx,ny,nxe,nye,kx,ky,k,M,Im,If,Dm,Df,alpha,z0);
mapper(x,y,Tcalc,'nT','Calculated Anomaly From Depth Model',1)
defp={num2str(If*180/pi);num2str(Df*180/pi);num2str(Im*180/pi);...
    num2str(Dm*180/pi);num2str(M);num2str(z0)};
save('magbinvForw.mat','Tcalc','defp','-append')
set(findobj(gcf,'Tag','swiftm'),'Enable','on');
set(findobj(gcf,'Tag','savf'),'Enable','on');
set(findobj(gcf,'Tag','Action'),'Pointer','arrow')
end
end

function [pwr,wn]= raps_data(T,dx) 
% get radially averaged spectrum vs frequency 
[ny,nx]=size(T); 
nrow=2*floor(ny/2); ncol=2*floor(nx/2); 
maxdim=max([nrow ncol]);
np2=2^nextpow2(maxdim);  
rowdiff=round((np2-nrow)/2); 
coldiff=round((np2-ncol)/2);
T=T(1:nrow,1:ncol);
T=T-mean(T(:)); 
wf=tukeywin(nrow,.05)*tukeywin(ncol,.05)'; %truncation window 5% from edges
dw=T.*wf;
TT=zeros(np2); 
TT(rowdiff+1:rowdiff+nrow,coldiff+1:coldiff+ncol)=dw;  
spec =(abs(fftshift(fft2(TT))));
spec(spec<1.0e-13)=1.0e-13;
spec=log(spec);
% create cartesian coordinate system (unit) and
% shifting the origin to zero frequency.
[xo,yo]=meshgrid((1:np2)-np2/2,(1:np2)-np2/2); 
[~,L]=cart2pol(xo,yo);L=round(L);% calculate distances to the center
halfspec=(np2/2)-1; %considering the half of the spectrum
wn=(1:halfspec)/(np2*dx); % freq axis
for i=1:halfspec
pwr(i)=mean(spec(L==i));%find data of equal distances and get mean  
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PLOTTER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rmsplot(~,~) 
% displays the rms variation obtained during the iteration steps. 
load('magbinv.mat','rmstor');
plot(rmstor,'-ko','MarkerFaceColor','r','MarkerSize',5);
set(gca,'FontSize',10,'FontWeight','bold')
sent1=['RMS(1)=' num2str(rmstor(1))];
sent2=['RMS(end)=' num2str(rmstor(end))];
title ({'RMS-Graph'; [sent1 '    ' sent2]})
xlabel('iteration step');ylabel('RMS');
xlim([1 numel(rmstor)])
      grid on
      pbaspect([1 1 1]);axis square
rotate3d off
set(findobj(gcf,'Tag','uyari'),'string',...
'Click on screen to switch between output maps','ForeGroundColor','b')
set(findobj(gcf,'Tag','plc1'),'enable','off')% set buttons enable or dissable
set(findobj(gcf,'Tag','plc2'),'enable','off')%
set(findobj(gcf,'Tag','plc3'),'enable','off')%
set(findobj(gcf,'Tag','plc4'),'enable','off')%
end

function mapper(x,y,matrix,unitt,tit,index) 
% displays  a color filled contour map of gridded data
if index==1;
contourf(x,y,matrix,18);shading flat;
rotate3d off;
end
set(gca,'FontSize',10,'FontWeight','bold')
h=colorbar('eastoutside');title(h,unitt,'FontWeight','bold');
set(h,'FontSize',10,'FontWeight','bold')
xlabel('X (Km)');ylabel('Y (Km)');title(tit)
axis equal
axis tight
end

function map3v(x,y,z,tit) 
% displays  3D contour map of basement depth.
surf(x,y,z,'FaceColor','interp','EdgeColor','none')
rotate3d on
title(tit)
box on
pbaspect([1 1 .5])
end

function plotcross(PROFX,TC,TCC,ZC) 
% plots cross-section data of observed and calculated anomalies 
% and basement relief
fig3 = figure('MenuBar','none','Name','Cross-Section',...
'NumberTitle','off','Resize','off',...
'Color','w',...
'units','normalized','outerposition',[0 .1 .49 .9],...
'DockControls','off','WindowButtonDownFcn',@posfig1);
ax1=axes('Parent',fig3,'units','normalized','Position',[0.1 0.6 0.8 0.35]);
axis off
ax2=axes('Parent',fig3,'units','normalized','Position',[0.1 0.1 0.8 0.35]);
axis off
uimenu('Parent',fig3,'Label','Export Profile Data','CallBack',@savcross);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ax1)
plot(PROFX,TC,'+r',PROFX,TCC,'-k','linewidth',1);grid on;box on;
xlim([min(PROFX) max(PROFX)]);ylabel('nT');
legend('Observed','Calculated');
%%%%%%%%%%%%%%%
axes(ax2)
dfz=max(ZC)+(max(ZC)/5);
ZC(end+1:end+3)=[dfz dfz ZC(1)]; 
PROFX(end+1:end+3)=[PROFX(end) PROFX(1) PROFX(1)];
plot(PROFX(1:end-3),-ZC(1:end-3),'-k','linewidth',2);
patch(PROFX,-ZC,[.1 .1 .1],'FaceAlpha',.8);
ylabel('Depth (Km)')
xlabel('Distance (Km)')
xlim([min(PROFX) max(PROFX)]);
ylim([min(-ZC) 0]);
box on
grid on
end

function rapsplot(src,~)
% plots the radially averaged spectrum of input data and two vertical lines 
% (in blue) indicating the current roll-off frequencies 
% of the filter design (SH and Wh).
sourL=get(src,'string');
  switch sourL
    case 'Filter' 
load('magbinv.mat','T','dx');
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
[pwr,wn]= raps_data(T,dx); %get radially averaged spectrum
plot(wn,pwr,'k','LineWidth',2);
line([setto(7) setto(7)],[min(pwr) max(pwr)],'Color','b','LineStyle','--','LineWidth',3,...
     'ButtonDownFcn',@setLin,'Tag','L1')
line([setto(8) setto(8)],[min(pwr) max(pwr)],'Color','b','LineStyle','--','LineWidth',3,...
     'ButtonDownFcn',@setLin,'Tag','L2')
pbaspect([1 1 1])
grid on
title ({'Radially Averaged Pow-Spectrum Graph';...
'Set SH and WH frequency band array for lowpass filtering'})
xlabel('k/2*pi');ylabel('Log(P)');
set(src,'string','Observed Anomaly');
    case 'Observed Anomaly'
    load('magbinv.mat','x','y','T')
    mapper(x,y,T,'nT','Observed Mag.',1)
    set(src,'string','Filter');
   end
end

function setLin(src,~)
% enables the user to set the frequency band array (SW-WH) interactively 
% by mouse control
set(src,'Color',[.5 .5 .5])
new_p=ginput(1);
if new_p(1)<0;new_p(1)=0.01;end
set(src,'XData',[new_p(1) new_p(1)],'Color','b')
get_SHWH
end

function get_SHWH
% updates the cells related to SH and WH in inputs table according 
% the interactive selection of these parameters on raps plot
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
Lx1 = get(findobj(gcf,'Tag','L1'), 'XData');
Lx2 = get(findobj(gcf,'Tag','L2'), 'XData');
whsh=unique([Lx1 Lx2]);
setto(8)=whsh(1);
setto(7)=whsh(2);
set(findobj(gcf,'Tag','tabl1'),'Data',setto);
xlabel(['SH = ' num2str(setto(7)) '    WH = ' num2str(setto(8))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% GUI WINDOW MANAGERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function posfig1(~,~) 
% reset position of GUI-1
set(gcf,'outerposition',[0 .1 .49 .9])
end

function fig2=createwindow2  
% create empty secondary figure window
fig2=figure('MenuBar','none','Name','Action',...
'NumberTitle','off','Resize','off',...
'Color','w',...
'units','normalized','outerposition',[.5 .1 .49 .9],...
'DockControls','off');
end

function createMainwindow   
% create Main gui window
fig1 = figure('MenuBar','none','Name','MagB_inv',...
'NumberTitle','off','Resize','off',...
'Tag','mainfig','Color','w',...
'units','normalized','outerposition',[0 .1 .49 .9],...
'DockControls','off','WindowButtonDownFcn',@posfig1);
s1='MagB_inv : Magnetic basement estimation';
s2='2019/9';
uicontrol('Parent',fig1,'Style','text','units','normalized','Position',...
[0.65 0.91 0.33 0.08],'FontWeight','normal','HorizontalAlignment','center',...
'BackGroundColor','w','ForeGroundColor',[.2 .2 .2],...
'string',{s1;s2})
%%%% Create Menu items...
%%%% Import Data
menu1 = uimenu('Parent',fig1,'Label','IMPORT DATA'); 
uimenu('Parent',menu1,'Label','New 2-D Grid (*.grd)','CallBack',@impgrd);
uimenu('Parent',menu1,'Label','New XYZ grid (*.dat, *.txt, *.csv)','CallBack',@impgrdat);
uimenu('Parent',menu1,'Label','OPEN MAGBinv Project','CallBack',@openproj);
uimenu('Parent',fig1,'Label','Forward-Calc','CallBack',@createFWgui);
%%%% Add a listbox for Grid info 
uicontrol('Parent',fig1,'Style','listbox','units','normalized','Position',[.02,.9,.6,.09],...
'String',{'Grid Info:'},'Tag','listb1','BackGroundColor','w')
%%% Create an Axis for plotting of the input map 
axes('Parent',fig1,'units','normalized','Position',...
    [0.2 0.1 0.6 0.6],'Tag','axsig');
axis off
%%% Create a table for parameter inputs
tbl1=uitable('Parent',fig1,'units','normalized','Position',[.02,.79,.97,.11],...
'ColumnName',{'IncF','DecF','IncM','DecM','M','AvgZo','SH','WH','RMS crit.','Maxiter'},...
'Tag','tabl1','Rowname',{'Parameter'},...
'ColumnWidth',{40,40,40,40,30,80,70,70,70,50},'ColumnEditable',true,...
'CellEditCallback',@tbl_selectshwh);
uicontrol('Parent',fig1,'Style','pushbutton','units','normalized',...
'Position',[.7,.795,.28,.04],...
'String','Start Iteration','Tag','startbut','ForeGroundColor','b',...
'FontWeight','bold','CallBack',@startiter,'enable','off')
uicontrol('Parent',fig1,'Style','pushbutton','units','normalized',...
'Position',[.03,.795,.2,.04],...
'String','Filter','Tag','rpsbut','ForeGroundColor','k',...
'FontWeight','bold','CallBack',@rapsplot,'enable','off')
text(.30,.0,'MagB_ inv','FontWeight','bold','FontSize',20,'Color',[.75 .75 .75])
axis off
%%%%%%%%%%%%%%%%% Load Settings from magbinv.mat (if exists)
%%%% otherwise print table with default values and create a magbinv.mat file
%%%% to store in use for next run of code....
if exist('magbinv.mat','file')==2;
load('magbinv.mat','setto');
    if numel(setto)==10;
    set(tbl1,'Data',setto);
    else setto=[90 0 90 0 2 10 0.06 0.03 10^(-4) 20];
        save('magbinv.mat','setto');
        set(tbl1,'Data',setto);
    end
else
setto=[90 0 90 0 2 10 0.06 0.03 10^(-4) 20]; %default settings 
%%%=[IncF,DecF,IncM,DecM,M,AvgZo,SH,WH,DVRG,CNVG,Maxiter]
save('magbinv.mat','setto');
set(tbl1,'Data',setto);
end
end

function tbl_selectshwh(~,event)
% updates the positions of the vertical lines related to the roll-off 
% frequencies SH and Wh in raps plot according the manual input in inputs table
if event.Indices(2)==7 || event.Indices(2)==8
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
set(findobj(gcf,'Tag','L1'), 'XData',[setto(7) setto(7)]);
set(findobj(gcf,'Tag','L2'), 'XData',[setto(8) setto(8)]); 
ordershwh=unique([setto(7) setto(8)]);
setto(7:8)=fliplr(ordershwh);
set(findobj(gcf,'Tag','tabl1'),'Data',setto);
xlabel(['SH = ' num2str(setto(7)) '    WH = ' num2str(setto(8))])
end    
end


function createmenu2(fig2) 
% CREATE GUI-2 OUTPUTS WINDOW
set(fig2,'WindowButtonDownFcn', @selectoutmap)
axes('Parent',fig2,'units','normalized','Position',[0.15 0.02 0.75 0.75]);
axis off
menu1 = uimenu('Parent',fig2,'Label','SAVE/Export'); 
uimenu('Parent',menu1,'Label','Save as MAGBinv Project (*mat)','CallBack',@savbut);
menu2=uimenu('Parent',menu1,'Label','Export');
uimenu('Parent',menu2,'Label','RMS List (*.dat)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Calc-Z grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Calc-T grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Diff. (T-CalcT) grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Export All','CallBack',@savbut);
uimenu('Parent',fig2,'Label','Save Plots as Images','CallBack',@savbut);
uitable('Parent',fig2,'units','normalized','Position',[.02,.87,.97,.07],...
'ColumnName',{'Max Depth','Min Depth','Mean Depth','Iteration Step'},...
'Rowname',{'STAT'},...
'ColumnWidth',{100,100,100,100},...
'Enable','off','Data',[0 0 0 0],'Tag','tabl2');
uicontrol('Parent',fig2,'Style','text','units','normalized',...
'Position',[.02,.81,.4,.04],'String','Click on screen to switch between output maps',...
'BackGroundColor','w','ForeGroundColor','b','Tag','uyari',...
'HorizontalAlignment','left')
uicontrol('Parent',fig2,'Style','pushbutton','units','normalized',...
'Position',[.85,.83,.13,.03],'String','Color-Map','Tag','plc1',...
'ForeGroundColor','k','CallBack',@plotcontrol)
uicontrol('Parent',fig2,'Style','pushbutton','units','normalized',...
'Position',[.71,.83,.13,.03],'String','View-3','Tag','plc2',...
'ForeGroundColor','k','CallBack',@plotcontrol)
uicontrol('Parent',fig2,'Style','pushbutton','units','normalized',...
'Position',[.55,.83,.15,.03],'String','View-2 (default)','Tag','plc3',...
'ForeGroundColor','k','CallBack',@plotcontrol)
uicontrol('Parent',fig2,'Style','pushbutton','units','normalized',...
'Position',[.8,.88,.18,.05],'String','Cross-Section','Tag','plc4',...
'ForeGroundColor','k','CallBack',@plotcontrol)
end

function createFWgui (~,~)
% creates menu items related to forward-calc procedure
closwindows
fig2 = createwindow2;
uimenu('Parent',fig2,'Label','Load Depth Grid','CallBack',@impordep);
uimenu('Parent',fig2,'Label','Settings & Calc','CallBack',@calcmag);
sm=uimenu('Parent',fig2,'Label','Plot','Tag','swiftm','Enable','off');
uimenu('Parent',sm,'Label','Depth Model','CallBack',@switcmap);
uimenu('Parent',sm,'Label','Calc-Mag','CallBack',@switcmap);
savf=uimenu('Parent',fig2,'Label','Save','Tag','savf','Enable','off');
uimenu('Parent',savf,'Label','Get Screen as Image File(*.png)','CallBack',@saveforward);
uimenu('Parent',savf,'Label','Save Mag-Data File (*.grd)','CallBack',@saveforward);
axes('Parent',fig2,'units','normalized','Position',[0.15 0.1 0.75 0.75]);
axis off
end

function expofig (~,~)
% creates a GUI for taking screen shots of desired output plots.
closwindows
fig2 = createwindow2;
uimenu('Parent',fig2,'Label','PRINT','CallBack',@ploprint);
uimenu('Parent',fig2,'Label','RETURN','CallBack',@return2menu2);
set(fig2,'WindowButtonDownFcn', @selectoutmapimag)
axes('Parent',fig2,'units','normalized','Position',[0.15 0.1 0.75 0.75]);
text(0.2,0.6,{'1. Click on screen to select an output map/graphic';...
              '2. Use Print menu to save as image file...'})
axis off
end

function selectoutmap(~,~) 
% displays a dialog box for selecting the desired output to plot.  
set(gcf,'outerposition',[.5 .1 .49 .9])
str = {'Basement Depth';'Inverted Anomaly';'Anomaly Difference';
    'RMS Plot'};
      [s,v] = listdlg('PromptString','Select Output Map:',...
                      'SelectionMode','single',...
                      'ListString',str,'ListSize',[160 80]);

load('magbinv.mat','x','y','Tcalc','Zcalc','Tdiff');
 if v==1;    
  switch s
      case 1
      mapper(x,y,-Zcalc,'Km','Calculated Depth',1)
      set(findobj(gcf,'Tag','plc1'),'enable','on')
      set(findobj(gcf,'Tag','plc2'),'enable','on')
      set(findobj(gcf,'Tag','plc3'),'enable','on')
      set(findobj(gcf,'Tag','plc4'),'enable','on')
      case 2
      mapper(x,y,Tcalc,'nT','Calculated Mag.',1)
      set(findobj(gcf,'Tag','plc1'),'enable','on')
      set(findobj(gcf,'Tag','plc2'),'enable','off')
      set(findobj(gcf,'Tag','plc3'),'enable','off')
      set(findobj(gcf,'Tag','plc4'),'enable','on')
      case 3
      mapper(x,y,Tdiff,'nT','Anomaly Difference',1)
      set(findobj(gcf,'Tag','plc1'),'enable','on')
      set(findobj(gcf,'Tag','plc2'),'enable','off')
      set(findobj(gcf,'Tag','plc3'),'enable','off')
      set(findobj(gcf,'Tag','plc4'),'enable','on')
      case 4
      rmsplot
      
  end
 end
                  
end

function selectoutmapimag(~,~) 
% displays a dialog box for selecting the desired output to save as an image (screen shot).
set(gcf,'outerposition',[.5 .1 .49 .9])
str = {'Basement Depth';'Inverted Anomaly';'Anomaly Difference';
    'RMS Plot';'Observed Anomaly'};
      [s,v] = listdlg('PromptString','Select Output Map:',...
                      'SelectionMode','single',...
                      'ListString',str,'ListSize',[200 120]);
load('magbinv.mat','x','y','T','Tcalc','Zcalc','Tdiff','rmstor');
 if v==1;    
  switch s
      case 1
      mapper(x,y,-Zcalc,'Km','Calculated Depth',1)
      case 2
      mapper(x,y,Tcalc,'nT','Calculated Mag.',1)
      case 3
      mapper(x,y,Tdiff,'nT','Anomaly Difference',1)
      case 4
      plot(rmstor,'-ko','MarkerFaceColor','r','MarkerSize',5);
      set(gca,'FontSize',10,'FontWeight','bold')
      sent1=['RMS(1)=' num2str(rmstor(1))];
      sent2=['RMS(end)=' num2str(rmstor(end))];
      title ({'RMS-Graph'; [sent1 '    ' sent2]})
      xlabel('Iteration number');ylabel('RMS');
      xlim([1 numel(rmstor)])
      grid on
      pbaspect([1 1 1]);axis square
      case 5
      mapper(x,y,T,'nT','Observed Mag.',1)
  end
 end
end

function return2menu2(~,~) 
% closes GUI created by expofig and returns to GUI created by createmenu2
clc;clear all;clear global;
load('magbinv.mat','Zcalc','rmstor');
%%%%%%%%%%%%%%close existing open figures
closwindows
%%%%%%%%%%%%%%%%%%%%
%return to Second window
fig2 = createwindow2;
drawnow
createmenu2(fig2)
set(findobj(gcf,'Tag','tabl2'),'Data',[max(Zcalc(:)) min(Zcalc(:)) ...
mean(Zcalc(:)) numel(rmstor)])
end

function closwindows 
% closes all windows opened (except the main GUI)
delete(findobj('Type','figure','name','Action'));
delete(findobj('Type','figure','name','Cross-Section'));
delete(findobj('Type','figure','name','err-m'));
end

function errormess (e) 
% displays an error message window explaining the kind of error.  
set(findobj(gcf,'Tag','mainfig'),'Pointer','arrow')    
figure('MenuBar','none','Name','err-m',...
'NumberTitle','off','Resize','off',...
'Color','k',...
'units','normalized','outerposition',[0.05 .8 .4 .1],...
'DockControls','off');
switch e
    case 1
text(0.2,.5,'equal spaced grid dx=dy is required..','Color','w')
    case 2
text(0,.5,'File format not supported..Load Surfer 6 text grid / Surfer 7 Binary Grid','Color','w')
    case 3
text(0.3,.5,'Blanked grid not supported..','Color','w')
    case 4
text(0.1,.5,'XYZ data input is not equal spaced..(check intervals of X or Y colums)','Color','w')
end
axis off
end

function plotcontrol(src,~) 
% checks and directs to the plot style that user has selected  
% (2D, 3D, cross-section). 
srclabel=get(src,'String');
switch srclabel
    case 'Color-Map'
        colormapeditor
    case 'View-3'
    ax = gca;
    ch=get(ax,'Children');
    %data at current ax
    x=get(ch,'XData');y=get(ch,'YData');z=get(ch,'ZData');
    titax=get(ax,'title');
    tit=get(titax,'String');% title of current ax 
    map3v(x,y,z,tit)
    set(findobj(gcf,'Tag','uyari'),'string','Rotate ON','ForeGroundColor','r')
    set(findobj(gcf,'Tag','plc4'),'enable','off')
    case  'View-2 (default)'
    ax = gca;
    ch=get(ax,'Children');
    %data at current ax
    x=get(ch,'XData');y=get(ch,'YData');z=get(ch,'ZData');
    titax=get(ax,'title');
    tit=get(titax,'String');% title of current ax 
    mapper(x,y,z,'Km',tit,1)
    set(findobj(gcf,'Tag','uyari'),'string',...
    'Click on screen to switch between output maps','ForeGroundColor','b')
    set(findobj(gcf,'Tag','plc4'),'enable','on')
    case 'Cross-Section'
    getcross
end
end

function switcmap (src,~)
% toggle between contour maps of imported depth model and the calculated
% magnetic anomaly 
srclabel=get(src,'Label');
drawnow;set(findobj(gcf,'Tag','Action'),'Pointer','watch');
switch srclabel
    case 'Depth Model'
    load('magbinvForw.mat','x','y','z');mapper(x,y,-z,'Km','Depth Model',1)
    case 'Calc-Mag'
    load('magbinvForw.mat','x','y','Tcalc');mapper(x,y,Tcalc,'nT','Mag. Model',1)    
end
drawnow;set(findobj(gcf,'Tag','Action'),'Pointer','arrow');
end


function getcross 
% creates a new GUI for cross-sectional view, interpolates data 
% (observed anomaly, calculated anomaly and basement) between interactively 
% selected two coordinates   
drawnow
delete(findobj(gcf,'Tag','crossLin'))
set(findobj(gcf,'Tag','plc2'),'enable','off')%% view2 and view3 dissabled
set(findobj(gcf,'Tag','plc3'),'enable','off')
xL=xlabel('Start Line: left mouse click    End Line: right mouse click');
set(xL,'Color','w','FontWeight','Bold','BackGroundColor','k')
figure3 = findobj('Type','figure','name','Cross-Section');
delete (figure3);
load('magbinv.mat','x','y','T','Tcalc','Zcalc');
[xxx,yyy]=getline;
if length(xxx)>1
set(xL,'Color','k','FontWeight','bold','BackGroundColor','w',...
    'String','X (Km)')
xp1=xxx(1);yp1=yyy(1);xp2=xxx(2);yp2=yyy(2);
xc=linspace(xp1,xp2,100);
yc=linspace(yp1,yp2,100);
deltx=abs(xp1-xp2);delty=abs(yp1-yp2);PL=hypot(deltx,delty);
PROFX=linspace(0,PL,100);% profile X
TC=interp2(x,y,T,xc,yc); % profile value Obs Mag.
TCC=interp2(x,y,Tcalc,xc,yc);%profile calc Mag.
ZC=interp2(x,y,Zcalc,xc,yc);% profile basement
hold on;line(xc,yc,'Color','k','linewidth',2,'Tag','crossLin');hold off;
plotcross(PROFX,TC,TCC,ZC)
save('magbinv.mat','TC','TCC','ZC','xc','yc','PROFX','-append')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SAVING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savbut(src,~) 
% enables the user to select the data to be exported.  
srclabel=get(src,'Label');
switch srclabel
    case 'Save as MAGBinv Project (*mat)' 
    savmat
    case 'RMS List (*.dat)'
    savout(1)    
    case 'Calc-Z grid (*.grd)'
    savout(2)    
    case 'Calc-T grid (*.grd)'
    savout(3)    
    case 'Diff. (T-CalcT) grid (*.grd)'
    savout(4)
    case 'Export All'
    savout(5)
    case 'Save Plots as Images'
    expofig    
end
end

function savout(typ) 
% exports the desired output as data file. Maps are stored as 2D data 
% in *.grd format where as the RMS variation is stored in ascii format (*.dat).
load('magbinv.mat');
[filename, pathname] = uiputfile([filnam '.mgb'],...
'Set filename of data: an extension will be determined automatically according the data type');
kk=[pathname filename];
if ischar(kk)
[~,name,~] = fileparts(kk);
kk=[pathname name];  
switch typ
    case 1
        rmstor=rmstor';n=1:numel(rmstor);n=n';matrix=[n rmstor];
        save([kk '-rms.dat'],'matrix','-ascii')
    case 2
        grdout(-Zcalc,xmin,xmax,ymin,ymax,[kk '-Zcalc.grd'])
    case 3
        grdout(Tcalc,xmin,xmax,ymin,ymax,[kk '-Tcalc.grd'])
    case 4
        grdout(Tdiff,xmin,xmax,ymin,ymax,[kk '-Tdiff.grd'])
    case 5
        rmstor=rmstor';n=1:numel(rmstor);n=n';matrix=[n rmstor];
        save([kk '-rms.dat'],'matrix','-ascii')
        grdout(-Zcalc,xmin,xmax,ymin,ymax,[kk '-Zcalc.grd'])
        grdout(Tcalc,xmin,xmax,ymin,ymax,[kk '-Tcalc.grd'])
        grdout(Tdiff,xmin,xmax,ymin,ymax,[kk '-Tdiff.grd'])
end
end
end

function savmat
% exports all inputs/outputs as a single matlab binary file (*.mat) 
load('magbinv.mat');
[filename, pathname] = uiputfile([filnam '.mat'], ' Set a Project Name');
[~,filename,~]=fileparts(filename);
filename=[filename '_mgb.mat'];
kk=[pathname filename];
if ischar(kk)
save(kk,'x','y','T','Tcalc','Tdiff','Zcalc','rmstor',...
    'nx','ny','dx','dy','xmin','xmax','ymin','ymax','setto','filnam')
msgbox({'Project is Stored as *mgd.mat file...';...
    ' can be re-Open by Import Data > OPEN MAGBinv Project)'})
end
end

function savcross(~,~) 
% exports cross-section data to an ascii file (*.dat)  
[filename, pathname] = uiputfile(['profile' '.dat'], 'Set a Profile Name');
kk=[pathname filename];
if ischar(kk)
load('magbinv.mat','TC','TCC','ZC','xc','yc','PROFX');    
fidd=fopen(kk,'wt');
matrix=[xc' yc' PROFX' TC' TCC' -ZC'];
fprintf(fidd,'%s %s %s %s %s %s','Long','Lat','PROFX','Mag','CalcMag','CalcZ');
fprintf(fidd,'\n');
for jj=1:numel(xc); % Write matrix
fprintf(fidd,'%f %f %f %f %f %f',matrix(jj,:));
fprintf(fidd,'\n');
end
fclose(fidd);
title('PROFILE DATA STORED','BackGroundColor','g');
end
end

function ploprint(src,event) 
% exports plot on screen as *png image format with 500 dpi in resolution
clc;clear all; clear global;
[filename, pathname] = uiputfile('*.png','Set a filename');
kk=[pathname filename];
if ischar(kk)
print(kk, '-dpng', '-r500');
drawnow
plot(.5,.5,'k+')
text(0.1,0.6,{'1. Image saved....Click on screen for new selection.';...
              '2. Use Print menu to save as image file...'})
axis off
end
end

function saveforward(src,~) 
% exports the calculated magnetic anomaly grid due to the depth model. 
srclabel=get(src,'Label');
switch srclabel
    case 'Get Screen as Image File(*.png)'
        ext='*.png';
    case 'Save Mag-Data File (*.grd)'
        ext='*.grd';
end
[filename, pathname] = uiputfile(ext,'Set a filename');
fsourc=[pathname filename];
if ischar(fsourc)
if strcmp(srclabel,'Get Screen as Image File(*.png)');
   print(fsourc, '-dpng', '-r500');
   title('Saved as image file...')
else
    set(findobj(gcf,'Tag','Action'),'Pointer','watch')
    load('magbinvForw.mat','Tcalc','x','y','xmin','xmax','ymin','ymax');
    mapper(x,y,Tcalc,'nT','Mag-Model is saved as grid file...',1)
    [~,fnam,~]=fileparts(fsourc);
    grdout(Tcalc,xmin,xmax,ymin,ymax,[pathname fnam '-Tcalc.grd']);
    set(findobj(gcf,'Tag','Action'),'Pointer','arrow')
end
end    
end

function grdout(matrix,xmin,xmax,ymin,ymax,namefile)
% saves output as grid format (*.grd)
% Get grid dimensions
aux=size(matrix);
nx=aux(2);ny=aux(1);
grdfile=fopen(namefile,'w');                % Open file
fprintf(grdfile,'%c','DSAA');               % Header code
fprintf(grdfile,'\n %i %i',[nx ny]);        % Grid size
fprintf(grdfile,'\n %f %f',[xmin xmax]);    % X limits
fprintf(grdfile,'\n %f %f',[ymin ymax]);    % Y limits
fprintf(grdfile,'\n %f %f',[min(min(matrix)) max(max(matrix))]); % Z limits
fprintf(grdfile,'\n');
for jj=1:ny                                 % Write matrix
    for ii=1:nx
        fprintf(grdfile,'%g %c',matrix(jj,ii),' ');
    end
    fprintf(grdfile,'\n');
end
fclose(grdfile);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Table print tasks %%%%%%%%%%%%%%%%%%%%%%%%%
function list1w(T,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k) 
% wites the mesh information of loaded data, as well as root of file 
% to the listbox located at the left upper part of the main GUI window. 
set(findobj(gcf,'Tag','listb1'),'value',1);
s0='Grid Info:';
s1=['Source : ' k];
s2=['NX : ' num2str(nx) '   NY : ' num2str(ny)];
s3=['Grid Spacing     dx: ' num2str(dx) '   dy: ' num2str(dy)];
s4=['Xmin - Xmax:   ' num2str(xmin) '  /   ' num2str(xmax)];
s5=['Ymin - Ymax:   ' num2str(ymin) '  /   ' num2str(ymax)];
s6=['Zmin - Zmax:   ' num2str(min(min(T))) '  /  ' num2str(max(max(T)))];
str={s0;s1;s2;s3;s4;s5;s6};
set(findobj(gcf,'Tag','listb1'),'string',str);
end

function instat(h_old,z0,nxe,nye,nx0,ny0,istep) 
% instant display of the maximum/minimum/mean depth and the current 
% step of iteration during the on-going iterative procedure 
% (in the table in outputs GUI)  
z=h_old+z0;
z=z(nye/2+1:nye/2+ny0,nxe/2+1:nxe/2+nx0);
drawnow
set(findobj(gcf,'Tag','tabl2'),'Data',[max(z(:)) min(z(:)) mean(z(:)) istep])
drawnow
end
