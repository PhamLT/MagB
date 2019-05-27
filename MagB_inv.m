function MagB_inv 
clc;clear all;clear global;
%Code by 
%Oksum. E. : eroksum@gmail.com / erdincoksum@sdu.edu.tr
%Pham L.T. : luanpt@hus.edu.vn

%%%%%%%%%%%%%%%%%%%%%%%%%%% Structure of main figure window
%%%MainFig
fig1 = figure('MenuBar','none','Name','MagB_inv',...
'NumberTitle','off','Resize','off',...
'Tag','mainfig','Color','w',...
'units','normalized','outerposition',[0 .1 .49 .9],...
'DockControls','off','WindowButtonDownFcn',@posfig1);
s1='MagB_inv : Magnetic basement estimation';
s3='2019';
uicontrol('Parent',fig1,'Style','text','units','normalized','Position',...
[0.65 0.91 0.33 0.08],'FontWeight','normal','HorizontalAlignment','center',...
'BackGroundColor','w','ForeGroundColor',[.2 .2 .2],...
'string',{s1;s3})
%%%% Menu items
%%%% Import Data
menu1 = uimenu('Parent',fig1,'Label','IMPORT DATA'); 
uimenu('Parent',menu1,'Label','New 2-D Grid (*.grd)','CallBack',@impgrd);
uimenu('Parent',menu1,'Label','New XYZ grid (*.dat)','CallBack',@impgrdat);
uimenu('Parent',menu1,'Label','OPEN MAGBinv Project','CallBack',@openproj);
uimenu('Parent',fig1,'Label','Forward-Calc','CallBack',@createAnom);
%%%%%%%%%%%%%%%%%
uicontrol('Parent',fig1,'Style','listbox','units','normalized','Position',[.02,.9,.6,.09],...
'String',{'Grid Info:'},'Tag','listb1','BackGroundColor','w')
%%%%%% Create an Axis 
axes('Parent',fig1,'units','normalized','Position',...
    [0.15 0.02 0.75 0.75],'Tag','axsig');
axis off
tbl1=uitable('Parent',fig1,'units','normalized','Position',[.02,.79,.97,.11],...
'ColumnName',{'IncF','DecF','IncM','DecM','M','AvgZo','DVRG','CNVG','Maxiter'},...
'Tag','tabl1','Rowname',{'Parameter'},...
'ColumnWidth',{55,55,55,55,55,85,55,55,55},...
'Enable','off');
c = uicontextmenu;
uicontrol(fig1,'UIContextMenu',c,'string','Initial-Settings (right click for update)',...
'units','normalized','ForeGroundColor','b',...
'Position',[.022,.795,.45,.04],'Tag','setto');
c1=uimenu('Parent',c,'Label','Set Initial Parameters','Tag','c1',...
    'Separator','on','Callback',@setact1,'Separator','on');
c2=uimenu('Parent',c,'Label','Set Iteration-Stop','Checked','on','Separator','on');
uimenu('Parent',c2,'Label','RMS Divergence (RMS(i+1)>RMS(i))','Checked','on','Tag','c21',...
    'Callback',@setact1,'Separator','on');
uimenu('Parent',c2,'Label','RMS Convergence (RMS<threshold)','Tag','c22',...
    'Callback',@setact1,'Separator','on');
uimenu('Parent',c2,'Label','MaxIter','Checked','on','Tag','c23',...
    'Callback',@setact1,'Separator','on');
c3=uimenu('Parent',c,'Label','Filter Observed Grid','Tag','c3',...
'Separator','on','Enable','off','ForeGroundColor','b','Callback',@setact1);
uicontrol('Parent',fig1,'Style','pushbutton','units','normalized',...
'Position',[.7,.795,.28,.04],...
'String','Start Iteration','Tag','startbut','ForeGroundColor','b',...
'FontWeight','bold','CallBack',@startiter,'enable','off')
%%%%%%%%%%%%%%%%% Default Settings
if exist('magbinvset.mat','file')==2;
a=load('magbinvset.mat');
set(tbl1,'Data',a.setto);
if a.setto(7)==0;set(findobj(gcf,'Tag','c21'),'Checked','off');
set(findobj(gcf,'Tag','c22'),'Checked','on');end
else
setto=[90 0 90 0 2 15 1 0 20]; %default settings
save('magbinvset.mat','setto');
set(tbl1,'Data',setto);
%fileattrib('magbinvset.mat','+h');%hidden temporary import info file
end
%%%%%%%%%%%%%%%%%
text(.30,.0,'MagB_ inv','FontWeight','bold','FontSize',20,'Color',[.75 .75 .75])
axis off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function openproj(src,event)   %%%Retriev a Magbinv project
clc;clear all;clear global;
[filename, pathname] = uigetfile('*_mgb.mat', 'Import MAGBinv Project file (*_mgb.)');
k=[pathname filename];
if ischar(k)
set(findobj(gcf,'Tag','mainfig'),'Pointer','watch')    
%%%%%%%%%%%%%%%%%%%%%Get previous interpretation info    
if exist('magbinv.mat','file')==2;delete('magbinv.mat');end
if exist('magbinvset.mat','file')==2;delete('magbinvset.mat');end
copyfile(k,'magbinv.mat');
load('magbinv.mat');
list1w(T,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k);% retrieve data info
set(findobj(gcf,'Tag','tabl1'),'Data',setto)%%%%retrieve table 1 info
save('magbinvset.mat','setto')
set(findobj(gcf,'Tag','mainfig'),'Pointer','arrow')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Set Windows and mappings
%%%%%%%%%%%%%%close existing open figures
figure2 = findobj('Type','figure','name','Action');
delete (figure2);
figure3 = findobj('Type','figure','name','Cross-Section');
delete (figure3);
drawnow
mapper(x,y,T,'nT','Observed Anomaly',1)
set(findobj(gcf,'Tag','c3'),'Enable','on')
set(findobj(gcf,'Tag','startbut'),'Enable','on')
%%%%%%%%%%%%%%%%%%%%
if SH==0 || WH==0;msgbox({'please set filter configuration';...
'SH and WH >0'}); setfilters;return;end
%Secondfig
fig2 = figure('MenuBar','none','Name','Action',...
 'NumberTitle','off','Resize','off',...
 'Color','w',...
 'units','normalized','outerposition',[.5 .1 .49 .9],...
 'DockControls','off');
 set(fig2,'WindowButtonDownFcn', @selectoutmap)
 axes('Parent',fig2,'units','normalized','Position',[0.15 0.02 0.75 0.75]);
 axis off
psbut=uicontrol(fig2,'Style','pushbutton','units','normalized',...
'BackGroundColor',[.8 .8 .8],'HorizontalAlignment','right',...
'Position',[0 .95 1 0.03],'String','RMS Plot','Callback',@rmsplot);
mapper(x,y,-Zcalc,'Km','Calculated Depth',1)
menu1 = uimenu('Parent',fig2,'Label','SAVE/Export'); 
uimenu('Parent',menu1,'Label','Save as MAGBinv Project (*mat)','CallBack',@savbut);
menu2=uimenu('Parent',menu1,'Label','Export');
uimenu('Parent',menu2,'Label','RMS List (*.dat)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Calc-Z grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Calc-T grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Diff. (T-CalcT) grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Export All','CallBack',@savbut);
uimenu('Parent',fig2,'Label','Save Plots as Images','CallBack',@savbut);
maxd=max(max(Zcalc));
 mind=min(min(Zcalc));
 meand=mean(mean(Zcalc));
 rs=numel(rmstor);
 uitable('Parent',fig2,'units','normalized','Position',[.02,.87,.97,.07],...
 'ColumnName',{'Max Depth','Min Depth','Mean Depth','Iteration Step'},...
 'Rowname',{'STAT'},...
 'ColumnWidth',{100,100,100,100},...
 'Enable','off','Data',[maxd mind meand rs]);
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
end

function posfig1(src,event)  % reset position of GUI
set(gcf,'outerposition',[0 .1 .49 .9])
end

function startiter(src,event) % start iterative procedure and create output GUI
load('magbinv.mat','SH','WH');
if SH==0 || WH==0;msgbox({'please set filter configuration';...
        'SH and WH >0'}); setfilters;return;end
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
%%%Secondfig
figure2 = findobj('Type','figure','name','Action');
delete (figure2);
fig2 = figure('MenuBar','none','Name','Action',...
'NumberTitle','off','Resize','off',...
'Color','w',...
'units','normalized','outerposition',[.5 .1 .49 .9],...
'DockControls','off');
set(fig2,'WindowButtonDownFcn', @selectoutmap)
axes('Parent',fig2,'units','normalized','Position',[0.15 0.02 0.75 0.75]);
axis off
pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALL MAIN Inversion
maininv(setto,fig2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

menu1 = uimenu('Parent',fig2,'Label','SAVE/Export'); 
uimenu('Parent',menu1,'Label','Save as MAGBinv Project (*mat)','CallBack',@savbut);
menu2=uimenu('Parent',menu1,'Label','Export');
uimenu('Parent',menu2,'Label','RMS List (*.dat)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Calc-Z grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Calc-T grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Diff. (T-CalcT) grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Export All','CallBack',@savbut);
menu2 = uimenu('Parent',fig2,'Label','Save Plots as Images','CallBack',@savbut);
load('magbinv.mat','Zcalc','Tdiff','rmstor');
maxd=max(max(Zcalc));
mind=min(min(Zcalc));
meand=mean(mean(Zcalc));
rs=numel(rmstor);
tbl2=uitable('Parent',fig2,'units','normalized','Position',[.02,.87,.97,.07],...
'ColumnName',{'Max Depth','Min Depth','Mean Depth','Iteration Step'},...
'Rowname',{'STAT'},...
'ColumnWidth',{100,100,100,100},...
'Enable','off','Data',[maxd mind meand rs]);
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

function plotcontrol(src,event) % manage output plots
nedir=get(src,'String');
switch nedir
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

function getcross % extract profile data
xL=xlabel('Start Line: left mouse click    End Line: right mouse click');
set(xL,'Color','w','FontWeight','Bold','BackGroundColor','k')
figure3 = findobj('Type','figure','name','Cross-Section');
delete (figure3);
load('magbinv.mat','x','y','T','Tcalc','Zcalc');
[xxx,yyy]=getline;
if length(xxx)>1
set(xL,'Color','k','FontWeight','normal','BackGroundColor','w',...
    'String','Easting (Km)')
xp1=xxx(1);yp1=yyy(1);xp2=xxx(2);yp2=yyy(2);
xc=linspace(xp1,xp2,100);
yc=linspace(yp1,yp2,100);
deltx=abs(xp1-xp2);delty=abs(yp1-yp2);PL=hypot(deltx,delty);
PROFX=linspace(0,PL,100);% profile X
TC=interp2(x,y,T,xc,yc); % profile value Obs Mag.
TCC=interp2(x,y,Tcalc,xc,yc);%profile calc Mag.
ZC=interp2(x,y,Zcalc,xc,yc);% profile basement
plotcross(PROFX,TC,TCC,ZC)
save('magbinv.mat','TC','TCC','ZC','xc','yc','PROFX','-append')
end
end

function plotcross(PROFX,TC,TCC,ZC) % plot of extracted profile
%%%thirdfig
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
xlim([min(PROFX) max(PROFX)]);ylabel('nT');xlabel('Distance (Km)')
legend('Observed','Calculated');
%%%%%%%%%%%%%%%
axes(ax2)
dfz=max(ZC)+(max(ZC)/5);
ZC(end+1:end+3)=[dfz dfz ZC(1)]; 
PROFX(end+1:end+3)=[PROFX(end) PROFX(1) PROFX(1)];
plot(PROFX(1:end-3),-ZC(1:end-3),'-k','linewidth',2);
patch(PROFX,-ZC,[.1 .1 .1],'FaceAlpha',.8);
ylabel('Depth (Km)')
xlim([min(PROFX) max(PROFX)]);
ylim([min(-ZC) 0]);
box on
grid on
end

function savcross(src,event) % storing profile data
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

function map3v(x,y,z,tit) % 3D view option
surf(x,y,z,'FaceColor','interp','EdgeColor','none')
rotate3d on
title(tit)
box on
pbaspect([1 1 .5])
end


function savbut(src,event) % manage data storings 
nedir=get(src,'Label');
switch nedir
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

function expofig (src,event) % manage image storing
figure2 = findobj('Type','figure','name','Action');
delete (figure2);
fig2 = figure('MenuBar','none','Name','Action',...
'NumberTitle','off','Resize','off',...
'Color','w',...
'units','normalized','outerposition',[.5 .1 .49 .9],...
'DockControls','off');
uimenu('Parent',fig2,'Label','PRINT','CallBack',@ploprint);
uimenu('Parent',fig2,'Label','RETURN','CallBack',@returproj);
set(fig2,'WindowButtonDownFcn', @selectoutmapimag)
axes('Parent',fig2,'units','normalized','Position',[0.15 0.1 0.75 0.75]);
text(0.2,0.6,{'1. Click on screen to select an output map/graphic';...
              '2. Use Print menu to save as image file...'})
axis off
end
function selectoutmapimag(src,event) % popup window for output selection
set(gcf,'outerposition',[.5 .1 .49 .9])
str = {'Basement Depth';'Inverted Anomaly';'Anomaly Difference';
    'Filter Construction';'RMS Plot';'Observed Anomaly'};
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
      showfiltercons
      case 5
      plot(rmstor,'-ko','MarkerFaceColor','r','MarkerSize',5);
      set(gca,'FontSize',10,'FontWeight','bold')
      sent1=['RMS(1)=' num2str(rmstor(1))];
      sent2=['RMS(end)=' num2str(rmstor(end))];
      title ({'RMS-Graph'; [sent1 '    ' sent2]})
      xlabel('Iteration number');ylabel('RMS');
      xlim([1 numel(rmstor)])
      grid on
      pbaspect([1 1 1]);axis square
      case 6
      mapper(x,y,T,'nT','Observed Mag.',1)
  end
 end
                  
end



function ploprint(src,event) % export as png image format
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

function returproj(src,event) % close and return to output GUI window
clc;clear all;clear global;
load('magbinv.mat');
%%%%%%%%%%%%%%close existing open figures
figure2 = findobj('Type','figure','name','Action');
delete (figure2);
figure3 = findobj('Type','figure','name','Cross-Section');
delete (figure3);
%%%%%%%%%%%%%%%%%%%%
%Secondfig
fig2 = figure('MenuBar','none','Name','Action',...
 'NumberTitle','off','Resize','off',...
 'Color','w',...
 'units','normalized','outerposition',[.5 .1 .49 .9],...
 'DockControls','off');
 set(fig2,'WindowButtonDownFcn', @selectoutmap)
 axes('Parent',fig2,'units','normalized','Position',[0.15 0.02 0.75 0.75]);
 axis off
psbut=uicontrol(fig2,'Style','pushbutton','units','normalized',...
'BackGroundColor',[.8 .8 .8],'HorizontalAlignment','right',...
'Position',[0 .95 1 0.03],'String','RMS Plot','Callback',@rmsplot);
mapper(x,y,-Zcalc,'Km','Calculated Depth',1)
menu1 = uimenu('Parent',fig2,'Label','SAVE/Export'); 
uimenu('Parent',menu1,'Label','Save as MAGBinv Project (*mat)','CallBack',@savbut);
menu2=uimenu('Parent',menu1,'Label','Export');
uimenu('Parent',menu2,'Label','RMS List (*.dat)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Calc-Z grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Calc-T grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Diff. (T-CalcT) grid (*.grd)','CallBack',@savbut);
uimenu('Parent',menu2,'Label','Export All','CallBack',@savbut);
uimenu('Parent',fig2,'Label','Save Plots as Images','CallBack',@savbut);
maxd=max(max(Zcalc));
 mind=min(min(Zcalc));
 meand=mean(mean(Zcalc));
 rs=numel(rmstor);
 uitable('Parent',fig2,'units','normalized','Position',[.02,.87,.97,.07],...
 'ColumnName',{'Max Depth','Min Depth','Mean Depth','Iteration Step'},...
 'Rowname',{'STAT'},...
 'ColumnWidth',{100,100,100,100},...
 'Enable','off','Data',[maxd mind meand rs]);
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


function savout(typ) % store data
load('magbinv.mat');
[filename, pathname] = uiputfile([filnam '.mgb'],...
'Set filename of data: an extension will be determined automatically according the data type');
kk=[pathname filename];
if ischar(kk)
kk=kk(1:end-4);    
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

function savmat %store as mat file
clc;clear all;clear global;
load('magbinv.mat');
load('magbinvset.mat');
[filename, pathname] = uiputfile([filnam '.mat'], ' Set a Project Name');
filename=filename(1:end-4);filename=[filename '_mgb.mat'];
kk=[pathname filename];
if ischar(kk)
save(kk,'x','y','T','Tcalc','Tdiff','Zcalc','rmstor','SH','WH','p','px',...
    'nx','ny','dx','dy','xmin','xmax','ymin','ymax','setto','filnam')
msgbox({'Project is Stored as *mgd.mat file...';...
    ' can be re-Open by Import Data > OPEN MAGBinv Project)'})
end
end

function selectoutmap(src,event) % output plot manager
set(gcf,'outerposition',[.5 .1 .49 .9])
str = {'Basement Depth';'Inverted Anomaly';'Anomaly Difference';
    'Filter Construction'};
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
      showfiltercons
      set(findobj(gcf,'Tag','plc1'),'enable','off')
      set(findobj(gcf,'Tag','plc2'),'enable','off')
      set(findobj(gcf,'Tag','plc3'),'enable','off')
      set(findobj(gcf,'Tag','plc4'),'enable','off')
  end
 end
                  
end

function showfiltercons % filter design plotter
load('magbinv.mat','p','px','SH','WH');
plot(px,p,'k','LineWidth',2);
line([WH WH],[min(p) max(p)],'Color','b','LineStyle','--','LineWidth',3)
line([SH SH],[min(p) max(p)],'Color','r','LineStyle','--','LineWidth',3)
pbaspect([1 .75 1])
grid on
title(['SH=' num2str(SH) '   WH=' num2str(WH)])
xlabel('k/2*pi');ylabel('Log(P)');
end

function rmsplot(src,event) % RMS plotter
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
set(findobj(gcf,'Tag','plc1'),'enable','off')
set(findobj(gcf,'Tag','plc2'),'enable','off')
set(findobj(gcf,'Tag','plc3'),'enable','off')
set(findobj(gcf,'Tag','plc4'),'enable','off')
end

function maininv(setto,fig2) %%% Kernel of main inversion
load('magbinv.mat','T','x','y','nx','ny','dx','dy','SH','WH');
freqorder=sort([SH WH]);WH=freqorder(1);SH=freqorder(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%get initials 
T0=T;
T=T./(10^9);
M=setto(5);
Im=setto(3)*pi/180;
If=setto(1)*pi/180;
Dm=setto(4)*pi/180;
Df=setto(2)*pi/180;
z0=setto(6);
alpha=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx0=nx;ny0=ny;
T(1,nx+floor(nx/2))=0;
T(ny+floor(ny/2),1)=0;
T=rot90(rot90(T));
T(1,nx+2*floor(nx/2))=0;
T(ny+2*floor(ny/2),1)=0;
T=rot90(rot90(T));
if (mod(nx,2)~=0) nx=nx-1; T(:,end)=[]; end
if (mod(ny,2)~=0) ny=ny-1; T(end,:)=[]; end
nxm=2*nx; nym=2*ny;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx= 2.*pi./((nxm-1).*dx);
dky= 2.*pi./((nym-1).*dy);
nyqx= (nxm/2)+1;
nyqy= (nym/2)+1; 
 for j=1:nxm
        if j <= nyqx
          kx(j)=(j-1)*dkx;
        else
          kx(j)=(j-(nxm+1))*dkx;
        end 
 for i=1:nym
        if i <= nyqy
          ky(i)=(i-1)*dky;
        else
          ky(i)=(i-(nym+1))*dky;
        end 
        k(i,j)=sqrt(kx(j).^2+ky(i).^2);
        if(k(i,j)==0) k(i,j)=0.0000000001; end;
 end
 end
 [kx ky]=meshgrid(kx,ky);
 cm=10^-7;
 mx=cos(Im)*sin(Dm-alpha);
 my=cos(Im)*cos(Dm-alpha);
 mz=sin(Im);
 fx=cos(If).*sin(Df-alpha);
 fy=cos(If).*cos(Df-alpha);
 fz=sin(If);

for m=1:nym
for n=1:nxm
    thetam(m,n)=mz+(1i).*((mx.*kx(m,n)+my.*ky(m,n))./abs(k(m,n)));
    thetaf(m,n)=fz+(1i).*((fx.*kx(m,n)+fy.*ky(m,n))./abs(k(m,n)));
    hs(m,n)=M*2.*pi.*cm.*thetam(m,n).*thetaf(m,n).*exp(-k(m,n).*z0);
end
end
%%%%%%%%%%%%%%%%%%%%            
ktotal=k./(2*pi); 
for f=1:nym;
   for g=1:nxm;
      if ktotal(f,g)<WH
      filter(f,g)=1;  
elseif ktotal(f,g)<SH
      filter(f,g)=0.5.*(1+cos((((2*pi)*ktotal(f,g))-(2*pi*WH))/(2*(SH-WH))));
else
filter(f,g)=0;
end
   end;
end;
FT=fft2(T);
Fh=FT./(hs.*(-k));
Fh=Fh.*filter;
Fh(1,1)=0;
h=real(ifft2(Fh)); h0=h;
h_old=h; %%%first approximation
%%%%%%%%%%%%%%%%%%%%% Iteration RMS Convergence (Stop RMS < threshold)
if setto(7)==0;
rms=1000;
rmstor=1000; % 
criterio=setto(8);% threshold value
m=2;
%%%%%%%%%%%create a waitbarr
steps = setto(9);
psbut=uicontrol(fig2,'Style','pushbutton','units','normalized',...
'BackGroundColor','r','HorizontalAlignment','right',...
'Position',[0 .95 1/steps 0.03],'Callback',@rmsplot);
%%%%%%%%%%%
 while rms>criterio % h_old needs to change...
     drawnow
     set(psbut,'Position',[0 .95 (m-1)/steps 0.03],'String',num2str(m-1));
     Fh=Fh-((-k).^(m-1)).*fft2(h.^m)./factorial(m);
     Fh=Fh.*filter;
     h=real(ifft2(Fh)); %new h in step
     dh=h-h_old; 
     dh2=dh.^2;
     rms=sqrt(sum(sum(dh2))./(nxm*nym)); % rms from new h
     rmstor(m)=rms;
     h_old=h; % if rms < criterio then this new h is our best result 
     if m-1==setto(9);break;end % step is equal maxiter then quit...
     m=m+1; %% non of criteria is achieved..new iteration starts
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RMS Divergence (Stop when RMSi>RMSi-1)
if setto(7)==1;
rms=1000;
rmstor(1)=1000; 
steps = setto(9);
%%%%%%%%%%% create a waitbarr
psbut=uicontrol(fig2,'Style','pushbutton','units','normalized',...
'BackGroundColor','r','Position',[0 .95 1/steps 0.03],...
'Callback',@rmsplot);
%%%%%%%%%%%
 for istep=2:steps+1
     m=istep;
     drawnow
     set(psbut,'Position',[0 .95 (m-1)/steps 0.03],'String',num2str(m-1))
     Fh=Fh-((-k).^(m-1)).*fft2(h.^m)./factorial(m);
     Fh=Fh.*filter;
     h=real(ifft2(Fh)); %new h in step m
     dh=h-h_old; 
     dh2=dh.^2;
     rms=sqrt(sum(sum(dh2))./(nxm*nym)); % new rms from new h
     if rms>rmstor(istep-1); break;end % % new h is not suitable..use h_old of best approximation
     rmstor(istep)=rms; %new rms is stored for next iteration
     h_old=h;
     if istep-1==setto(9);break;end % step is equal maxiter then quit...
     %%% h_old is updated 
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Plot results obtained from Convergence or divergence
%%%%%%%%%%%%%%%%%%%%%% inversion.
z=h_old+z0;
z=z(ny/2+1:ny/2+ny0,nx/2+1:nx/2+nx0);
%%%%%%%%%%% calculate forward mag from calculated z
Tcalc=forwardmag(z,nx0,ny0,dx,dy,M,Im,If,Dm,Df,alpha,z0);
Zcalc=z;
Tdiff=T0-Tcalc;
rmstor(1)=[];% delete the no comparing status
save('magbinv.mat','Tcalc','Zcalc','Tdiff','rmstor','-append');
drawnow
set(psbut,'String','preparing outputs..please wait','BackGroundColor',[.8 .8 .8],...
    'Position',[0 .95 1 0.03],'Enable','off');
drawnow
mapper(x,y,-z,'Km','Calculated Depth',1)
drawnow
set(psbut,'String','Plot RMS graph','BackGroundColor','g',...
    'Position',[0 .95 1 0.03],'Enable','on');
end

%%%%%%%%%%%%%%%%%%%%%%FORWARD MAG CALC %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tcalc=forwardmag(z,nx,ny,dx,dy,M,Im,If,Dm,Df,alpha,z0)
order=10;
z=z-z0;
%% extending
nx0=nx;ny0=ny;
z(1,nx+floor(nx/2))=0;
z(ny+floor(ny/2),1)=0;
z=rot90(rot90(z));
z(1,nx+2*floor(nx/2))=0;
z(ny+2*floor(ny/2),1)=0;
z=rot90(rot90(z));
if (mod(nx,2)~=0) nx=nx-1; z(:,end)=[]; end
if (mod(ny,2)~=0) ny=ny-1; z(end,:)=[]; end 
nxm=2*nx; nym=2*ny;
%% 
dkx= 2.*pi./((nxm-1).*dx);
dky= 2.*pi./((nym-1).*dy);
nyqx= (nxm/2)+1;
nyqy= (nym/2)+1; 
for j=1:nxm
       if j <= nyqx
         kx(j)=(j-1)*dkx;
       else
         kx(j)=(j-(nxm+1))*dkx;
       end 
for i=1:nym
       if i <= nyqy
         ky(i)=(i-1)*dky;
       else
         ky(i)=(i-(nym+1))*dky;
       end 
       k(i,j)=sqrt(kx(j).^2+ky(i).^2);
       if(k(i,j)==0) k(i,j)=0.00000001; end;
end
end
[kx ky]=meshgrid(kx,ky);
%%
cm=10^-7;
mx=cos(Im)*sin(Dm-alpha);
my=cos(Im)*cos(Dm-alpha);
mz=sin(Im);
fx=cos(If).*sin(Df-alpha);
fy=cos(If).*cos(Df-alpha);
fz=sin(If);
%%
for m=1:nym
for n=1:nxm
    thetam(m,n)=mz+(1i).*((mx.*kx(m,n)+my.*ky(m,n))./abs(k(m,n)));
    thetaf(m,n)=fz+(1i).*((fx.*kx(m,n)+fy.*ky(m,n))./abs(k(m,n)));
    hs(m,n)=M*2.*pi.*cm.*thetam(m,n).*thetaf(m,n).*exp(-k(m,n).*z0);
end
end
%%
SumF=0;
for m=1:order;
     SumF=SumF+((-k).^(m))./(factorial(m)).*fft2(z.^m);
end;
FT=hs.*SumF;
%FT(1,1)=0;
T=(ifft2(FT));
T=real(T);
T=T.*10^9;
Tcalc=T(ny/2+1:ny/2+ny0,nx/2+1:nx/2+nx0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function setact1(src,event) % activate Settings menu
nedir=get(src,'Label');
switch nedir
    case 'Set Initial Parameters'
    parametw
    case 'RMS Divergence (RMS(i+1)>RMS(i))'
    set(findobj(gcf,'Tag','c21'),'Checked','on')
    set(findobj(gcf,'Tag','c22'),'Checked','off')
    setto=get(findobj(gcf,'Tag','tabl1'),'Data');
    setto(7:8)=[1 0];
    save('magbinvset.mat','setto','-append');
    set(findobj(gcf,'Tag','tabl1'),'Data',setto);
    case 'RMS Convergence (RMS<threshold)'
    set(findobj(gcf,'Tag','c21'),'Checked','off')
    set(findobj(gcf,'Tag','c22'),'Checked','on')
    setcriterion
    case 'MaxIter'
    setmaxiter
    case 'Filter Observed Grid'
    setfilters    
end
end

function setfilters %%% activate filter design window
%%%Secondfig
figure2 = findobj('Type','figure','name','Action');
delete (figure2);
fig2 = figure('MenuBar','none','Name','Action',...
'NumberTitle','off','Resize','off',...
'Color','w',...
'units','normalized','outerposition',[.5 .1 .49 .9],...
'DockControls','off');
uicontrol('Parent',fig2,'Style','pushbutton','units','normalized',...
'Position',[.89,.96,.1,.03],...
'String','Help',...
'FontWeight','bold','CallBack',@positionshwh)
uicontrol('Parent',fig2,'Style','pushbutton','units','normalized',...
'Position',[0,.96,.15,.03],...
'String','Set SH',...
'FontWeight','bold','ForeGroundColor','r','CallBack',@positionshwh)
uicontrol('Parent',fig2,'Style','pushbutton','units','normalized',...
'Position',[.16,.96,.15,.03],...
'String','Set WH',...
'FontWeight','bold','ForeGroundColor','b','CallBack',@positionshwh)

uicontrol('Parent',fig2,'Style','pushbutton','units','normalized',...
'Position',[.32,.96,.15,.03],...
'String','Set Manually',...
'FontWeight','bold','ForeGroundColor','k','CallBack',@manushwh)

uicontrol('Parent',fig2,'Style','text','units','normalized',...
'Position',[.02,.89,.6,.04],'String',...
'Click on dashed lines to determine the current position of SH and WH',...
'BackGroundColor','w','ForeGroundColor','b','Tag','uyari',...
'HorizontalAlignment','left')

axes('Parent',fig2,'units','normalized','Position',[0.2 0.1 0.6 0.6]);
%%%%%current SHWH
load('magbinv.mat','p','px','SH','WH');
plot(px,p,'k','LineWidth',2);
line([WH WH],[min(p) max(p)],'Color','b','LineStyle','--','LineWidth',3,...
    'ButtonDownFcn',@shwhfig)
line([SH SH],[min(p) max(p)],'Color','r','LineStyle','--','LineWidth',3,...
    'ButtonDownFcn',@shwhfig)
pbaspect([1 1 1])
grid on
title ({'Radially Averaged Pow-Spectrum Graph';...
'Set SH and WH frequency band array for lowpass filtering'})
xlabel('k/2*pi');ylabel('Log(P)');
%%%%%
end
function shwhfig(src,event) % set current SH WH values as title
load('magbinv.mat','SH','WH');
title(['SH=' num2str(SH) '   WH=' num2str(WH)])
end

function manushwh(src,event) % Set SH WH by manually entering
load('magbinv.mat','p','px','SH','WH');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'SH'; 'WH'};
dlg_title = 'parameters';
num_lines = 1;
def = {num2str(SH);num2str(WH)};
answer= inputdlg(prompt,dlg_title,num_lines,def);
if numel(answer)>0
SH=str2num(cell2mat(answer(1)));if SH<0;SH=0;end
WH=str2num(cell2mat(answer(2)));if WH<0;WH=0;end
plot(px,p,'k','LineWidth',2);
line([WH WH],[min(p) max(p)],'Color','b','LineStyle','--','LineWidth',3,...
    'ButtonDownFcn',@shwhfig)
line([SH SH],[min(p) max(p)],'Color','r','LineStyle','--','LineWidth',3,...
    'ButtonDownFcn',@shwhfig)
pbaspect([1 1 1])
grid on
xlabel('k/2*pi');ylabel('Log(P)');
save('magbinv.mat','SH','WH','-append');
shwhfig
fg1 = findobj('Type','figure','name','MagB_inv');
set(findobj(fg1,'Tag','startbut'),'Enable','on')
end
end


function positionshwh(src,event) %%% set SH WH interactively
nedir=get(src,'String');
switch nedir
    case 'Set SH'
    set (gcf, 'WindowButtonMotionFcn', @mouseMove);
    set (gcf, 'WindowButtonDownFcn', @getmouseloc);
    set(gcf,'Pointer','fullcrosshair');
    load('magbinv.mat','p','px','SH','WH');
    plot(px,p,'k','LineWidth',2);
    line([WH WH],[min(p) max(p)],'Color','b','LineStyle','--','LineWidth',3)
    pbaspect([1 1 1])
    grid on
    title ('Set SH')
    fg1 = findobj('Type','figure','name','MagB_inv');
    set(findobj(fg1,'Tag','startbut'),'Enable','on')
    case 'Set WH'
    set (gcf, 'WindowButtonMotionFcn', @mouseMove);
    set (gcf, 'WindowButtonDownFcn', @getmouseloc);
    set(gcf,'Pointer','fullcrosshair');
    load('magbinv.mat','p','px','SH','WH');
    plot(px,p,'k','LineWidth',2);
    line([SH SH],[min(p) max(p)],'Color','r','LineStyle','--','LineWidth',3)
    pbaspect([1 1 1])
    grid on
    title ('Set WH')
    fg1 = findobj('Type','figure','name','MagB_inv');
    set(findobj(fg1,'Tag','startbut'),'Enable','on')
    case 'Help'
    set (gcf, 'WindowButtonMotionFcn', '');
    set (gcf, 'WindowButtonDownFcn', '');
    set(gcf,'Pointer','arrow');
    k=1:50;f=ones(1,50);sh=30;wh=20;
    f(20:30)=0.5.*(1+cos((((2*pi)*k(20:30))-(2*pi*wh))/(2*(sh-wh))));
    f(31:end)=0;
    plot(k,f,'--k','linewidth',3)
    pbaspect([1 1 1])
    grid on
    title('Filter configuration')
    hold on
    plot([20 20],[-0.5 1.5],'--b','Linewidth',2);
    text(16,-0.45,'WH','Color','b')
    plot([30 30],[-0.5 1.5],'--r','Linewidth',2)
    text(31,-0.45,'SH','Color','r')
    hold off
    xlabel('frequency')
    ylim([-0.5 1.5])
    set(gca,'xtick',[])
    text(10,1.2,'Full-Pass ','Color','b')
    text(35,0.2,'Reject','Color','r')
end
end

function getmouseloc(src,event) % mouse location finder regarding x axis
C = get (gca, 'CurrentPoint');
xlabel(['Selected frequency=' num2str(C(1,1))])
tit=get (gca, 'Title');
stit=get(tit,'String');
switch stit
    case 'Set SH'
    SH=C(1,1); if SH<0;SH=0;end
    save('magbinv.mat','SH','-append');
    load('magbinv.mat','p');
    line([SH SH],[min(p) max(p)],'Color','r','LineStyle','--','LineWidth',3)
    case 'Set WH'
    WH=C(1,1);if WH<0;WH=0;end
    save('magbinv.mat','WH','-append');
    load('magbinv.mat','p');
    line([WH WH],[min(p) max(p)],'Color','b','LineStyle','--','LineWidth',3)
end
set (gcf, 'WindowButtonMotionFcn', '');
set (gcf, 'WindowButtonDownFcn', '');
set(gcf,'Pointer','arrow');
shwhfig
end
function mouseMove(src,event) %% instant location info
C = get (gca, 'CurrentPoint');
xlabel(gca, ['Current-freq = ', num2str(C(1,1)) '  (k/2pi)']);
end


function setmaxiter %%% set of maximum iteration number
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
prompt = {'Maximum number of iteration...(in case criterion is not achieved)'};
dlg_title = 'Maxiter';
num_lines = 1;
def = {num2str(setto(9))};
answer= inputdlg(prompt,dlg_title,num_lines,def);
if numel(answer)>0
maxiter=str2num(cell2mat(answer(1)));maxiter=abs(maxiter);
setto(9)=round(maxiter);
save('magbinvset.mat','setto','-append');
set(findobj(gcf,'Tag','tabl1'),'Data',setto);
end
end

function setcriterion %%% manage termination criteria
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
prompt = {'set threshold RMS (Terminating if RMS<= threshold)'};
dlg_title = 'Criterion';
num_lines = 1;
def = {num2str(setto(8))};
answer= inputdlg(prompt,dlg_title,num_lines,def);
if numel(answer)>0
criteria=str2num(cell2mat(answer(1)));
setto(7:8)=[0 abs(criteria)];
save('magbinvset.mat','setto','-append');
set(findobj(gcf,'Tag','tabl1'),'Data',setto);
end
end

function parametw %%% Grid info panel
prompt = {'Inc. angle of ambient field (Incf)?'; 'Dec. angle of ambient field (Df)?';...
    'Inc. angle of magnetization(Im) ?'; 'Dec. angle of magnetization(Dm)?';...
'Magnetization Contrast (M)?';
'Average interface depth Z0 (unit in Km)?'};
dlg_title = 'parameters';
num_lines = 1;
setto=get(findobj(gcf,'Tag','tabl1'),'Data');
def = {num2str(setto(1));num2str(setto(2));num2str(setto(3));...
    num2str(setto(4));num2str(setto(5));num2str(setto(6))};
answer= inputdlg(prompt,dlg_title,num_lines,def);
if numel(answer)>0
If=str2num(cell2mat(answer(1)));
Df=str2num(cell2mat(answer(2)));
Im=str2num(cell2mat(answer(3)));
Dm=str2num(cell2mat(answer(4)));
M=str2num(cell2mat(answer(5)));
z0=abs(str2num(cell2mat(answer(6))));
setto(1:6)=[If Df Im Dm M z0];
save('magbinvset.mat','setto','-append');
set(findobj(gcf,'Tag','tabl1'),'Data',setto);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%radially averaged power spectrum
function [p,xlimits]=getraps(T,nx,ny,dx,dy);
T(1,nx+floor(nx/2))=0;
T(ny+floor(ny/2),1)=0;
T=rot90(rot90(T));
T(1,nx+2*floor(nx/2))=0;
T(ny+2*floor(ny/2),1)=0;
T=rot90(rot90(T));
if (mod(nx,2)~=0) nx=nx-1; T(:,end)=[]; end
if (mod(ny,2)~=0) ny=ny-1; T(end,:)=[]; end 
[p,~]=calcrad2(T);% radially averaged power spectrum
nyq = 1/(2*dx); % frequency axis 
p = log(smooth(p));
sca = length(p);
xlimits = linspace(0,nyq,sca);
end
%%%%%%%%%%%%%% New data import
%%%%%action menu 1 IMPORT 2-D (*.grd) Golden Software Binary/or Text grid 
function impgrd(src,event);
clc;
[filename, pathname] = uigetfile('*.grd', 'Import Golden Software Binary/Text grid (*.grd)');
k=[pathname filename];
if ischar(k)
fidc=fopen(k);
header= fread(fidc,4,'*char' )';
fclose(fidc);
c1=strcmp(header,'DSAA');
c2=strcmp(header,'DSRB');
sumc=sum([c1 c2]);
if sumc>0
%%delete temporary import file 
if exist('magbinv.mat','file')==2;delete('magbinv.mat');end
set(findobj(gcf,'Tag','mainfig'),'Pointer','watch')
switch c1
    case 1
[T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k);
    case 0
[T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd7bin(k);        
end
if dx~=dy;msgbox('grid interval dx=dy required...');return;end
[p,px]=getraps(T,nx,ny,dx,dy);
SH=0;
WH=0;
filnam=filename(1:end-4);
list1w(T,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k);
save('magbinv.mat','T','x','y','nx','ny','xmin','xmax','ymin','ymax','dx','dy',...
    'p','px','SH','WH','filnam');
%fileattrib('magbinv.mat','+h');%hidden temporary import info file
drawnow
mapper(x,y,T,'nT','Observed Anomaly',1)
set(findobj(gcf,'Tag','mainfig'),'Pointer','arrow')
figure2 = findobj('Type','figure','name','Action');
delete (figure2);
figure3 = findobj('Type','figure','name','Cross-Section');
delete (figure3);
set(findobj(gcf,'Tag','c3'),'Enable','on')
set(findobj(gcf,'Tag','startbut'),'Enable','off')
else
msgbox('File format not supported..Load Surfer 6 text grid / Surfer 7 Binary Grid')
end
end
end
%%%%%action menu 1 IMPORT a XYZ (*.dat) equal spaced columnwise grid file 
function impgrdat(src,event);
clc;
[filename, pathname] = uigetfile('*.dat', 'Import XYZ grid data(*.dat)');
k=[pathname filename];
if ischar(k)
if exist('magbinv.mat','file')==2;delete('magbinv.mat');end
set(findobj(gcf,'Tag','mainfig'),'Pointer','watch')
[T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrdatfile_Callback(k);
[p,px]=getraps(T,nx,ny,dx,dy);
SH=0;
WH=0;
filnam=filename(1:end-4);
list1w(T,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k)
save('magbinv.mat','T','x','y','nx','ny','xmin','xmax','ymin','ymax','dx','dy',...
    'p','px','SH','WH','filn');
%fileattrib('magbinv.mat','+h');%hidden temporary import info file
drawnow
mapper(x,y,T,'nT','Observed Anomaly',1)
set(findobj(gcf,'Tag','mainfig'),'Pointer','arrow')
figure2 = findobj('Type','figure','name','Action');
delete (figure2);
figure3 = findobj('Type','figure','name','Cross-Section');
delete (figure3);
set(findobj(gcf,'Tag','c3'),'Enable','on')
set(findobj(gcf,'Tag','startbut'),'Enable','off')
end
end

%%%%%read Surfer 6 text grid(*.grd)
function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k)
surfergrd=fopen(k,'r'); % Open *.grid file
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read Surfer 7 Binary grid
function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy] = lodgrd7bin(filename)
fid= fopen(filename);
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
%%%%%%%%%%read XYZ gridded data (*.dat)
function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrdatfile_Callback(k)
data=load(k);
x=data(:,1);y=data(:,2);T=data(:,3); 
unix=unique(x);xmin=unix(1);xmax=unix(end);dx=unix(2)-unix(1);
uniy=unique(y);ymin=uniy(1);ymax=uniy(end);dy=uniy(2)-uniy(1);
nx=numel(unix);
ny=numel(uniy);
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

function mapper(x,y,matrix,unitt,tit,index) %%map viewer
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


%%%%%%%%%%% Function for output of a GRID
function grdout(matrix,xmin,xmax,ymin,ymax,namefile)
%Get grid dimensions
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

%%%Grid info List Box
function list1w(T,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k)
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

function [Prad,w]=calcrad2(X)
% Calculates the radially-averaged power spectrum
% (RAPS) of the data in the array X
% Notes: Uses pergram2.m and fcoef
% By Thomas A. Ridsdill-Smith (March 2000)
% Ridsdill-Smith, T.A., 2000. The application of the wavelet transform to the processing
% of aeromagnetic data. Ph. D. Dissertation, The University of Western Australia,
% Australia, 197pp.
sz=size(X);
nmin=min(sz);
fX=fft2(X);
% Calculate periodogram
P=pergram2(fX,sz(1),sz(2));
% Radially average
w=fcoef1(nmin);
w=w(1:floor(nmin/2)+1);
wbin=round((nmin*w)./(2*pi));
for k=0:(0.5*nmin) % Don't want to go past Nyquist
   Prad(k+1)=mean(mean(P(find(wbin==k))));
end

end

function [P,r]=pergram2(fx,M,N)
% Calculates the periodogram P of the discrete Fourier 
% transform fx at the frequencies fi and fj (2D Algorithm)
% By Thomas A. Ridsdill-Smith (March 2000)
% Ridsdill-Smith, T.A., 2000. The application of the wavelet transform to the processing
% of aeromagnetic data. Ph. D. Dissertation, The University of Western Australia,
% Australia, 197pp.

% Setup frequencies
fi=[0:(floor(M/2))].*(2*pi/M);
fj=[0:(floor(N/2)) (-(floor(N/2))+1):-1].*(2*pi/N);
[xi,yi]=meshgrid(fj,fi);
r=sqrt(xi.^2+yi.^2);
% Initialise P
P=zeros((floor(M/2))+1,N);
% Calculate power
modfx2=(abs(fx)).^2;
for j=1:N
   P(1,j)=modfx2(1,j);             % Zero wavenumber (row)
   P(floor(M/2)+1,j)=modfx2(floor(M/2)+1,j); % Nyquist
   for i=2:(floor(M/2))
      P(i,j)=modfx2(i,j)+modfx2(M-i+2,j);
      % Sum positive and negative frequencies
   end
end

%P=(1/((N*M)^2))*P; %normalize or not
P=P';
P = P(:,2:end);
end

function [w]=fcoef1(n)
% Calculates the 1D Fourier angular frequencies
% (in radians) for a signal of length n
% (eg w=[(0:(n/2))(((-n/2)+1):-1)].*(2*pi/n) )
% By Thomas A. Ridsdill-Smith (March 2000)
% Ridsdill-Smith, T.A., 2000. The application of the wavelet transform to the processing
% of aeromagnetic data. Ph. D. Dissertation, The University of Western Australia,
% Australia, 197pp.

if (rem(n,2)==0)
   % Even signal length
   w=[(0:(n/2)) (((-n/2)+1):-1)].*(2*pi/n);
else
   % Odd signal length
   w=[(0:((n-1)/2)) ((-(n-1)/2):-1)].*(2*pi/(n-1));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createAnom (src,event)
figure2 = findobj('Type','figure','name','Action');
delete (figure2);
fig2 = figure('MenuBar','none','Name','Action',...
'NumberTitle','off','Resize','off',...
'Color','w',...
'units','normalized','outerposition',[.5 .1 .49 .9],...
'DockControls','off');
uimenu('Parent',fig2,'Label','Load Depth Grid','CallBack',@impordep);
uimenu('Parent',fig2,'Label','Settings & Calc','CallBack',@calcmag);
sm=uimenu('Parent',fig2,'Label','Plot','Tag','swiftm','Enable','off');
uimenu('Parent',sm,'Label','Depth Model','CallBack',@switcmap1);
uimenu('Parent',sm,'Label','Calc-Mag','CallBack',@switcmap2);
set(fig2,'WindowButtonDownFcn', @saveforward)
axes('Parent',fig2,'units','normalized','Position',[0.15 0.1 0.75 0.75]);
drawnow
plot(.5,.5,'w+')
text(0.0,0.8,{'1. Load depth model grid file (*.grd)';' ';' ';...
              '2. Set model parameters....'; '  ';' ';...
              '3. Click on desired maps for exporting as image file and data'})
axis off
axis off
end

function impordep(src,event)
[filename, pathname] = uigetfile('*.grd', 'Import Golden Software Binary/Text grid (*.grd)');
k=[pathname filename];
if ischar(k)
fidc=fopen(k);
header= fread(fidc,4,'*char' )';
fclose(fidc);
c1=strcmp(header,'DSAA');
c2=strcmp(header,'DSRB');
sumc=sum([c1 c2]);
if sumc>0
switch c1
    case 1
[z,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k);
    case 0
[z,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd7bin(k);        
end
z=abs(z);
if dx~=dy;msgbox('grid interval dx=dy required...');return;end
mpstyle=1;
save('magbinvForw.mat','z','x','y','nx','ny','dx','dy','xmin','xmax','ymin','ymax','mpstyle');
drawnow
plot(.5,.5,'k+')
axis off
mapper(x,y,-z,'Km','Depth Model',1)
set(findobj(gcf,'Tag','swiftm'),'Enable','off');
end
end
end

function calcmag(src,event)
prompt = {'Inc. angle of ambient field (Incf)?'; 'Dec. angle of ambient field (Df)?';...
    'Inc. angle of magnetization(Im) ?'; 'Dec. angle of magnetization(Dm)?';...
'Magnetization Contrast (M)?';
'Average interface depth Z0 (unit in Km)?'};
dlg_title = 'parameters';
num_lines = 1;
def = {'';'';'';'';'2';'15'};
answer= inputdlg(prompt,dlg_title,num_lines,def);
if numel(answer)>0
If=str2num(cell2mat(answer(1)))*pi/180;
Df=str2num(cell2mat(answer(2)))*pi/180;
Im=str2num(cell2mat(answer(3)))*pi/180;
Dm=str2num(cell2mat(answer(4)))*pi/180;
M=str2num(cell2mat(answer(5)));
z0=abs(str2num(cell2mat(answer(6))));
load('magbinvForw.mat');
alpha=0;
Tcalc=forwardmag(z,nx,ny,dx,dy,M,Im,If,Dm,Df,alpha,z0);
drawnow
plot(.5,.5,'k+')
axis off
mapper(x,y,Tcalc,'nT','Calculated Anomaly From Depth Model',1)
mpstyle=2;
save('magbinvForw.mat','Tcalc','mpstyle','-append')
set(findobj(gcf,'Tag','swiftm'),'Enable','on');
end
end

function saveforward(src,event)
clc;clear all;clear global
[filename, pathname] = uiputfile('*.png','Set a filename');
kk=[pathname filename];
if ischar(kk)
print(kk, '-dpng', '-r500');
drawnow
plot(.5,.5,'w+')
axis off
load('magbinvForw.mat');

if mpstyle==2;
grdout(Tcalc,xmin,xmax,ymin,ymax,[kk '-Tcalc.grd']);
text(.1,.6,'Calculated magnetic anomaly is exported as image and *.grd file',...
'Color','b')
else
text(.1,.6,' Depth model is exported as image file','Color','b')
end
end
end

function switcmap1 (src,event)
drawnow
plot(.5,.5,'k+')
axis off
load('magbinvForw.mat');
mapper(x,y,-z,'Km','Depth Model',1)
mpstyle=1;
save('magbinvForw.mat','mpstyle','-append')
end

function switcmap2 (src,event)
drawnow
plot(.5,.5,'k+')
axis off
load('magbinvForw.mat');
mapper(x,y,Tcalc,'nT','Calculated Anomaly From Depth Model',1)    
mpstyle=2;
save('magbinvForw.mat','mpstyle','-append')
end