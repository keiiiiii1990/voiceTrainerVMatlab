function fig=training_gui(nbfiles,nbexpe,run_all)

% TRAINING_GUI Creates the GUI for the training phase
%
% fig=training_gui(nbfile,nbexpe,run_all) creates the GUI given the number
% of tested files nbfile (including the reference), the number of
% experiments nbexpe and the flag run_all describing whether the evaluation
% phase is to be run afterwards
%
% nbfile: �����ł�1�����ɂ�����t�@�C����
% nbexpe: ������
% run_all:
% 
% 
% returns a handle to the figure if succeeds (nbfile and nbexpe small
% enough), returns false otherwise

%%%getting the screen size and computing the figure width/height
%screen size

nbfile = max(nbfiles);

set(0,'Units','characters');
siz=get(0,'Screensize');
maxwidth=siz(3);
maxheight=siz(4)*.92;
%width of each item
%iwidth=min(11,(maxwidth-48)/(nbfile-1));
iwidth=min(11,(maxwidth-48)/nbfile);
%{
if iwidth < 9.2,
    fig=false;
    return;
end
%}
%figure width
%width=48+iwidth*(nbfile-1);
width=48+iwidth*(nbfile);
%height of each item
iheight=min(2.8,(maxheight-12.3)/nbexpe);
if iheight < 2,
    fig=false;
    return;
end
%figure height
height=12.3+iheight*nbexpe;

%%%opening the figure
fig=figure('Name','MUSHRAM - Training phase','NumberTitle','off','MenuBar','none','Resize','off','Color',[0.701961 0.701961 0.701961],'Units','characters','Position',[0 0 width height]);

%%%displaying fixed items
%titles
uicontrol(fig,'Style','Text','Units','characters','Position',[0 9.3+nbexpe*iheight width 1.5],'FontSize',12,'String','Training','Tag','training');
%uicontrol(fig,'Style','Text','Units','characters','Position',[26.5 6.3+nbexpe*iheight 11 1],'FontSize',10,'String','Reference','Tag','reference');
%uicontrol(fig,'Style','Text','Units','characters','Position',[46.5 6.3+nbexpe*iheight width-49 1],'FontSize',10,'String','Test','Tag','test');
uicontrol(fig,'Style','Text','Units','characters','Position',[26.5 6.3+nbexpe*iheight width-49 1],'FontSize',10,'String','Test','Tag','test');
%exit or proceed button
proceed=uicontrol(fig,'Style','Pushbutton','Units','characters','FontSize',10,'Tag','proceed','Callback','training_callbacks(''proceed'',guidata(gcbo))');
if run_all,
    set(proceed,'Position',[width/2-13.5 1.5 27 1.8],'String','Proceed to evaluation');
else,
    set(proceed,'Position',[width/2-4.5 1.5 9 1.8],'String','Exit');
end

%%%displaying items depending on the number of experiments and test signals
%handles = guidata(gcbo);
for e=1:nbexpe,
    %experiment number
    uicontrol(fig,'Style','Text','Units','characters','Position',[1.5 6.1+(nbexpe-e)*iheight 15.5 1],'FontSize',10,'String',['Experiment ' int2str(e)],'Tag',['title' int2str(e)]);
    %reference play buttons
    %uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[22.5 5.7+(nbexpe-e)*iheight 19 1.8],'FontSize',10,'String','Play reference','Callback',['training_callbacks(''play'',guidata(gcbo),' int2str(e) ',0)'],'Tag',['play' int2str(e) '_0']);
    for f=1:nbfile,
        %play buttons
        %uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[46.5+(f-1)*iwidth 5.7+(nbexpe-e)*iheight 9 1.8],'FontSize',10,'String','Play','Callback',['training_callbacks(''play'',guidata(gcbo),' int2str(e) ',' int2str(f) ')'],'Tag',['play' int2str(e) '_' int2str(f)]);
        %{
        if ~isempty(handles.files{handles.expe_order(e),handles.file_order(e,f)})
            uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[26.5+(f-1)*iwidth 5.7+(nbexpe-e)*iheight 9 1.8],'FontSize',10,'String','Play','Callback',['training_callbacks(''play'',guidata(gcbo),' int2str(e) ',' int2str(f) ')'],'Tag',['play' int2str(e) '_' int2str(f)]);
        else
            uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[26.5+(f-1)*iwidth 5.7+(nbexpe-e)*iheight 9 1.8],'FontSize',10,'String','Play','Callback',['training_callbacks(''play'',guidata(gcbo),' int2str(e) ',' int2str(f) ')'],'Tag',['play' int2str(e) '_' int2str(f)],'enable','off');
        end
        %}
        uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[26.5+(f-1)*iwidth 5.7+(nbexpe-e)*iheight 9 1.8],'FontSize',10,'String','Play','Callback',['training_callbacks(''play'',guidata(gcbo),' int2str(e) ',' int2str(f) ')'],'Tag',['play' int2str(e) '_' int2str(f)]);
    end
end

%%%moving onscreen
movegui(fig,'onscreen');