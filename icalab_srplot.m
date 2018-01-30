%Draw Estimated soueces (independent components) figure 
%Improved by T.Watanabe 07.05.2003

function h = icalab_srplot( par, action)


global sFontSize backGroundColor
global plotoptionb
% sFontSize = 8;
ICALAB_PlotTreshold = 10e-3;
ICALAB_tmp_post = '.      ';


defposi=[0.1 0.1 0.65 0.8];
if nargin == 1,
    if ~isfield(par, 'x_limit')
        par.x_limit = [];
    end
    
    if ~isfield(par, 'x_title')
        par.x_title = [];
    end
    
    NoOfSign_in_par = size(par.signal,1); % liczba sygnalow w par.signal
      

        
        for kk=1:NoOfSign_in_par   % EEG_params.ch_num,
            par.labels{kk}=[par.sname int2str(kk)];
        end

    %----------------------------------------------------
    [sensorsno,T] = size(par.signal); % par.plotno - sensors number, T - signals length
    par.siglength=T;
    %Initialize xwindowsize
    if ~isfield(par,'xwindowsize')
        if par.siglength > 2000
            par.xwindowsize=ceil(par.siglength/2);
        else
            par.xwindowsize=par.siglength;
        end
    end
    %---------------------
    %Window Xlim Size
    if par.siglength > par.xwindowsize
        par.xlimsize=par.xwindowsize;
    else
        par.xlimsize=par.siglength;
    end
    %---------------------
    main_fig = findobj( 'Tag', par.tag );
    if ~isempty( main_fig )
        close(main_fig);
    end
    %---------------------
    % Create Figure
    %---------------------
    h = figure( 'Visible', 'off');
    par.window_number = h;
    if strcmp( par.sposition, 'auto' ) == 1,
        set( h, 'Units', par.sunits );
    else
        set( h, 'Units', par.sunits,...
            'OuterPosition', par.sposition );
    end
    set( h, 'Tag', par.tag,...
        'Name', par.stitle,...
        'NumberTitle', 'off',...
        'ToolBar','figure',...
        'Color', backGroundColor );
   
         
    %---------------------
    % Plot Signals
    %---------------------
    sstart=1;
    if sensorsno <= par.num,
        par.plotno=sensorsno;
    else
        par.plotno=par.num;
    end
  
    
    %-------------------------------------------------   
    if par.plotno < sensorsno
        mslider = uicontrol( 'Parent', h, ...
            'Units', 'normalized', ...
            'BackgroundColor', backGroundColor, ...
            'ListboxTop', 0, ...
            'Position', [defposi(1)+defposi(3)+0.02   defposi(2)   0.03   defposi(4)], ...
            'Callback', 'icalab_srplot( [], ''slider_moved'' );', ...
            'Style', 'slider', ...
            'Tag', ['MojSlider_' int2str(par.window_number)] );
        
        set( mslider, 'Min', 0,...
            'Max', sensorsno-par.plotno ,...
            'Value', sensorsno-par.plotno );
        if par.plotno == 1
            set( mslider, 'SliderStep', [1/(sensorsno-par.plotno ) (1+par.plotno )/(sensorsno-par.plotno )] );
        else
            set( mslider, 'SliderStep', [1/(sensorsno-par.plotno ) par.plotno /(sensorsno-par.plotno )] );
        end
    end
    %-----------------------------------------------------
    %  Indicator for Paging
    %-----------------------------------------------------
    par.page=1;
    par.allpage=1;
    %-------------------------
    %indicator
    wide=0.05;
    vm.h_axes= subplot('position',[defposi(1) defposi(2)+defposi(4)+0.045 defposi(3) wide-0.01]);
    vm.map=image([35]);
    set(vm.map,'xdata',par.page-0.5);
    par.allpage=ceil(par.siglength/par.xlimsize);
    vm.text=text(par.page-0.5, 1,[num2str(par.page) '/' num2str(par.allpage)]);
    set(vm.text,'FontWeight','bold',...
        'HorizontalAlignment','center')
    
    if par.siglength >= par.xlimsize
        
        hold on
        for k=1:par.allpage-1,
            plot([k k],[0.5 1.5],'k');
        end
        axis([0 par.siglength/par.xlimsize 0.5 1.5 ]);
        set(vm.h_axes,...
            'XTickLabel', [1 par.siglength],...
            'YTickLabel', [] ,...
            'XTick',[0 par.siglength/par.xlimsize],...
            'YTick',[],...
            'Fontsize',sFontSize,...
            'ButtonDownFcn','ha=get(gca,''Currentpoint'');icalab_srplot(ha(1,1),''indicator_move'');');
        hold off
        %-------------------------
        %push button
        sp1=0.075;
        hi=0.05;
        hispace=0.08;
        wide2=0.1;
        defposi1=[defposi(1)+defposi(3)+sp1 defposi(2)+defposi(4) 0.9-defposi(1)+defposi(3)+sp1];
        h_text_paging=uicontrol('style','text',...
            'Units','normalized',...
            'String','Paging',...
            'fontweight','bold',...
            'position',[defposi1(1) defposi1(2) wide2 hi]);
        par.push_pagedw=uicontrol('style','pushbutton',...
            'Units','normalized',...
            'FontWeight','bold',...
            'position',[defposi1(1) defposi1(2)-hispace wide2/2 hi],...
            'string','<<',...
            'callback','icalab_srplot(''down'',''indicator_move'');');
        par.push_pageup=uicontrol('style','pushbutton',...
            'Units','normalized',...
            'FontWeight','bold',...
            'position',[defposi1(1)+wide2/2 defposi1(2)-hispace wide2/2 hi],...
            'string','>>',...
            'callback','icalab_srplot(''up'',''indicator_move'');');
        h_text=uicontrol('style','text',...
            'Units','normalized',...
            'String','pt',...
            'position',[defposi1(1)+wide2 defposi1(2)-hispace*2 wide2/2 hi]);
        
        par.edit_paging=uicontrol('style','edit',...
            'Units','normalized',...
            'String',num2str(par.xlimsize),...
            'backgroundcolor',[1 1 1],...
            'position',[defposi1(1) defposi1(2)-hispace*2 wide2 hi],...
            'callback','icalab_srplot([],''change_xwindowsize'');');
        sc.h_text_scaling=uicontrol('style','text',...
            'Units','normalized',...
            'String','Scaling',...
            'fontweight','bold',...
            'position',[defposi1(1) defposi1(2)-hispace*4 wide2 hi]);
        sc.h_radio_scaling_auto=uicontrol('style','radiobutton',...
            'Units','normalized',...
            'String','automatic',...
            'position',[defposi1(1) defposi1(2)-hispace*5 wide2*1.5 hi],...
            'value',1,...
            'callback',[...
                'sc=getappdata(gcf,''ScalingData'');'...
                'set(sc.h_radio_scaling_auto,''value'',1);'...
                'set(sc.h_radio_scaling_fixed,''value'',0);'...
                'set(sc.h_text_scaling_min,''enable'',''off'');'...
                'set(sc.h_text_scaling_max,''enable'',''off'');'...
                'set(sc.h_edit_scaling_min,''enable'',''off'',''string'',''N/A'');'...
                'set(sc.h_edit_scaling_max,''enable'',''off'',''string'',''N/A'');'...
                'sc.type=''auto'';'...
                'setappdata(gcf,''ScalingData'',sc);'...
                'icalab_srplot([],''slider_moved'');'...
            ]);
        sc.h_radio_scaling_fixed=uicontrol('style','radiobutton',...
            'Units','normalized',...
            'String','fixed',...
            'position',[defposi1(1) defposi1(2)-hispace*6 wide2*1.5 hi],...
            'callback',[...
                'sc=getappdata(gcf,''ScalingData'');'...
                'set(sc.h_radio_scaling_auto,''value'',0);'...
                'set(sc.h_radio_scaling_fixed,''value'',1);'...
                'signalylimdata=sc.signalylimdata;'...
                'sc.ylim1=signalylimdata(1);sc.ylim2=signalylimdata(2);'...
                'set(sc.h_text_scaling_min,''enable'',''on'');'...
                'set(sc.h_text_scaling_max,''enable'',''on'');'...
                'set(sc.h_edit_scaling_min,''enable'',''on'',''string'',sc.ylim1);'...
                'set(sc.h_edit_scaling_max,''enable'',''on'',''string'',sc.ylim2);'...
                'sc.type=''fixed'';'...
                'setappdata(gcf,''ScalingData'',sc);'...
                'icalab_srplot([],''slider_moved'');'...
            ]);
        sc.h_text_scaling_min=uicontrol('style','text',...
            'enable','off',...
            'Units','normalized',...
            'String','min:',...
            'position',[defposi1(1) defposi1(2)-hispace*7 wide2 hi]);
        sc.h_edit_scaling_min=uicontrol('style','edit',...
            'enable','off',...
            'Units','normalized',...
            'String','N/A',...
            'backgroundcolor',[1 1 1],...
            'position',[defposi1(1) defposi1(2)-hispace*8 wide2 hi],...
            'callback','icalab_srplot([],''edit_set_min_scale'');');
        sc.h_text_scaling_max=uicontrol('style','text',...
            'enable','off',...
            'Units','normalized',...
            'String','max:',...
            'position',[defposi1(1) defposi1(2)-hispace*9 wide2 hi]);
        sc.h_edit_scaling_max=uicontrol('style','edit',...
            'enable','off',...
            'Units','normalized',...
            'String','N/A',...
            'backgroundcolor',[1 1 1],...
            'position',[defposi1(1) defposi1(2)-hispace*10 wide2 hi],...
            'callback','icalab_srplot([],''edit_set_max_scale'');');
        sc.ylim1=[];
        sc.ylim2=[];
        sc.type='auto';
        setappdata(h,'ScalingData',sc);
    end
    if isequal(par.siglength,par.xlimsize)
        set(par.push_pagedw,'enable','off');
        set(par.push_pageup,'enable','off');
        set(vm.h_axes,'visible','off');
        set(vm.map,'visible','off');
        set(vm.text,'visible','off');
    else
        set(par.push_pagedw,'enable','on');
        set(par.push_pageup,'enable','on');
        set(vm.h_axes,'visible','on');
        set(vm.map,'visible','on');
        set(vm.text,'visible','on');
    end
    
    setappdata(h,'IndicatorData',vm);
    %-----------------------------------------------------
    set( h, 'Visible', 'on' );
    set( h,'UserData', par,'Visible', 'on' );
    icalab_srplot([],'slider_moved');
    %----------------------------------------------------------
    % Actions CallbackFcn
    %----------------------------------------------------------
    
   
else
       
    
    page=par;
    h = figure( gcf );
    par = get( h, 'UserData' );
    %-----------------------------------------------------
    %ylim setting
    %-----------------------------------------------------
    if strcmp(action,'edit_set_min_scale')
        sc=getappdata(gcf,'ScalingData');
        val=get(sc.h_edit_scaling_min,'string');
        val=str2double(val);
        if isnan(val)
            set(sc.h_edit_scaling_min,'string',sc.ylim1);
            return
        else
            if val>=sc.ylim2
                val=sc.ylim2-abs(sc.ylim2-sc.ylim1)*0.1;
            end
            set(sc.h_edit_scaling_min,'string',val);
            sc.ylim1=val;
            setappdata(gcf,'ScalingData',sc);
            icalab_srplot([],'slider_moved');
        end            
    elseif strcmp(action,'edit_set_max_scale')
        sc=getappdata(gcf,'ScalingData');
        val=get(sc.h_edit_scaling_max,'string');
        val=str2double(val);
        if isnan(val)
            set(sc.h_edit_scaling_max,'string',sc.ylim2);
            return
        else
            if val<=sc.ylim2
                val=sc.ylim1+abs(sc.ylim2-sc.ylim1)*0.1;
            end
            set(sc.h_edit_scaling_max,'string',val);
            sc.ylim2=val;
            setappdata(gcf,'ScalingData',sc);
            icalab_srplot([],'slider_moved');
        end            
        
        %-----------------------------------------------------
        % Callback 'plotoption' Option Button on main icalab figure
        %-----------------------------------------------------
    elseif strcmp(action,'plotoption');
        figdata.h=h;
        optiondata=get(gcbo,'userdata');
        if isempty(optiondata)
            check1=1;
            check2=0;
            str1='[0,1]';
            str2='Normalized Frequency';
        else
            check1=0;
            check2=1;
            str1=optiondata.str1;
            str2=optiondata.str2;
        end
        optionh=dialog('name','Plot Option','Units','characters','visible','off');
        posi1=get(optionh,'position');
        posi1=[posi1(1)+posi1(3)/2-30 posi1(2)+posi1(4)/2-15 60 15];
        strhi=1.3;
        edithi=1.7;
        wi=50;
        
        set(optionh,'position',posi1);
        
        figdata.checkhdl1=uicontrol(optionh,...
            'units','characters',...
            'style','radiobutton',...
            'position',[5 12 wi edithi],...
            'string','point',...
            'value',check1,...
            'callback',['datas=get(gcbo,''userdata'');'...
                'dataa=datas{2};set(dataa(1),''value'',0);set(dataa(2:5),''enable'',''off'');']);
        figdata.checkhdl2=uicontrol(optionh,'style','radiobutton',...
            'units','characters',...
            'position',[5 10 wi edithi],...
            'string','manual',...
            'value',check2,...
            'callback',['datas=get(gcbo,''userdata'');dataa=datas{1};set(dataa(2:5),''enable'',''on'');' ...
                'set(datas{2},''value'',0);']);
        figdata.texthdl1=uicontrol(optionh,'style','text',...
            'units','characters',...
            'position',[5 8 15 strhi],...
            'string','X axis limit:',...
            'HorizontalAlignment','left');
        figdata.texthdl2=uicontrol(optionh,'style','text',...
            'units','characters',...
            'position',[5 6 15 strhi],...
            'string','X axis label:',...
            'HorizontalAlignment','left');
        figdata.edithdl1=uicontrol(optionh,'style','edit',...
            'units','characters',...
            'position',[20 8 35 edithi],...
            'string',str1,...
            'HorizontalAlignment','left',...
            'BackgroundColor',[1 1 1]);
        figdata.edithdl2=uicontrol(optionh,'style','edit',...
            'units','characters',...
            'position',[20 6 35 edithi],...
            'string',str2,...
            'HorizontalAlignment','left',...
            'BackgroundColor',[1 1 1]);
        pushhdl1=uicontrol(optionh,'style','pushbutton',...
            'units','characters',...
            'position',[20 2 12 2],...
            'string','OK',...
            'FontWeight','bold',...
            'callback','icalab_srplot([],''plotoption_ok'');');
        pushhdl1=uicontrol(optionh,'style','pushbutton',...
            'units','characters',...
            'position',[43 2 12 2],...
            'string','Cancel',...
            'FontWeight','bold',...
            'callback','delete(gcbf);');
        dataa={[figdata.checkhdl2;figdata.texthdl1;...
                    figdata.texthdl2;figdata.edithdl1;figdata.edithdl2]};
        datab={[figdata.checkhdl1]};
        set(figdata.checkhdl1,'userdata',[datab;dataa]);
        set(figdata.checkhdl2,'userdata',[dataa;datab]);
        set(optionh,'visible','on','userdata',figdata);
        
        if check1 > 0
            set(figdata.edithdl1,'enable','off');
            set(figdata.edithdl2,'enable','off');
            set(figdata.texthdl1,'enable','off');
            set(figdata.texthdl2,'enable','off');
        end
        
        %-----------------------------------------------------
        % Callback 'plotoption' Option Button on main icalab figure
        %-----------------------------------------------------
    elseif strcmp(action,'plotoption_ok');
        figdata=get(gcbf,'userdata');
        check1=get(figdata.checkhdl1,'value');
        check2=get(figdata.checkhdl2,'value');
        str1=get(figdata.edithdl1,'string');
        str2=get(figdata.edithdl2,'string');
        if check1>0
            delete(gcbf);
            figure(figdata.h);
            set(plotoptionb,'userdata',[]);
        else
            if isempty(str1)
                errordlg('X axis limit is incorrect');
                return
            else
                val1=str2num(str1);
                if ~(isequal(size(val1),[1 2]) |isequal(size(val1),[2 1]))
                    errordlg({'X axis limit is incorrect';...
                            'currect xample: [3.4,20],[1 30.5], 1 30.5.'});
                    return
                elseif val1(1) >= val1(2)
                    errordlg({'X axis limit is incorrect',...
                            'xample: [val1,val2] (where val1 < val2 ~=0)'});
                    return     
                elseif val1(1)<0 | val1(2) < 0
                    errordlg({'X axis limit is incorrect';...
                            'values shold be positive'});
                    return     
                end
                delete(gcbf);
                figure(figdata.h);
                datas.str1=str1;
                datas.str2=str2;
                set(plotoptionb,'userdata',datas);
            end
        end
        
        
        %-----------------------------------------------------
        % Callback 'change_xwindowsize' when xwindowsize EditBox
        %-----------------------------------------------------
    elseif strcmp(action,'change_xwindowsize');
        strsize=str2num(get(gcbo,'String'));
        strsize=fix(strsize);
        %-------------------------------
        %Error
        if isempty(strsize) |( strsize < 1 ),
            set(gcbo,'String',num2str(par.xwindowsize));
            return
        elseif ( strsize > par.siglength )
            strsize=par.siglength;
            set(gcbo,'String',num2str(strsize));
        end
        %-------------------------------
        par.xwindowsize=strsize;
        icalab_srplot(par);
        %-----------------------------------------------------
        % Callback 'indicator_move' when indicator push
        %-----------------------------------------------------
    elseif strcmp(action,'indicator_move');
        %-----------------------------
        % callback from pushdutton
        if isstr(page)
            if strcmp(page,'down')
                if par.page > 1,
                    page=par.page-1;
                else;return;
                end
            elseif strcmp(page,'up')
                if par.page < par.allpage,
                    page=par.page+1;
                else;return;
                end
            end
        end
        page=ceil(page);      
        if isequal(par.page,page),return;end;
        %-----------------------------
        set(par.push_pageup,'enable','on')
        set(par.push_pagedw,'enable','on')
        if isequal(page,1)
            set(par.push_pagedw,'enable','off')
        end
        if isequal(page,par.allpage)
            set(par.push_pageup,'enable','off')
        end
        par.page=page;
        vm=getappdata(gcf,'IndicatorData');
        
        set(vm.map,'xdata',page-0.5);
        if isequal(page,par.allpage)
            set(vm.text,'position',[(page-1+(par.siglength/par.xlimsize))/2, 1],'string',[num2str(page) '/' num2str(par.allpage)]);
        else
            set(vm.text,'position',[page-0.5, 1],'string',[num2str(page) '/' num2str(par.allpage)]);
        end
        set(vm.text,'FontWeight','bold')
        set(h,'Userdata',par)
        icalab_srplot([],'slider_moved');
        
        %-----------------------------------------------------
        % Callback 'slider__moved' when mslider push
        %-----------------------------------------------------
    elseif strcmp( action, 'slider_moved' );
        signal=par.signal;
        [sensorsno,T] = size(par.signal); % par.plotno - sensors number, T - signals length
        
        if isequal(par.page,par.allpage)
            xlim=[1 par.siglength-par.xwindowsize*(par.page-1)];
        else
            xlim=[1 par.xlimsize];
        end
        slid = findobj( 'Tag', ['MojSlider_' int2str(par.window_number)] );
        if ~isempty( slid )
            val = get( slid, 'Value' );
            val = round( val );
            set( slid, 'Value', val );			
            sstart = sensorsno-par.plotno-val+1;
        else
            sstart = sensorsno-par.plotno+1;      
        end	
    
        
  % Title of the first plot
    if ~isempty( par.title_first_plot )
        title_disp = par.title_first_plot;
    else
        title_disp = [];
    end
     
        
        %ylim
        sc=getappdata(gcf,'ScalingData');
        if strcmp(sc.type,'auto');
            ylimstring='';
        else
            ylimstring='axis([1 2 sc.ylim1 sc.ylim2]);';
        end
        sc.signalylimdata=[min(min(signal(sstart:sstart:par.plotno-1,:))) max(max(signal(sstart:sstart:par.plotno-1,:)))];
        setappdata(gcf,'ScalingData',sc);
        
        i=1;
        high=(defposi(4)-0.1)/par.plotno;
        %     spacehigh=0.1/(par.plotno-1)+high;
        if par.plotno == 1
            spacehigh = high;
        else
            spacehigh=0.1/(par.plotno-1)+high;
        end
        defhigh=defposi(2)+defposi(4)-high;
        for j=sstart:sstart+par.plotno-1,
            s = subplot('position',[defposi(1) defhigh-spacehigh*(i-1) defposi(3) high]);
            plot( real(signal(j,:)) )
            eval(ylimstring);
            set( s, 'GridLineStyle', ':',...
                'XColor', [0 0 0],...
                'XGrid', 'on',...
                'YColor', [0 0 0],...
                'YGrid', 'on',...
                'FontSize', sFontSize,...
                'XLim', xlim+par.xlimsize*(par.page-1),...
                'Tag', ['ICALAB_loc_subplot_' int2str(i), '_' int2str(par.window_number)] );
 
              % Threshold for small signals   
%             if max( abs(max(signal(j,:))), abs(min(signal(j,:))) ) < ICALAB_PlotTreshold,
%                 set( s, 'YLim', [-1 1] );
%             end   
            
            % Title of the first plot
            if (j == 1) & ~isempty(title_disp), tit_first = title(title_disp); set(tit_first,'FontSize',10); end
                
            %-----------------------------------------------------
            % Axis label
            if ~isempty(par.x_limit)
                if i==1
                    xlimstart=1;
                    temp=par.xlimsize/(par.siglength-1)*(par.x_limit(2)-par.x_limit(1))/10;
                    nn=temp;
                    kk=0;
                    if nn>=10
                        while nn>1
                            nn=nn/10;
                            kk=kk+1;
                        end
                    else
                        while nn<=1
                            nn=nn*10;
                            kk=kk-1;
                        end
                    end
                    nn=ceil(nn);
                    tick=nn*10^kk;
                    xticklabel=[0:tick:par.x_limit(2)];
                    xticklabel=xticklabel(find(xticklabel>=par.x_limit(1)));
                    xtick=(xticklabel-par.x_limit(1))*(par.siglength-1)/(par.x_limit(2)-par.x_limit(1))+1;
                    xticklabel=num2str(xticklabel');
                end
                set(s,'Xtick',xtick);
                set(s,'XtickLabel',xticklabel);
            end
            if i == par.plotno
                xlabel(par.x_title);
            else
                set( s, 'XTickLabel', [] );
            end
            %--------------------------------------------------------
            
            yl = ylabel( char(par.labels{j}) );
            set( yl, 'FontWeight', 'bold' );

            set( yl, 'Rotation', 0,...
                'FontSize', sFontSize+1 );
            %			set( yl, 'FontWeight', 'bold' );
            
            i = i + 1;
        end
    end
end
