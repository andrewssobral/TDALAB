function [s_new, ch_num] = selectchn( s, NoOfSubplots, action )
% function [s_new, fs_new, ch_num] = selectchn( s, NoOfSubplots, action )
% Select Channels Figure
%
% written by T.Watanabe and many others
%
TDALAB_select_sources_window_name = 'Select channels analyzing window';
%initialize
s_new =[];
fs_new =[];
ch_num =[];


global backGroundColor
global colorMap
MaxI=size(s,2);

persistent channels_to_load
persistent RealNoOfSubplots
persistent n
persistent T
persistent min_zakr
persistent max_zakr
persistent x_min_display
persistent x_max_display
persistent start_select
persistent stop_select
persistent s_snew
persistent fs_snew
persistent local_MaxI
persistent min_zakr_ax
persistent max_zakr_ax
persistent l_min_zakr
persistent l_max_zakr
persistent all_min_zakr
persistent all_max_zakr

small_font_size_7=7;

h0=[];

global ScreenWidth;
SCR_800x600 = 70; 
SCR_1024x748 = 83;
SCR_1280x1024 = 100;

if ScreenWidth < 900,
    dd = 100/SCR_800x600;
elseif ScreenWidth < 1100,
    dd = 100/SCR_1024x748;
else
    dd = 100/SCR_1280x1024;
end
% some colors

NUM_RANDOM_POINTS = 2000;

[NumberOfColors, p_col] = size(colorMap);

zoom_option = 0;
%%===============================================================
%%=======================


if isempty(s),
    axes_all = findobj( 'Tag', 'TDALAB_select_channels_AxesAll' );
    if ~isempty(axes_all),
        s = get( axes_all, 'UserData' );
    end
    
    main_fig = findobj( 'Tag', 'TDALAB_select_channel_window' );
end

% -############################################################################
if nargin > 2,
    % callback functions
    
    % --------------------------------------------------------------------------
    % OK button
    if strcmp( action, 'button_ok' ) == 1,
        
        is_some_ch = sum(sum(channels_to_load));
        
        if is_some_ch == 0,
%             errordlg( 'There are no selected channels', 'There are no selected channels', 'modal' );
            choice = questdlg('Delete this mode?', 'Confirm', 'Confirm','Cancel','Cancel');
            switch choice
                case 'Confirm'
                    s_snew='delete';
                    close(main_fig);                    
                case 'Cancel'
                    return;
            end
    
        else
            num_ch = [];
            for i=1:n,
                if channels_to_load(i,1) == 1,
                    num_ch = [num_ch i];
                end
            end
            
            s_snew = s( num_ch, start_select:stop_select );
            fs_snew = s_snew;
            ch_num = channels_to_load;
            
            close( main_fig );
        end
        
        
        % --------------------------------------------------------------------------
        % shuffle button 
    elseif strcmp( action, 'button_shuffle' ) == 1,
        
        is_some_ch = sum(sum(channels_to_load));
        
        if is_some_ch == 0,
            errordlg( 'There are no selected channels', 'There are no selected channels', 'modal' );
        else
            num_ch = [];
            for i=1:n,
                if channels_to_load(i,1) == 1,
                    num_ch = [num_ch i];
                end
            end
            
            fs_snew = s( num_ch, start_select:stop_select );
            
            for i=1:length(num_ch),
                s_rand(i,:) = fs_snew(i,randperm(length(fs_snew)));
                
                if NUM_RANDOM_POINTS > length(s_rand(i,:))
                    NUM_RANDOM_POINTS = round(length(s_rand(i,:))*0.9);
                end
                
                s_snew(i,1:NUM_RANDOM_POINTS) = s_rand(i,1:NUM_RANDOM_POINTS);   
            end
            
            close( main_fig );
        end
        
        
        % --------------------------------------------------------------------------
        % Cancel button
    elseif strcmp( action, 'button_cancel' ) == 1,
        
        s_snew = [];
        
%         figg = findobj('tag','s_TDALAB_HEAD_PLOT');
%         if figg,
%             delete(figg);
%         end
        
        close( main_fig );
        
        
        % --------------------------------------------------------------------------
        % Save button
    elseif strcmp( action, 'button_save' ) == 1,
        
        is_some_ch = sum(sum(channels_to_load));
        
        if is_some_ch == 0,
            errordlg( 'There is nothing to save', 'There is nothing to save', 'modal' );
        else
            num_ch = [];
            for i=1:n,
                if channels_to_load(i,1) == 1,
                    num_ch = [num_ch i];
                end
            end
            
            selected_signals = s( num_ch, start_select:stop_select );
            
            [sfile spath] = uiputfile( 'selected_signals.mat', 'Save as' );
            if sfile == 0,
                return
            end
            
            save( [spath sfile], 'selected_signals' );
        end
        
        
        % --------------------------------------------------------------------------
        % Help button
    elseif strcmp( action, 'button_help' ) == 1,
        TDALAB_help('select_channels');
        
        %      figure( main_fig );        
        
        % --------------------------------------------------------------------------
        % Horizontal move
        % --------------------------------------------------------------------------
    elseif strcmp( action, 'horizontal_slider' ) == 1,
        
        local_MaxI = x_max_display - x_min_display;
        
        mslider = findobj( 'Tag', 'SelectChnHorizSlider' );
        val = get( mslider, 'Value' );
        val = round( val );
        set( mslider, 'Value', val );			
        
        local_MaxI = x_max_display - x_min_display;
        x_min_display = val;
        x_max_display = val + local_MaxI;
        
        edit_start = findobj( 'Tag', 'DisplayRangeStartEdit' );
        edit_stop  = findobj( 'Tag', 'DisplayRangeStopEdit'  );
        
        set( edit_start, 'Value', x_min_display );
        set( edit_start, 'String', num2str(x_min_display) );
        set( edit_stop, 'Value', x_max_display );
        set( edit_stop, 'String', num2str(x_max_display) );
        drawnow
        
        axesy = findobj( 'Tag', 'ChannelAxes' );
        set( axesy, 'XLim', [x_min_display x_max_display] ); 
        axesy = findobj( 'Tag', 'TDALAB_select_channels_AxesAll' );
        set( axesy, 'XLim', [x_min_display x_max_display] ); 
        
        l1 = findobj( 'Tag', 'StartDisplayLine' );
        l2 = findobj( 'Tag', 'StopDisplayLine' );
        
        set( l1, 'XData', [x_min_display x_min_display] );
        set( l1, 'YData', [l_min_zakr l_max_zakr] );
        
        set( l2, 'XData', [x_max_display x_max_display] );
        set( l2, 'YData', [l_min_zakr l_max_zakr] );
        
        
        % --------------------------------------------------------------------------
        % set start range to display
    elseif strcmp( action, 'set_start_display_range' ) == 1,
        
        new_val_str = get( gcbo, 'String' );
        new_val = str2double( new_val_str );
        if isnan( new_val )
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            fprintf( '\a' );
        else
            set( gcbo, 'Value', new_val );
            
            if new_val > T,
                new_val = T;
            end
            if new_val < 1,
                new_val = 1;
            end
            
            if new_val < x_max_display,
                x_min_display = new_val;      
            elseif new_val > x_max_display,
                x_min_display = x_max_display;
                x_max_display = new_val;      
            else
                x_min_display = new_val;
                nv = x_min_display + local_MaxI;
                if nv > T
                    nv = T;
                end
                x_max_display = nv;
            end
            
            if x_min_display == 1 & x_max_display == 1,
                nv = x_min_display + local_MaxI;
                if nv > T
                    nv = T;
                end
                x_max_display = nv;
            end
            
            if x_min_display == T & x_min_display == T,
                nv = x_min_display - local_MaxI;
                if nv < 1
                    nv = 1;
                end
                x_min_display = nv;
            end
            
            edit_start = findobj( 'Tag', 'DisplayRangeStartEdit' );
            edit_stop  = findobj( 'Tag', 'DisplayRangeStopEdit'  );
            
            set( edit_start, 'String', num2str(x_min_display) );
            set( edit_stop, 'String', num2str(x_max_display) );
            
            set( edit_start, 'Value', x_min_display );
            set( edit_stop, 'Value', x_max_display );
            
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            
            local_MaxI = x_max_display - x_min_display;
            
            mslider = findobj( 'Tag', 'SelectChnHorizSlider' );
            set( mslider, 'Value', x_min_display );			         
            if x_max_display == T,
                set( mslider, 'Visible', 'off' );
            else            
                set( mslider, 'Visible', 'on' );
                set( mslider, 'Max', T-local_MaxI );
                set( mslider, 'SliderStep', [local_MaxI/(T-local_MaxI)/100 local_MaxI/(T-local_MaxI)] );
            end
            
            axesy = findobj( 'Tag', 'ChannelAxes' );
            set( axesy, 'XLim', [x_min_display x_max_display] ); 
            gaxes = findobj( 'Tag', 'TDALAB_select_channels_AxesAll' );
            set( gaxes, 'XLim', [x_min_display x_max_display] ); 
        end
        
        l1 = findobj( 'Tag', 'StartDisplayLine' );
        l2 = findobj( 'Tag', 'StopDisplayLine' );
        
        set( l1, 'XData', [x_min_display x_min_display] );
        set( l1, 'YData', [l_min_zakr l_max_zakr] );
        
        set( l2, 'XData', [x_max_display x_max_display] );
        set( l2, 'YData', [l_min_zakr l_max_zakr] );
        
        
        % --------------------------------------------------------------------------
        % set stop range to display
    elseif strcmp( action, 'set_stop_display_range' ) == 1,
        
        new_val_str = get( gcbo, 'String' );
        new_val = str2double( new_val_str );
        if isnan( new_val )
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            fprintf( '\a' );
        else
            set( gcbo, 'Value', new_val );
            
            if new_val > T,
                new_val = T;
            end
            if new_val < 1,
                new_val = 1;
            end
            
            if new_val > x_min_display,
                x_max_display = new_val;      
            elseif new_val < x_min_display,
                x_max_display = x_min_display;
                x_min_display = new_val;      
            else
                x_max_display = new_val;
                nv = x_max_display - local_MaxI;
                if nv < 1
                    nv = 1;
                end
                x_min_display = nv;
            end
            
            if x_min_display == 1 & x_max_display == 1,
                nv = x_min_display + local_MaxI;
                if nv > T
                    nv = T;
                end
                x_max_display = nv;
            end
            
            if x_min_display == T & x_min_display == T,
                nv = x_min_display - local_MaxI;
                if nv < 1
                    nv = 1;
                end
                x_min_display = nv;
            end
            
            edit_start = findobj( 'Tag', 'DisplayRangeStartEdit' );
            edit_stop  = findobj( 'Tag', 'DisplayRangeStopEdit'  );
            
            set( edit_start, 'String', num2str(x_min_display) );
            set( edit_stop, 'String', num2str(x_max_display) );
            
            set( edit_start, 'Value', x_min_display );
            set( edit_stop, 'Value', x_max_display );
            
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            
            local_MaxI = x_max_display - x_min_display;
            
            mslider = findobj( 'Tag', 'SelectChnHorizSlider' );
            set( mslider, 'Value', x_min_display );			
            %         set( mslider, 'SliderStep', [local_MaxI/(T)/100 local_MaxI/(T)] );
            
            %         local_MaxI/(T-local_MaxI)/100 
            %         local_MaxI/(T-local_MaxI)
            
            if x_max_display == T,
                set( mslider, 'Visible', 'off' );
            else            
                %            [local_MaxI/(T-local_MaxI)/100 local_MaxI/(T-local_MaxI)]
                set( mslider, 'Visible', 'on' );
                set( mslider, 'Max', T-local_MaxI );
                p1 = local_MaxI/T/100;
                p2 = local_MaxI/T;
                %            if p1> 1
                %               p1 = 1-1/p1;
                %            end
                %            if p2> 1
                %               p2 = 1-1/p2;
                %            end
                if p1 > p2
                    tmp = p2;
                    p2 = p1;
                    p1 = tmp;
                end
                %            p1
                %            p2
                
                %            set( mslider, 'SliderStep', [local_MaxI/(T-local_MaxI)/100 local_MaxI/(T-local_MaxI)] );
                set( mslider, 'SliderStep', [p1 p2] );
            end
            
            axesy = findobj( 'Tag', 'ChannelAxes' );
            set( axesy, 'XLim', [x_min_display x_max_display] ); 
            gaxes = findobj( 'Tag', 'TDALAB_select_channels_AxesAll' );
            set( gaxes, 'XLim', [x_min_display x_max_display] ); 
        end
        
        l1 = findobj( 'Tag', 'StartDisplayLine' );
        l2 = findobj( 'Tag', 'StopDisplayLine' );
        
        set( l1, 'XData', [x_min_display x_min_display] );
        set( l1, 'YData', [min_zakr max_zakr] );
        
        set( l2, 'XData', [x_max_display x_max_display] );
        set( l2, 'YData', [min_zakr max_zakr] );
        
        
        % --------------------------------------------------------------------------
        % Plots rectangle on selecting part of signals
    elseif strcmp( action, 'plot_rectangle' ) == 1,
        
        mslider = findobj( 'Tag', 'Slider_channels' );
        if ~isempty(mslider)
            val = get( mslider, 'Value' );
            start = n-RealNoOfSubplots-val;
        else
            start = 0;
        end
        
        for i=1:RealNoOfSubplots,
            num_sig = i + start;
            
            l1 = findobj( 'Tag', ['StartLimitLine' int2str(i)] );
            l2 = findobj( 'Tag', ['StopLimitLine' int2str(i)] );
            
            ax_ch = get( l1, 'Parent' );
            zakr_y = get( ax_ch, 'YLim' );
            
            set( l1, 'XData', [start_select start_select] );
            %         set( l1, 'YData', [min_zakr_ax(num_sig) max_zakr_ax(num_sig)] );
            set( l1, 'YData', zakr_y );
            
            set( l2, 'XData', [stop_select stop_select] );
            %         set( l2, 'YData', [min_zakr_ax(num_sig) max_zakr_ax(num_sig)] );
            set( l2, 'YData', zakr_y );
            
            if channels_to_load(num_sig,1) == 1,
                set( l1, 'Visible', 'on' );
                set( l2, 'Visible', 'on' );
            else
                set( l1, 'Visible', 'off' );
                set( l2, 'Visible', 'off' );
            end   
        end
        
        al1 = findobj( 'Tag', 'ALLStartLimitLine' );
        al2 = findobj( 'Tag', 'ALLStopLimitLine' );
        bl1 = findobj( 'Tag', 'StartDisplayLine' );
        bl2 = findobj( 'Tag', 'StopDisplayLine' );
        
        all_ax = findobj( 'Tag', 'TDALAB_select_channels_AxesAll' );
        
        chall = findobj( 'Tag', 'Check_Plot_All' );
        act = get( chall, 'Value' );
        
        set( al1, 'XData', [start_select start_select] );
        set( al1, 'YData', [all_min_zakr all_max_zakr] );
        
        set( al2, 'XData', [stop_select stop_select] );
        set( al2, 'YData', [all_min_zakr all_max_zakr] );
        
        set( bl1, 'XData', [x_min_display x_min_display] );
        set( bl1, 'YData', [l_min_zakr l_max_zakr] );
        
        set( bl2, 'XData', [x_max_display x_max_display] );
        set( bl2, 'YData', [l_min_zakr l_max_zakr] );
        
        if act == 1,
            set( al1, 'Visible', 'on' );
            set( al2, 'Visible', 'on' );
            set( bl1, 'Visible', 'on' );
            set( bl2, 'Visible', 'on' );
        else
            set( al1, 'Visible', 'off' );
            set( al2, 'Visible', 'off' );
            set( bl1, 'Visible', 'off' );
            set( bl2, 'Visible', 'off' );
        end   
        
        if zoom_option == 1,
            set( all_ax, 'XLim', [start_select stop_select] );
            set( all_ax, 'YLim', [all_min_zakr all_max_zakr] );
        end
        set( all_ax, 'YLim', [all_min_zakr all_max_zakr] );
        
        
        % --------------------------------------------------------------------------
        % vertical slider
    elseif strcmp( action, 'slider_show_channels' ) == 1,
        
        mslider = findobj( 'Tag', 'Slider_channels' );
        if ~isempty( mslider )
            val = get( mslider, 'Value' );
            val = round( val );
            set( mslider, 'Value', val );			
        end	
        
        start = n-RealNoOfSubplots-val+1;
        for i=1:RealNoOfSubplots,
            idx_col = rem(i+start-1,NumberOfColors);
            if idx_col == 0,
                idx_col = NumberOfColors;
            end
            
            axes_line = findobj( 'Tag', ['AxesLine' int2str(i)] );
            set( axes_line, 'YData', s(i+start-1,:) );
            set( axes_line, 'Color', colorMap(idx_col,:) );
            
            ax_ch = get( axes_line, 'Parent' );
            set(ax_ch, 'YLim', [min_zakr_ax(i+start-1) max_zakr_ax(i+start-1)] );
            
            axes_check = findobj( 'Tag', ['Checkbox' int2str(i)] );
            set( axes_check, 'Value', channels_to_load(i+start-1,1) );
            

            str_tmp = ['s' int2str(i+start-1)];

            
            set( axes_check, 'String', str_tmp );
        end                   
        
        selectchn( [], [], 'plot_rectangle' );
        
        
        % --------------------------------------------------------------------------
        % check/uncheck of signals
    elseif strcmp( action, 'check_signal' ) == 1,
        
        mslider = findobj( 'Tag', 'Slider_channels' );
        if ~isempty(mslider)
            val = get( mslider, 'Value' );
            start = n-RealNoOfSubplots-val;
        else
            start = 0;
        end
        
        num_sig = get( gcbo, 'UserData' ) + start;
        
        curr_val = get( gcbo, 'Value' );
        if curr_val == 1,
            channels_to_load( num_sig, 1 ) = 1;
        else
            channels_to_load( num_sig, 1 ) = 0;
        end
        
        min_zakr =  1e200;
        max_zakr = -1e200;
        
        for i=1:n,
            if channels_to_load(i,1) == 1,
                lmaxs = max( s(i,:) );
                lmins = min( s(i,:) );
                if lmaxs > max_zakr,
                    max_zakr = lmaxs;
                end
                if lmins < min_zakr,
                    min_zakr = lmins;
                end
            end
        end
        
        % zapas 5%
        min_zakr = min_zakr - 0.05*abs(min_zakr);
        max_zakr = max_zakr + 0.05*abs(max_zakr);
        
        all_min_zakr = min_zakr;
        all_max_zakr = max_zakr;
        
        if all_min_zakr == all_max_zakr,
            all_max_zakr = all_min_zakr + 1;
        end
        
        selectchn( [], [], 'check_plot_all_signals' );
        
        
        % --------------------------------------------------------------------------
        % plot/unplot all signals
    elseif strcmp( action, 'check_plot_all_signals' ) == 1,
        
        chall = findobj( 'Tag', 'Check_Plot_All' );
        act = get( chall, 'Value' );
        
        for i=1:n,
            linia = findobj( 'Tag', ['AxesAllLine' int2str(i)] );
            
            if act == 1,
                if channels_to_load(i,1) == 1,
                    set( linia, 'Visible', 'on' );           
                else
                    set( linia, 'Visible', 'off' );           
                end              
            else   
                set( linia, 'Visible', 'off' );
            end
        end
        
        selectchn( [], [], 'plot_rectangle' );
        
        
        % --------------------------------------------------------------------------
        % guzik select all
    elseif strcmp( action, 'select_all_channels' ) == 1,
        
        for i=1:n,
            chch = findobj( 'Tag', ['Checkbox' int2str(i)] );
            if ~isempty(chch),
                set( chch, 'Value', 1 );
            end
            channels_to_load(i,1) = 1;
        end
        
        selectchn( [], [], 'check_plot_all_signals' );
        
        % --------------------------------------------------------------------------
        % guzik unselect all
    elseif strcmp( action, 'unselect_all_channels' ) == 1,
        
        for i=1:n,
            chch = findobj( 'Tag', ['Checkbox' int2str(i)] );
            if ~isempty(chch),
                set( chch, 'Value', 0 );
            end
            channels_to_load(i,1) = 0;
        end
        selectchn( [], [], 'check_plot_all_signals' );
        
        
        % --------------------------------------------------------------------------
        % radio set scale - auto/fixed
    elseif strcmp( action, 'radio_set_scale' ) == 1,
        
        type_scale = get( gcbo, 'UserData' );
        
        radio_auto_scale = findobj( 'Tag', 'Radiobutton_automatic_scale' );
        radio_fixed_scale = findobj( 'Tag', 'Radiobutton_fixed_scale' );
        edit_min_zakr = findobj( 'Tag', 'Edit_set_min_scale' );
        edit_max_zakr = findobj( 'Tag', 'Edit_set_max_scale' );
        
        mslider = findobj( 'Tag', 'Slider_channels' );
        if ~isempty( mslider )
            val = get( mslider, 'Value' );
            val = round( val );
            set( mslider, 'Value', val );			
        else
            val = 0;
        end	
        
        if type_scale == 1, 
            % auto
            set( radio_auto_scale, 'Value', 1 );
            set( radio_fixed_scale, 'Value', 0 );
            set( edit_max_zakr, 'Enable', 'off' );
            set( edit_min_zakr, 'Enable', 'off' );
                        
            set( edit_min_zakr, 'String', 'N/A' );
            set( edit_max_zakr, 'String', 'N/A' );
            
            start = n-RealNoOfSubplots-val;
            for i=1:RealNoOfSubplots,
                axes_line1 = findobj( 'Tag', ['StartLimitLine' int2str(i)] );
                axes_line2 = findobj( 'Tag', ['StopLimitLine' int2str(i)] );
                ax_ch = get( axes_line1, 'Parent' );
                set(ax_ch, 'YLim', [min_zakr_ax(i+start) max_zakr_ax(i+start)] );
                set(axes_line1, 'YData', [min_zakr_ax(i+start) max_zakr_ax(i+start)] );
                set(axes_line2, 'YData', [min_zakr_ax(i+start) max_zakr_ax(i+start)] );
            end
            
            gaxes = findobj( 'Tag', 'TDALAB_select_channels_AxesAll' );
            set( gaxes, 'YLim', [all_min_zakr all_max_zakr] ); 
        else
            % fixed
            set( radio_auto_scale, 'Value', 0 );
            set( radio_fixed_scale, 'Value', 1 );
            
            axesy = findobj( 'Tag', 'ChannelAxes' );
            set( axesy, 'YLim', [min_zakr max_zakr] ); 
            
            set( edit_min_zakr, 'String', num2str(min_zakr) );
            set( edit_max_zakr, 'String', num2str(max_zakr) );
            set( edit_max_zakr, 'Enable', 'on' );
            set( edit_min_zakr, 'Enable', 'on' );
            
            start = n-RealNoOfSubplots-val;
            for i=1:RealNoOfSubplots,
                axes_line1 = findobj( 'Tag', ['StartLimitLine' int2str(i)] );
                axes_line2 = findobj( 'Tag', ['StopLimitLine' int2str(i)] );
                ax_ch = get( axes_line1, 'Parent' );
                set(ax_ch, 'YLim', [min_zakr max_zakr] );
                set(axes_line1, 'YData', [min_zakr max_zakr] );
                set(axes_line2, 'YData', [min_zakr max_zakr] );
            end
        end
        
        
        % --------------------------------------------------------------------------
        % edit set min scale 
    elseif strcmp( action, 'edit_set_min_scale' ) == 1,
        
        new_val_str = get( gcbo, 'String' );
        new_val = str2double( new_val_str );
        if isnan( new_val )
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            fprintf( '\a' );
        else
            set( gcbo, 'Value', new_val );
            min_zakr = new_val;
            
            if min_zakr > max_zakr,
                tmp = min_zakr;
                min_zakr = max_zakr;
                max_zakr = tmp;
                
                edit_min_zakr = findobj( 'Tag', 'Edit_set_min_scale' );
                edit_max_zakr = findobj( 'Tag', 'Edit_set_max_scale' );
                
                set( edit_min_zakr, 'String', num2str(min_zakr) );
                set( edit_max_zakr, 'String', num2str(max_zakr) );
                
                set( edit_min_zakr, 'Value', min_zakr );
                set( edit_max_zakr, 'Value', max_zakr );
            end
            
            if min_zakr == max_zakr,
                max_zakr = max(max(s));
                min_zakr = min(min(s));
                
                % zapas 5%
                min_zakr = min_zakr - 0.05*abs(min_zakr);
                max_zakr = max_zakr + 0.05*abs(max_zakr);
                
                edit_min_zakr = findobj( 'Tag', 'Edit_set_min_scale' );
                edit_max_zakr = findobj( 'Tag', 'Edit_set_max_scale' );
                
                set( edit_min_zakr, 'String', num2str(min_zakr) );
                set( edit_max_zakr, 'String', num2str(max_zakr) );
                
                set( edit_min_zakr, 'Value', min_zakr );
                set( edit_max_zakr, 'Value', max_zakr );
            end
            
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            
            axesy = findobj( 'Tag', 'ChannelAxes' );
            set( axesy, 'YLim', [min_zakr max_zakr] ); 
            
            %         gaxes = findobj( 'Tag', 'TDALAB_select_channels_AxesAll' );
            %         set( gaxes, 'YLim', [min_zakr max_zakr] ); 
        end
        
        selectchn( [], [], 'plot_rectangle' );
        
        
        % --------------------------------------------------------------------------
        % edit set max scale 
    elseif strcmp( action, 'edit_set_max_scale' ) == 1,
        
        new_val_str = get( gcbo, 'String' );
        new_val = str2double( new_val_str );
        if isnan( new_val )
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            fprintf( '\a' );
        else
            set( gcbo, 'Value', new_val );
            max_zakr = new_val;
            
            if min_zakr > max_zakr,
                tmp = min_zakr;
                min_zakr = max_zakr;
                max_zakr = tmp;
                
                edit_min_zakr = findobj( 'Tag', 'Edit_set_min_scale' );
                edit_max_zakr = findobj( 'Tag', 'Edit_set_max_scale' );
                
                set( edit_min_zakr, 'String', num2str(min_zakr) );
                set( edit_max_zakr, 'String', num2str(max_zakr) );
                
                set( edit_min_zakr, 'Value', min_zakr );
                set( edit_max_zakr, 'Value', max_zakr );
            end
            
            if min_zakr == max_zakr,
                max_zakr = max(max(s));
                min_zakr = min(min(s));
                
                % zapas 5%
                min_zakr = min_zakr - 0.05*abs(min_zakr);
                max_zakr = max_zakr + 0.05*abs(max_zakr);
                
                edit_min_zakr = findobj( 'Tag', 'Edit_set_min_scale' );
                edit_max_zakr = findobj( 'Tag', 'Edit_set_max_scale' );
                
                set( edit_min_zakr, 'String', num2str(min_zakr) );
                set( edit_max_zakr, 'String', num2str(max_zakr) );
                
                set( edit_min_zakr, 'Value', min_zakr );
                set( edit_max_zakr, 'Value', max_zakr );
            end
            
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            
            axesy = findobj( 'Tag', 'ChannelAxes' );
            set( axesy, 'YLim', [min_zakr max_zakr] ); 
            
            %         gaxes = findobj( 'Tag', 'TDALAB_select_channels_AxesAll' );
            %         set( gaxes, 'YLim', [min_zakr max_zakr] ); 
        end
        
        selectchn( [], [], 'plot_rectangle' );
        
        
        %%   % --------------------------------------------------------------------------
        %%   % slider set min zakr 
        %%   elseif strcmp( action, 'slider_set_start_pos' ) == 1,
        %%      
        %%      val = get( gcbo, 'Value' );
        %%      val = round(val);
        %%         
        %%      if val < stop_select,
        %%         start_select = val;      
        %%      elseif val > stop_select,
        %%         start_select = stop_select;
        %%         stop_select = val;      
        %%      else
        %%         start_select = val;
        %%         nv = start_select + round(T/10);
        %%         if nv > T
        %%            nv = T;
        %%         end
        %%         stop_select = nv;
        %%      end
        %%      
        %%      if start_select == 1 & stop_select == 1,
        %%         nv = start_select + round(T/10);
        %%         if nv > T
        %%            nv = T;
        %%         end
        %%         stop_select = nv;
        %%      end
        %%      
        %%      if start_select == T & stop_select == T,
        %%         nv = stop_select - round(T/10);
        %%         if nv < 1
        %%            nv = 1;
        %%         end
        %%         start_select = nv;
        %%      end
        %%      
        %%      sl_start = findobj( 'Tag','Slider_Set_Start' );
        %%      sl_stop = findobj( 'Tag','Slider_Set_Stop' );
        %%      set( sl_start, 'Value', start_select );
        %%      set( sl_stop, 'Value', stop_select );
        %%
        %%      ed_start = findobj( 'Tag', 'Edit_Set_Start' );
        %%      ed_stop  = findobj( 'Tag', 'Edit_Set_Stop' );
        %%      set( ed_start, 'Value', start_select );
        %%      set( ed_stop,  'Value', stop_select );
        %%      set( ed_start, 'String', num2str(start_select) );
        %%      set( ed_stop,  'String', num2str(stop_select) );
        %%      
        %%      selectchn( [], [], 'plot_rectangle' );
        
        
        %%   % --------------------------------------------------------------------------
        %%   % slider set max zakr 
        %%   elseif strcmp( action, 'slider_set_stop_pos' ) == 1,
        %%      
        %%      val = get( gcbo, 'Value' );
        %%      val = round(val);
        %%      set( gcbo, 'Value', val );
        %%      
        %%      if val > start_select,
        %%         stop_select = val;      
        %%      elseif val < start_select,
        %%         stop_select = start_select;
        %%         start_select = val;      
        %%      else
        %%         stop_select = val;
        %%         nv = start_select - round(T/10);
        %%         if nv < 1
        %%            nv = 1;
        %%         end
        %%         start_select = nv;
        %%      end
        %%      
        %%      if start_select == 1 & stop_select == 1,
        %%         nv = start_select + round(T/10);
        %%         if nv > T
        %%            nv = T;
        %%         end
        %%         stop_select = nv;
        %%      end
        %%      
        %%      if start_select == T & stop_select == T,
        %%         nv = stop_select - round(T/10);
        %%         if nv < 1
        %%            nv = 1;
        %%         end
        %%         start_select = nv;
        %%      end
        %%      
        %%      sl_start = findobj( 'Tag','Slider_Set_Start' );
        %%      sl_stop = findobj( 'Tag','Slider_Set_Stop' );
        %%      set( sl_start, 'Value', start_select );
        %%      set( sl_stop, 'Value', stop_select );
        %%      
        %%      ed_start = findobj( 'Tag', 'Edit_Set_Start' );
        %%      ed_stop  = findobj( 'Tag', 'Edit_Set_Stop' );
        %%      set( ed_start, 'Value', start_select );
        %%      set( ed_stop,  'Value', stop_select );
        %%      set( ed_start, 'String', num2str(start_select) );
        %%      set( ed_stop,  'String', num2str(stop_select) );
        %%      
        %%      selectchn( [], [], 'plot_rectangle' );
        
        % --------------------------------------------------------------------------
        % edit set min zakr 
    elseif strcmp( action, 'edit_set_start_pos' ) == 1,
        
        ed_start = findobj( 'Tag', 'Edit_Set_Start' );
        ed_stop  = findobj( 'Tag', 'Edit_Set_Stop' );
        %%      sl_start = findobj( 'Tag','Slider_Set_Start' );
        %%      sl_stop = findobj( 'Tag','Slider_Set_Stop' );
        
        val_str = get( gcbo, 'String' );
        val = str2double( val_str );
        if isnan( val )
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            fprintf( '\a' );
        else
            if val > T,
                val = T;
            end
            if val < 1,
                val = 1;
            end
            
            if val < stop_select,
                start_select = val;      
            elseif val > stop_select,
                start_select = stop_select;
                stop_select = val;      
            else
                start_select = val;
                nv = start_select + round(T/10);
                if nv > T
                    nv = T;
                end
                stop_select = nv;
            end
            
            if start_select == 1 & stop_select == 1,
                nv = start_select + round(T/10);
                if nv > T
                    nv = T;
                end
                stop_select = nv;
            end
            
            if start_select == T & stop_select == T,
                nv = stop_select - round(T/10);
                if nv < 1
                    nv = 1;
                end
                start_select = nv;
            end
            
            %%         set( sl_start, 'Value', start_select );
            %%         set( sl_stop,  'Value', stop_select );
            set( ed_start, 'Value', start_select );
            set( ed_stop,  'Value', stop_select );
            set( ed_start, 'String', num2str(start_select) );
            set( ed_stop,  'String', num2str(stop_select) );
        end
        
        selectchn( [], [], 'plot_rectangle' );
        
        
        % --------------------------------------------------------------------------
        % edit set max zakr 
    elseif strcmp( action, 'edit_set_stop_pos' ) == 1,
        
        ed_start = findobj( 'Tag', 'Edit_Set_Start' );
        ed_stop  = findobj( 'Tag', 'Edit_Set_Stop' );
        %%      sl_start = findobj( 'Tag','Slider_Set_Start' );
        %%      sl_stop = findobj( 'Tag','Slider_Set_Stop' );
        
        val_str = get( gcbo, 'String' );
        val = str2double( val_str );
        if isnan( val )
            set( gcbo, 'String', num2str(get(gcbo,'Value')) );
            fprintf( '\a' );
        else
            if val > T,
                val = T;
            end
            if val < 1,
                val = 1;
            end
            
            if val > start_select,
                stop_select = val;      
            elseif val < start_select,
                stop_select = start_select;
                start_select = val;      
            else
                stop_select = val;
                nv = start_select - round(T/10);
                if nv < 1
                    nv = 1;
                end
                start_select = nv;
            end
            
            if start_select == 1 & stop_select == 1,
                nv = start_select + round(T/10);
                if nv > T
                    nv = T;
                end
                stop_select = nv;
            end
            
            if start_select == T & stop_select == T,
                nv = stop_select - round(T/10);
                if nv < 1
                    nv = 1;
                end
                start_select = nv;
            end
            
            %%         set( sl_start, 'Value', start_select );
            %%         set( sl_stop,  'Value', stop_select );
            set( ed_start, 'Value', start_select );
            set( ed_stop,  'Value', stop_select );
            set( ed_start, 'String', num2str(start_select) );
            set( ed_stop,  'String', num2str(stop_select) );
        end
        
        selectchn( [], [], 'plot_rectangle' );
        
        
        % --------------------------------------------------------------------------
   
        
    else
        % nieznane zdarzenie ?!
        fprintf( '\n\nTDALAB: Internal error. Please restart program.\n\n' );
        
        
        
    end
    
    
    
    
    % -############################################################################
    % -----------------------------------------------------------------------------
    % Creat Figure
    % -----------------------------------------------------------------------------
else
    if nargin == 1
        NoOfSubplots = 10;
    end
    
    [n,T] = size( s ); % n - sensors number, T - signals length
    RealNoOfSubplots = min( n, NoOfSubplots );
    
    ch_snum = ones( n, 1 );
    channels_to_load = ones( n, 1 );
    
    min_zakr = min(min(s));
    max_zakr = max(max(s));
    
    % zapas 5%
    min_zakr = min_zakr - 0.05*abs(min_zakr);
    max_zakr = max_zakr + 0.05*abs(max_zakr);
    
    l_min_zakr = min_zakr;
    l_max_zakr = max_zakr;
    
    all_min_zakr = min_zakr;
    all_max_zakr = max_zakr;
    
    if l_min_zakr == l_max_zakr,
        l_max_zakr = l_min_zakr + 1;
    end
    
    if all_min_zakr == all_max_zakr,
        all_max_zakr = all_min_zakr + 1;
    end
    
    min_zakr_ax = min( s' );
    max_zakr_ax = max( s' );
    temp=max_zakr_ax-min_zakr_ax;
    meantemp = mean(temp);
    tempidx= find(temp==0);
    if ~isempty(tempidx)
        min_zakr_ax(tempidx)=-meantemp;
        max_zakr_ax(tempidx)=meantemp;
    end
    
    
    % zapas 5%
    min_zakr_ax = min_zakr_ax - 0.05*abs(min_zakr_ax);
    max_zakr_ax = max_zakr_ax + 0.05*abs(max_zakr_ax);
    
    if T < MaxI
        x_min_display = 1;
        x_max_display = T;
        local_MaxI = T;
    else
        x_min_display = 1;
        x_max_display = MaxI;
        local_MaxI = MaxI;
    end
    
    start_select = 1;
    stop_select = T;
    
    s_snew = [];
    
    main_fig = findobj( 'Tag', 'TDALAB_select_channel_window' );
    if ~isempty( main_fig )
        close(main_fig);
    end
    
    h0 = figure(	'Units','characters', ...
        'Color', backGroundColor, ...
        'Position', [1/dd 1/dd 160/dd 66/dd], ...
        'Visible', 'off', ...
        'NumberTitle', 'off', ...
        'Resize', 'on', ...
        'Tag', 'TDALAB_select_channel_window', ...
        'Name', 'Modify/Update channels', ...
        'MenuBar', 'none', ...
        'ToolBar', 'none' );
    
    defposi=[0.09 0.18 0.86 0.68];
    
    %---------------------------------------------------
    % create Axis 'All'               
    axes_all = axes( 'Parent', h0, ...
        'DrawMode', 'fast', ...
        'Box', 'on', ...
        'Color', [1 1 1], ...
        'FontSize', small_font_size_7, ...
        'NextPlot','add', ...
        'Position', [defposi(1) defposi(2)+defposi(4)+0.02 defposi(3) 0.05], ...
        'Tag', 'TDALAB_select_channels_AxesAll', ...
        'XLim', [1 T], ...
        'YLim', [min_zakr max_zakr], ...
        'XGrid', 'on', ...
        'YGrid', 'on', ...
        'XColor', [0 0 0], ...
        'YColor', [0 0 0], ...
        'ZColor', [0 0 0], ...
        'UserData', s );
    %---------------------------------------------------
    % Lines of Axis 'All'                 
    for i=1:n,
        idx = rem(i,NumberOfColors);
        if idx == 0, 
            idx = NumberOfColors; 
        end
        
        h2 = line( 'Parent', axes_all, ...
            'Color', colorMap(idx,:), ...
            'Tag', ['AxesAllLine' int2str(i)], ...
            'XData', [1:T], ...
            'YData', s(i,:), ...
            'Visible', 'off' );
    end                       
    l1 = line( 'Parent', axes_all, ...
        'Tag', 'ALLStartLimitLine', ...
        'Color', [ 0 0 0 ], ... 
        'LineWidth', 2, ...
        'XData', [start_select start_select], ...
        'YData', [min_zakr max_zakr], ...
        'Visible', 'off' );
    l2 = line( 'Parent', axes_all, ...
        'Tag', 'ALLStopLimitLine', ...
        'Color', [ 0 0 0 ], ...
        'LineWidth', 2, ...
        'XData', [stop_select stop_select], ...
        'YData', [min_zakr max_zakr], ...
        'Visible', 'off' );
    l1 = line( 'Parent', axes_all, ...
        'Tag', 'StartDisplayLine', ...
        'Color', [ 0.501960784313725 0.501960784313725 0.501960784313725 ], ...
        'LineWidth', 4, ...
        'XData', [x_min_display x_min_display], ...
        'YData', [min_zakr max_zakr], ...
        'Visible', 'off' );
    l2 = line( 'Parent', axes_all, ...
        'Tag', 'StopDisplayLine', ...
        'Color', [ 0.501960784313725 0.501960784313725 0.501960784313725 ], ...
        'LineWidth', 4, ...
        'XData', [x_max_display x_max_display], ...
        'YData', [min_zakr max_zakr], ...
        'Visible', 'off' );
    %---------------------------------------------------
    % Check Box Of channel 'All'
    h1 = uicontrol( 'Parent', h0, ...
        'Units','normalized',...
        'BackgroundColor', backGroundColor, ...
        'Callback', 'selectchn( [], [], ''check_plot_all_signals'' );', ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [0.015 defposi(2)+defposi(4)+0.01+0.05/2 0.05 0.02], ...
        'String', 'All', ...
        'Style', 'checkbox', ...
        'Tag', 'Check_Plot_All', ...
        'TooltipString', 'Check to show all channels on one plot' );
    
    
    % edit do ustalania pozycji startowej sygnalow                
    h1 = uicontrol( 'Parent', h0, ...
        'Units','characters', ...
        'BackgroundColor', [1 1 1], ...
        'ListboxTop',0, ...
        'Callback', 'selectchn( [], [], ''edit_set_start_pos'' );', ...
        'Position',[30/dd 6.7/dd 15/dd 1.53846153846154/dd], ...
        'Style','edit', ...
        'Value', start_select, ...
        'String', num2str(start_select), ...
        'Tag','Edit_Set_Start');
    
    % edit do ustalania pozycji koncowej sygnalow                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', [1 1 1], ...
        'ListboxTop', 0, ...
        'Callback', 'selectchn( [], [], ''edit_set_stop_pos'' );', ...
        'Position', [55/dd 6.7/dd 15/dd 1.53846153846154/dd], ...
        'Style', 'edit', ...
        'Value', stop_select, ...
        'String', num2str(stop_select), ...
        'Tag', 'Edit_Set_Stop' );
    
    % text 'Start'                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [24/dd 6.5/dd 6/dd 1.53846153846154/dd], ...
        'String', 'start', ...
        'Style', 'text', ...
        'Tag', 'StaticText_Start' );
    
    % text 'Stop'                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position',[49/dd 6.5/dd 6/dd 1.53846153846154/dd], ...
        'String','end', ...
        'Style','text', ...
        'Tag','StaticText_Stop');
    
    % text 'Scale'                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [66/dd 4/dd 9/dd 1.53846153846154/dd], ...
        'String', 'Scale:', ...
        'Style','text', ...
        'Tag','StaticText_Scale');
    
    % Radio dla ustawiania auto scale
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'ListboxTop', 0, ...
        'Callback', 'selectchn( [], [], ''radio_set_scale'' );', ...
        'Position', [77/dd 4.2/dd 14/dd 1.53846153846154/dd], ...
        'String', 'automatic', ...
        'Style', 'radiobutton', ...
        'Tag', 'Radiobutton_automatic_scale', ...
        'TooltipString','Check for automatic scale for all plots', ...
        'Value', 1, ...
        'UserData', 1 );
    
    % Radio dla ustawiania fixed scale
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'ListboxTop', 0, ...
        'Callback', 'selectchn( [], [], ''radio_set_scale'' );', ...
        'Position', [94/dd 4.2/dd 14/dd 1.53846153846154/dd], ...
        'String','fixed', ...
        'Style','radiobutton', ...
        'Tag','Radiobutton_fixed_scale', ...
        'TooltipString','Check for fixed scale for all plots', ...
        'Value', 0, ...
        'UserData', 2 );
    
    % napis 'min:'                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [108/dd 4/dd 5/dd 1.53846153846154/dd], ...
        'String', 'min:', ...
        'Style','text', ...
        'Tag','StaticText_min');
    
    % napis 'max:'                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [133/dd 4/dd 5/dd 1.53846153846154/dd], ...
        'String', 'max:', ...
        'Style', 'text', ...
        'Tag','StaticText_max' );
    
    % edit do ustawiania min scale                 
    h1 = uicontrol( 'Parent', h0, ...
        'Units','characters', ...
        'BackgroundColor',[1 1 1], ...
        'ListboxTop',0, ...
        'Callback', 'selectchn( [], [], ''edit_set_min_scale'' );', ...
        'Position',[114/dd 4.2/dd 15/dd 1.53846153846154/dd], ...
        'Style','edit', ...
        'String', 'N/A', ...
        'Value', max_zakr, ...
        'Enable', 'off', ...
        'Tag', 'Edit_set_min_scale' );
    
    % edit do ustawiania min scale                 
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', [1 1 1], ...
        'ListboxTop', 0, ...
        'Callback', 'selectchn( [], [], ''edit_set_max_scale'' );', ...
        'Position', [139/dd 4.2/dd 15/dd 1.53846153846154/dd], ...
        'Style', 'edit', ...
        'String', 'N/A', ...
        'Value', min_zakr, ...
        'Enable', 'off', ...
        'Tag', 'Edit_set_max_scale' );
    
    % button selecting all channels to load                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'ListboxTop', 0, ...
        'Position', [6/dd 4/dd 18/dd 1.8/dd], ...
        'String', 'Select all', ...
        'Callback', 'selectchn( [], [], ''select_all_channels'' );', ...
        'Tag', 'Pushbutton_select_all', ...
        'TooltipString','Press for selecting all plots' );
    
    % button unselecting all channels to load                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'ListboxTop',0, ...
        'Callback', 'selectchn( [], [], ''unselect_all_channels'' );', ...
        'Position', [28/dd 4/dd 18/dd 1.8/dd], ...
        'String', 'Unselect all', ...
        'Tag', 'Pushbutton_unselect_all', ...
        'TooltipString', 'Press for unselecting all plots');
    
    % text 'Load range:'             
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [6/dd 6.5/dd 16/dd 1.53846153846154/dd], ...
        'String', 'Load range:', ...
        'Style', 'text', ...
        'Tag', 'LoadRangeText' );
    
    % text 'Display range:'             
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [88/dd 6.5/dd 18/dd 1.53846153846154/dd], ...
        'String', 'Display range:', ...
        'Style', 'text', ...
        'Tag', 'DisplayRangeText');
    
    % text 'start' - display range                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [108/dd 6.5/dd 5/dd 1.53846153846154/dd], ...
        'String', 'start:', ...
        'Style', 'text', ...
        'Tag', 'DisplayRangeStartText');
    
    % text 'end' - display range                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'HorizontalAlignment', 'left', ...
        'ListboxTop', 0, ...
        'Position', [133/dd 6.5/dd 5/dd 1.53846153846154/dd], ...
        'String', 'end:', ...
        'Style', 'text', ...
        'Tag', 'DisplayRangeStopText');
    
    % edit 'DisplayRangeStartEdit'             
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', [1 1 1], ...
        'ListboxTop', 0, ...
        'Callback', 'selectchn( [], [], ''set_start_display_range'' );', ...
        'Position', [114/dd 6.7/dd 15/dd 1.53846153846154/dd], ...
        'String', '1', ...
        'Value', 1, ...
        'Style', 'edit', ...
        'Tag', 'DisplayRangeStartEdit' );
    
    % edit 'DisplayRangeStopEdit'             
    h1 = uicontrol( 'Parent', h0, ...
        'Units','characters', ...
        'BackgroundColor',[1 1 1], ...
        'ListboxTop',0, ...
        'Callback', 'selectchn( [], [], ''set_stop_display_range'' );', ...
        'Position',[139/dd 6.7/dd 15/dd 1.53846153846154/dd], ...
        'String', num2str(local_MaxI), ...
        'Value', local_MaxI+1, ...
        'Style','edit', ...
        'Tag', 'DisplayRangeStopEdit' );
    
    % horizontal slider                
    if T <= MaxI
        h1 = uicontrol( 'Parent', h0, ...
            'Units', 'Normalized', ...
            'BackgroundColor', backGroundColor, ...
            'Callback', 'selectchn( [], [], ''horizontal_slider'' );', ...
            'ListboxTop',0, ...
            'Position',[defposi(1) defposi(2)-0.04 defposi(3) 0.022],...
            'Style','slider', ...
            'Min', 1, ...
            'Max', T, ...
            'Value', 1, ...
            'Visible', 'off', ...
            'Tag','SelectChnHorizSlider' );
    else
        h1 = uicontrol( 'Parent', h0, ...
            'Units', 'Normalized', ...
            'BackgroundColor', backGroundColor, ...
            'Callback', 'selectchn( [], [], ''horizontal_slider'' );', ...
            'ListboxTop',0, ...
            'Position',[defposi(1) defposi(2)-0.04 defposi(3) 0.022],...
            'Style','slider', ...
            'Min', 1, ...
            'Max', T-local_MaxI+1, ...
            'Value', 1, ...
            'SliderStep', [local_MaxI/(T-local_MaxI)/100 local_MaxI/(T-local_MaxI)], ...
            'Tag','SelectChnHorizSlider' );
    end   
    
    if n > NoOfSubplots, 
        % suwak do przesuwania wykresow                
        mslider = uicontrol( 'Parent', h0, ...
            'Units', 'normalized', ...
            'BackgroundColor', backGroundColor, ...
            'ListboxTop', 0, ...
            'Callback', 'selectchn( [], [], ''slider_show_channels'' );', ...
            'Position',[defposi(1)+defposi(3)+0.02 defposi(2) 0.02 defposi(4)],...
            'Style', 'slider', ...
            'Min', 0, ...
            'Max', n-RealNoOfSubplots, ...
            'Value', n-RealNoOfSubplots, ...
            'Tag', 'Slider_channels' );
        
        if RealNoOfSubplots == 1
            set( mslider, 'SliderStep', [1/(n-RealNoOfSubplots) (1+RealNoOfSubplots)/(n-RealNoOfSubplots)] );
        else
            set( mslider, 'SliderStep', [1/(n-RealNoOfSubplots) RealNoOfSubplots/(n-RealNoOfSubplots)] );
        end
    end
    
    % button 'OK'
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor' ,backGroundColor, ...
        'ListboxTop', 0, ...
        'Callback', 'selectchn( [], [], ''button_ok'' );', ...
        'Position', [10/dd 1/dd 20/dd 2/dd], ...
        'String', 'Load', ...
        'Tag', 'Pushbutton_OK' );
    
   
    
    % button 'Cancel'                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'Callback', 'selectchn( [], [], ''button_cancel'' );', ...
        'ListboxTop', 0, ...
        'Position', [40/dd 1/dd 20/dd 2/dd], ...
        'String', 'Cancel', ...
        'Tag', 'Pushbutton_Cancel' );
    
    % button 'Save'
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'Callback', 'selectchn( [], [], ''button_save'' );', ...
        'ListboxTop', 0, ...
        'Position', [70/dd 1/dd 20/dd 2/dd], ...
        'String', 'Save channels', ...
        'Tag', 'Pushbutton_Save', ...
        'TooltipString', 'Press to save selected channels' );
    
    % button 'analyzing window'

    
    
    % button 'Help'                
    h1 = uicontrol( 'Parent', h0, ...
        'Units', 'characters', ...
        'BackgroundColor', backGroundColor, ...
        'Callback', 'TDALAB_help', ...
        'ListboxTop', 0, ...
        'Callback', 'selectchn( [], [], ''button_help'' );', ...
        'Position', [130/dd 1/dd 20/dd 2/dd], ...
        'String', 'Help', ...
        'Tag', 'Pushbutton_help' );
    
    %-----------------------------------------------------
    % Creat Axes
    
    %dl = 45 / (10 * RealNoOfSubplots);             
    %wy = (45-(RealNoOfSubplots-1)*dl) / RealNoOfSubplots;
    
    high=(defposi(4)-0.07)/RealNoOfSubplots;
    spacehigh=0.07/(RealNoOfSubplots-1)+high;
    defhigh=defposi(2)+defposi(4)-high;
    start = 1;
    for i=start:start+RealNoOfSubplots-1,
        h1 = axes( 'Parent', h0, ...
            'Units','Normalized', ...
            'Box', 'on', ...
            'DrawMode', 'fast', ...
            'Color', [1 1 1], ...
            'FontSize', small_font_size_7, ...
            'Position',[defposi(1) defhigh-spacehigh*(i-1) defposi(3) high],...
            'Tag', 'ChannelAxes', ...
            'XLim', [1 T], ...
            'YLim', [min_zakr_ax(i) max_zakr_ax(i)], ...
            'XGrid', 'on', ...
            'YGrid', 'on', ...
            'XColor', [0 0 0], ...
            'YColor', [0 0 0], ...
            'ZColor', [0 0 0] );
        
        if i ~= start+RealNoOfSubplots-1,
            set( h1, 'XTickLabel', '' );
        end
        
        idx = rem(i,NumberOfColors);
        if idx == 0, 
            idx = NumberOfColors; 
        end
        
        h2 = line( 'Parent', h1, ...
            'Tag', ['AxesLine' int2str(i)], ...
            'Color', colorMap(idx+start-1,:), ...
            'XData', 1:T, ...
            'YData', s(i,:) );
        
        
        set( h1, 'XLim', [x_min_display x_max_display+1] );
        set( h1, 'YLim', [min_zakr_ax(i) max_zakr_ax(i)] );
        

         str_tmp = ['s' int2str(i)];

        %----------------------------------------------------------------
        % Creat Check Boxes
        h3 = uicontrol( 'Parent', h0, ...
            'Units', 'Normalized', ...
            'HorizontalAlignment', 'left', ...
            'ListboxTop', 0, ...
            'Callback', 'selectchn( [], [], ''check_signal'' );', ...
            'BackgroundColor', backGroundColor, ...
            'Position',[0.015 defhigh-spacehigh*(i-1)+0.05/2 0.05 0.02],...
            'String', str_tmp, ...
            'Style', 'checkbox', ...
            'Tag', ['Checkbox' int2str(i)], ...
            'Value', 1, ...
            'UserData', i, ...
            'TooltipString', ['Check to load signal number ' int2str(i)] );
        %----------------------------------------------------------------
        local_zakr = get( h1, 'YLim' );
        l1 = line( 'Parent', h1, ...
            'Tag', ['StartLimitLine' int2str(i)], ...
            'Color', [ 0 0 0 ], ...
            'LineWidth', 2, ...
            'XData', [start_select start_select], ...
            'YData', [local_zakr], ...
            'Visible', 'on' );
        l2 = line( 'Parent', h1, ...
            'Tag', ['StopLimitLine' int2str(i)], ...
            'Color', [ 0 0 0 ], ...
            'LineWidth', 2, ...
            'XData', [stop_select stop_select], ...
            'YData', [local_zakr], ...
            'Visible', 'on' );
    end                   
    set(findobj('tag','TDALAB_select_channels_AxesAll'), 'XLim', [x_min_display x_max_display+1] );
    %-----------------------------------------------------
    h_sliders = findobj( 'Style', 'slider' );
    h_axess = findobj( 'Type', 'axes' );
    h_buttons = findobj( 'Style', 'pushbutton' );
    h_rbuttons = findobj( 'Style', 'radiobutton' );
    h_texts = findobj( 'Style', 'text' );
    h_edits = findobj( 'Style', 'edit' );
    h_checks = findobj( 'Style', 'checkbox' );
    
    set( h0, 'Units', 'normalized' );
    set( h_sliders, 'Units', 'normalized' );
    set( h_axess, 'Units', 'normalized' );
    set( h_buttons, 'Units', 'normalized' );
    set( h_rbuttons, 'Units', 'normalized' );
    set( h_texts, 'Units', 'normalized' );
    set( h_edits, 'Units', 'normalized' );
    set( h_checks, 'Units', 'normalized' );
    
    pos = get( h0, 'Position' );
    
    %   pos(3) = pos(3)/dd;
    %   pos(4) = pos(4)/dd;
    
    if pos(4) >= 1
        pos(4) = 0.9;
    end
    
    if pos(3) >= 1
        pos(3) = 0.9;
    end
    
    set( h0, 'Position', [0.5-pos(3)/2 0.5-pos(4)/2 pos(3) pos(4)] );
    
    set( h0, 'Visible', 'on' );
    
    uiwait( h0 );
    
%     if findobj('tag','TDALAB')
%         
        s_new = s_snew;
        ch_num = channels_to_load;
%     end
end