function tcmgPlotFormat(varargin)  
    % Base version author: Lino Mediavilla
    % Contact: lino.mediavilla@estud.usfq.edu.ec
    % Last 29/09/2017
    %
    % Current Structure:
    % tcmgPlotFormat(tString,xString,yString,zString,style,fontSize,interpreter,legendObj, yRot, showGrid, showColorbar)
    %
    %** Basic Arguments ** 
    % tString:      plot title
    % xString:      x-axis label
    % yString:      y-axis label 
    % zString:      z-axis label (if plot is 2D, fill with '') 
    % fontSize:     font size for the text in labels (text size in legend = 0.8*fontSize)
    % interpreter:  interpreter for the text in labels and legend (latex, default, none, etc.)
    % legendObj:    legend object (first create, then pass to this function)
    % style:        plot background and foreground 'dark' or 'light'
    % yRot:         y-label rotation angle (0 = parallel to x-axis)
    % showGrid:     pass a 1 if you want to show the grid, or 0 otherwise.
    % showColorbar: pass a 1 if you want to show the colorbar, or 0 otherwise.
    % 
    % IMPORTANT: To leave intermediate positions in the arguments list 
    % unspecified just type ''
    %
    % %** 2D Example: **
    % close all;clc
    % timeData = 0:.01:10; angleData = sin(timeData); %Data to plot 
    %
    % %[FULLY SPECIFIED ARGUMENTS LIST]:
    %    figure(1); set(gcf,'name','[FULLY SPECIFIED ARGUMENTS LIST]')
    %    plot(timeData,angleData,'r') %Creation of plot 
    %    myLegend = legend('$\theta$');
    %    tcmgPlotFormat('Angular position', 't', '$\theta(t)$', '',  'dark', 20, 'latex', myLegend, 0, 1)
    %
    % %[DEFAULT INTERPRETER, DEFAULT STYLE, WITH GRID]:
    %    figure(2); set(gcf,'name','[DEFAULT INTERPRETER, DEFAULT STYLE, WITH GRID]');
    %    plot(timeData,angleData,'r') %Creation of plot 
    %    tcmgPlotFormat('', 't', '\theta(t)', '', '', 20, 'default','','',1)  
    %
    % %[JUST SPECIFYING TITLE, X-Y LABELS AND STYLE, NO INTERPRETER] 
    %    figure(3); set(gcf,'name', '[JUST SPECIFYING TITLE, X-Y LABELS AND STYLE, NO INTERPRETER]');
    %    plot(timeData,angleData,'r') %Creation of plot 
    %    tcmgPlotFormat('The Title', 't', 'theta(t)', '', 'dark')
    % 
    % %** 3D Example: ** 
    % xData = 0:.01:10; yData = cos(xData); zData = yData.^2;
    %
    % %[3D Example, FULLY SPECIFIED ARGUMENTS LIST]
    %    figure(4); set(gcf,'name', '[3D Example, FULLY SPECIFIED ARGUMENTS LIST]');
    %    plot3(xData,yData,zData,'r','linewidth',5) 
    %    tcmgPlotFormat('Z(x,$\theta(t)$)', 'x', 'y', 'Z', 'dark', 20, 'latex','',0,1,1)
    %
    % %[3D Example, FULLY SPECIFIED ARGUMENTS LIST]
    %    figure(5); set(gcf,'name', '[3D Example, FULLY SPECIFIED ARGUMENTS LIST]');
    %    peaks 
    %    tcmgPlotFormat('PEAKS', 'x', 'y', 'z', 'dark', 20, 'latex','',0,1,1)
    %    colormap(jet)
    %    alpha(0.5)
    %    axis manual
    %***********************************************************
    
    % CUSTOMIZABLE ARGUMENT LIST / FUNCTION OVERLOADING
    %********************************************************************** 
    maxArgs = 12;
    %Automatically fill in blanks and specified args in list with '~'
    for i = 1:nargin
        if(strcmp(varargin{i},''))
            varargin{i} = '~';
        end
    end
    for i = 1:(maxArgs-nargin)
        varargin = [varargin '~'];
    end  
    %Now to check whether a field was unspecified you just need to check if
    %the value in its corresponding position in the arguments list is equal 
    %to '~'
    
    %Store argumens in local variables (You can create new ones and
    %specify what to do with them, just make sure to update "maxArgs")
    %Current Structure:
    %tcmgPlotFormat(tString,xString,yString,zString,style,fontSize,interpreter,legendObj, yRot, showGrid, showColorbar)
    tString = varargin{1};   xString = varargin{2};  
    yString = varargin{3};   zString = varargin{4};     
    style = varargin{5};  fontSize = varargin{6};
    interpreter = varargin{7};  legendObj = varargin{8}; 
    yRot = varargin{9};  showGrid = varargin{10}; 
    showColorbar = varargin{11}; 
    %**********************************************************************
    
    % Background and foreground setup
    if(strcmp(style,'light')) 
        fg_color = 'k';
        bg_color = 'w';
    elseif(strcmp(style,'dark') ||  strcmp(style,'~'))  % Dark is the default style
        fg_color = 'w';
        bg_color = 'k';
    end
    set(gcf,'color',bg_color)  
    set(gca,'color',bg_color,'YColor',fg_color,'XColor',fg_color,'ZColor',fg_color)
 
    % Interpreter & Font Size check
    if (strcmp(interpreter,'~'))
        interpreter = 'latex';
        %interpreter = 'none';
    end   
    if (fontSize == '~')
        fontSize = 10;
        % 20 is the default size for latex text
        if(strcmp(interpreter, 'latex')) 
            fontSize = 20; 
        end
    end
    
    % Title and axis labels
    if(strcmp(tString,'~') == 0)
        title(tString,'interpreter',interpreter,'fontsize',fontSize,'Color',fg_color)
    end
    if(strcmp(xString,'~') == 0)
        xlabel(xString,'interpreter',interpreter,'fontsize',fontSize,'Color',fg_color)
    end
    if(strcmp(yString,'~') == 0)
        ylabel(yString,'interpreter',interpreter,'fontsize',fontSize,'Color',fg_color) 
    end
    if(strcmp(zString,'~') == 0)
        zlabel(zString,'interpreter',interpreter,'fontsize',fontSize,'Color',fg_color)   
    end
    %Grid
    if(showGrid ~= '~' && showGrid == 1)
        grid on
    end
    
     %Color bar
     if(showColorbar ~= '~' && showColorbar == 1)
        colormap(jet)
        C = colorbar;
        set(C,'Color',fg_color) 
     end
     
    %Legend
    if(legendObj ~= '~')
        set(legendObj,'interpreter',interpreter,'fontsize',1*fontSize,'Color', ...
            bg_color, 'EdgeColor',fg_color ,'TextColor',fg_color ) 
    end
    
    %Y-axis label rotation
    if(yRot ~= '~')
        set(get(gca,'YLabel'),'rotation',yRot) 
        %set(get(gca,'Ylabel'),'position',[-0.09562211981566827 9.215409155416902])
    end 
    
    %alpha(0.7)
    %Now add your own!
    % ... 
     
    
end

