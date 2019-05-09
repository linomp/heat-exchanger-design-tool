function varargout = Proyecto_transfer2(varargin)
% PROYECTO_TRANSFER2 MATLAB code for Proyecto_transfer2.fig
%      PROYECTO_TRANSFER2, by itself, creates a new PROYECTO_TRANSFER2 or raises the existing
%      singleton*.
%
%      H = PROYECTO_TRANSFER2 returns the handle to a new PROYECTO_TRANSFER2 or the handle to
%      the existing singleton*.
%
%      PROYECTO_TRANSFER2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROYECTO_TRANSFER2.M with the given input arguments.
%
%      PROYECTO_TRANSFER2('Property','Value',...) creates a new PROYECTO_TRANSFER2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Proyecto_transfer2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Proyecto_transfer2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Proyecto_transfer2

% Last Modified by GUIDE v2.5 25-Apr-2019 13:04:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Proyecto_transfer2_OpeningFcn, ...
                   'gui_OutputFcn',  @Proyecto_transfer2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Proyecto_transfer2 is made visible.
function Proyecto_transfer2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Proyecto_transfer2 (see VARARGIN)

% Choose default command line output for Proyecto_transfer2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Proyecto_transfer2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Proyecto_transfer2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function T_Tube_Callback(hObject, eventdata, handles)
% hObject    handle to T_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_Tube as text
%        str2double(get(hObject,'String')) returns contents of T_Tube as a double


% --- Executes during object creation, after setting all properties.
function T_Tube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function D_Tube_Callback(hObject, eventdata, handles)
% hObject    handle to D_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D_Tube as text
%        str2double(get(hObject,'String')) returns contents of D_Tube as a double


% --- Executes during object creation, after setting all properties.
function D_Tube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cp_Tube_Callback(hObject, eventdata, handles)
% hObject    handle to Cp_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cp_Tube as text
%        str2double(get(hObject,'String')) returns contents of Cp_Tube as a double


% --- Executes during object creation, after setting all properties.
function Cp_Tube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cp_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Q_Tube_Callback(hObject, eventdata, handles)
% hObject    handle to Q_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Q_Tube as text
%        str2double(get(hObject,'String')) returns contents of Q_Tube as a double


% --- Executes during object creation, after setting all properties.
function Q_Tube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function D_Shell_Callback(hObject, eventdata, handles)
% hObject    handle to D_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D_Shell as text
%        str2double(get(hObject,'String')) returns contents of D_Shell as a double


% --- Executes during object creation, after setting all properties.
function D_Shell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_Shell_Callback(hObject, eventdata, handles)
% hObject    handle to T_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_Shell as text
%        str2double(get(hObject,'String')) returns contents of T_Shell as a double


% --- Executes during object creation, after setting all properties.
function T_Shell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cp_Shell_Callback(hObject, eventdata, handles)
% hObject    handle to Cp_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cp_Shell as text
%        str2double(get(hObject,'String')) returns contents of Cp_Shell as a double


% --- Executes during object creation, after setting all properties.
function Cp_Shell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cp_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Q_Shell_Callback(hObject, eventdata, handles)
% hObject    handle to Q_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Q_Shell as text
%        str2double(get(hObject,'String')) returns contents of Q_Shell as a double


% --- Executes during object creation, after setting all properties.
function Q_Shell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Solve.
function Solve_Callback(hObject, eventdata, handles)
% hObject    handle to Solve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes1)
handles.Solve.Enable = 'off';
handles.statusLbl.String = 'Calculating...';

hold off
handles.resShellDiameter.String = '...';
handles.resTubeDiameter.String = '...';
handles.resBaffleSpacing.String = '...';
handles.resNumPasses.String = '...';
handles.resNumTubes.String = '...';
handles.resCost.String = '...';
handles.resQ.String = '...';


energy_cost = str2double(handles.energy_cost.String);
p_efficiency = str2double(handles.p_efficiency.String);
operating_time = str2double(handles.operating_time.String); %[hours per week]

operating_time = operating_time*12*4 %[hours per year]

M_Tube = str2double(handles.M_Tube.String);
M_Shell = str2double(handles.M_Shell.String);

T_Tube = str2double(handles.T_Tube.String);
T2_Tube = str2double(handles.T2_Tube.String);
workingTempTube = (T_Tube+T2_Tube)/2;
autoFillTubeProps(workingTempTube);

D_Tube = str2double(handles.D_Tube.String);
Cp_Tube = str2double(handles.Cp_Tube.String);
Q_Tube = str2double(handles.Q_Tube.String); %L/h

Visc_Tube = str2double(handles.Visc_Tube.String);
k_Tube = str2double(handles.k_Tube.String);

Q_Tube = Q_Tube/(3600*1000) %m3/s
massflow_Tube = Q_Tube * D_Tube

T_Shell = str2double(handles.T_Shell.String);
T2_Shell = str2double(handles.T2_Shell.String);
workingTempShell = (T_Shell+T2_Shell)/2;
autoFillTubeProps(workingTempShell);

D_Shell = str2double(handles.D_Shell.String);
Cp_Shell = str2double(handles.Cp_Shell.String);
Q_Shell = str2double(handles.Q_Shell.String);
Visc_Shell = str2double(handles.Visc_Shell.String);
k_Shell = str2double(handles.k_Shell.String);

Q_Shell = (Q_Shell*0.001)/3600 %m3/s
massflow_Shell = Q_Shell * D_Shell

shell_out_diameter_min = str2double(handles.shell_out_diameter_min.String);
shell_out_diameter_max = str2double(handles.shell_out_diameter_max.String);
tube_out_diameter_min = str2double(handles.tube_out_diameter_min.String);
tube_out_diameter_max = str2double(handles.tube_out_diameter_max.String);
baffle_spacing_min = str2double(handles.baffle_spacing_min.String);
baffle_spacing_max = str2double(handles.baffle_spacing_max.String);
num_passes_min = str2double(handles.num_passes_min.String);
num_passes_max = str2double(handles.num_passes_max.String);
num_tubes_min = str2double(handles.num_tubes_min.String);
num_tubes_max = str2double(handles.num_tubes_max.String);


parentAxis = handles.axes1; 

% User input
% -- shell
data.mass_flow_shell = massflow_Shell; %[kg/s]
data.T_shell_in = T_Shell+273; %[K] 
data.T_shell_out = T2_Shell+273; %[K]


%[rho, cp, visc, k] = getWaterProperties(workingTemp);

data.rho_shell = D_Shell;%850; %[kg/m^3]
data.cp_shell = Cp_Shell/1000;%2.47; % [KJ/kg]
data.shell_viscosity = Visc_Shell;% 0.0004; %[Pa-s]
data.k_shell = k_Shell; % 0.13; %[W/mK ?] 


% -- tube
data.mass_flow_tube = massflow_Tube; %[kg/s]
data.T_tube_in = T_Tube+273; %[K]
data.T_tube_out = T2_Tube+273; %[K]
data.rho_tube = D_Tube;%995; %[kg/m^3]
data.cp_tube = Cp_Tube/1000;%2.05; % [KJ/kg] 
data.tube_viscosity = Visc_Tube;%0.00358; %[Pa-s]
data.k_tube = k_Tube; %0.13; %[W/mK ?]
% -- energy
data.energyCost = energy_cost; % $/kWh energy cost
data.annualOperatingTimeInHours = operating_time;  % annual operating time in hours
data.pumpEfficiency = p_efficiency; % Pump effficiency  
% -- bounds
data.LB = [shell_out_diameter_min, tube_out_diameter_min, baffle_spacing_min, num_passes_min, num_tubes_min];       
data.UB = [shell_out_diameter_max, tube_out_diameter_max, baffle_spacing_max, num_passes_max, num_tubes_max]; 

pause(0.25);
[bestParameters]  = mainFunc(data, parentAxis);
handles.statusLbl.String = 'Done';

handles.resShellDiameter.String = sprintf('%.3f', bestParameters.r_shell*2);
handles.resTubeDiameter.String = sprintf('%.3f', bestParameters.r_tube*2);
handles.resBaffleSpacing.String = sprintf('%.3f', bestParameters.baffle);
handles.resNumPasses.String = bestParameters.n_passes;
handles.resNumTubes.String = bestParameters.n_tubes; 
handles.resCost.String = sprintf('%.3f', bestParameters.cost);
handles.resQ.String = sprintf('%.3f', bestParameters.Q);

handles.Solve.Enable = 'on';



function T2_Tube_Callback(hObject, eventdata, handles)
% hObject    handle to T2_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T2_Tube as text
%        str2double(get(hObject,'String')) returns contents of T2_Tube as a double


% --- Executes during object creation, after setting all properties.
function T2_Tube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T2_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T2_Shell_Callback(hObject, eventdata, handles)
% hObject    handle to T2_Shell (see GCBO)
% eventdata  reserved - to be defixned in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T2_Shell as text
%        str2double(get(hObject,'String')) returns contents of T2_Shell as a double


% --- Executes during object creation, after setting all properties.
function T2_Shell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T2_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in M_Tube.
function M_Tube_Callback(hObject, eventdata, handles)
% hObject    handle to M_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns M_Tube contents as cell array
%        contents{get(hObject,'Value')} returns selected item from M_Tube
contents = cellstr(get(handles.M_Tube,'String'));
fluido = contents{get(handles.M_Tube,'Value')}

if strcmp(fluido,'Water')
   
    autoFillTubeProps(25, handles);
    
elseif strcmp(fluido,'Crude Oil')
    
    handles.D_Tube.String = 995;
    handles.Cp_Tube.String = 2050;
    handles.Visc_Tube.String = 0.00358;
    handles.k_Tube.String = 0.13;
    handles.T_Tube.String = 37.8;
    handles.T2_Tube.String = 76.7;
    handles.Q_Tube.String = 68020.10;
    
elseif strcmp(fluido,'Kerosene')
   
    handles.D_Tube.String = 850;
    handles.Cp_Tube.String = 2470;
    handles.Visc_Tube.String = 0.0004;
    handles.k_Tube.String = 0.13;
    handles.T_Tube.String = 199;
    handles.T2_Tube.String = 93.3;
    handles.Q_Tube.String = 23378.82;
    
else
    
    limpiarTube(handles);
    % TO DO: Limpiar todos los campos de textos
    
end


% --- Executes during object creation, after setting all properties.
function M_Tube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in M_Shell.
function M_Shell_Callback(hObject, eventdata, handles)
% hObject    handle to M_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns M_Shell contents as cell array
%        contents{get(hObject,'Value')} returns selected item from M_Shell
contents = cellstr(get(handles.M_Shell,'String'));
fluido = contents{get(handles.M_Shell,'Value')}

if strcmp(fluido,'Water')
   
    autoFillShellProps(25,handles);
    
elseif strcmp(fluido,'Crude Oil')
    
    handles.D_Shell.String = 995;
    handles.Cp_Shell.String = 2050;
    handles.Visc_Shell.String = 0.00358;
    handles.k_Shell.String = 0.13;
    handles.T_Shell.String = 37.8;
    handles.T2_Shell.String = 76.7;
    handles.Q_Shell.String = 68020.10;
    
elseif strcmp(fluido,'Kerosene')
   
    handles.D_Shell.String = 850;
    handles.Cp_Shell.String = 2470;
    handles.Visc_Shell.String = 0.0004;
    handles.k_Shell.String = 0.13;
    handles.T_Shell.String = 199;
    handles.T2_Shell.String = 93.3;
    handles.Q_Shell.String = 23378.82;
        
else
    
    limpiarShell(handles);
    
end

% --- Executes during object creation, after setting all properties.
function M_Shell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Visc_Tube_Callback(hObject, eventdata, handles)
% hObject    handle to Visc_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Visc_Tube as text
%        str2double(get(hObject,'String')) returns contents of Visc_Tube as a double


% --- Executes during object creation, after setting all properties.
function Visc_Tube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Visc_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Visc_Shell_Callback(hObject, eventdata, handles)
% hObject    handle to Visc_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Visc_Shell as text
%        str2double(get(hObject,'String')) returns contents of Visc_Shell as a double


% --- Executes during object creation, after setting all properties.
function Visc_Shell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Visc_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function shell_out_diameter_min_Callback(hObject, eventdata, handles)
% hObject    handle to shell_out_diameter_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shell_out_diameter_min as text
%        str2double(get(hObject,'String')) returns contents of shell_out_diameter_min as a double


% --- Executes during object creation, after setting all properties.
function shell_out_diameter_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shell_out_diameter_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tube_out_diameter_min_Callback(hObject, eventdata, handles)
% hObject    handle to tube_out_diameter_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tube_out_diameter_min as text
%        str2double(get(hObject,'String')) returns contents of tube_out_diameter_min as a double


% --- Executes during object creation, after setting all properties.
function tube_out_diameter_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tube_out_diameter_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baffle_spacing_min_Callback(hObject, eventdata, handles)
% hObject    handle to baffle_spacing_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baffle_spacing_min as text
%        str2double(get(hObject,'String')) returns contents of baffle_spacing_min as a double


% --- Executes during object creation, after setting all properties.
function baffle_spacing_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baffle_spacing_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_passes_min_Callback(hObject, eventdata, handles)
% hObject    handle to num_passes_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_passes_min as text
%        str2double(get(hObject,'String')) returns contents of num_passes_min as a double


% --- Executes during object creation, after setting all properties.
function num_passes_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_passes_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_tubes_min_Callback(hObject, eventdata, handles)
% hObject    handle to num_tubes_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_tubes_min as text
%        str2double(get(hObject,'String')) returns contents of num_tubes_min as a double


% --- Executes during object creation, after setting all properties.
function num_tubes_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_tubes_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function shell_out_diameter_max_Callback(hObject, eventdata, handles)
% hObject    handle to shell_out_diameter_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shell_out_diameter_max as text
%        str2double(get(hObject,'String')) returns contents of shell_out_diameter_max as a double


% --- Executes during object creation, after setting all properties.
function shell_out_diameter_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shell_out_diameter_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tube_out_diameter_max_Callback(hObject, eventdata, handles)
% hObject    handle to tube_out_diameter_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tube_out_diameter_max as text
%        str2double(get(hObject,'String')) returns contents of tube_out_diameter_max as a double


% --- Executes during object creation, after setting all properties.
function tube_out_diameter_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tube_out_diameter_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baffle_spacing_max_Callback(hObject, eventdata, handles)
% hObject    handle to baffle_spacing_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baffle_spacing_max as text
%        str2double(get(hObject,'String')) returns contents of baffle_spacing_max as a double


% --- Executes during object creation, after setting all properties.
function baffle_spacing_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baffle_spacing_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_passes_max_Callback(hObject, eventdata, handles)
% hObject    handle to num_passes_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_passes_max as text
%        str2double(get(hObject,'String')) returns contents of num_passes_max as a double


% --- Executes during object creation, after setting all properties.
function num_passes_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_passes_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_tubes_max_Callback(hObject, eventdata, handles)
% hObject    handle to num_tubes_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_tubes_max as text
%        str2double(get(hObject,'String')) returns contents of num_tubes_max as a double


% --- Executes during object creation, after setting all properties.
function num_tubes_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_tubes_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function autoFillTubeProps(T,handles)

%handles.D_Tube, handles.Cp_Tube, handles.Q_Tube, handles.Visc_Tube
[rho, cp, visc, k] = getWaterProperties(T)
handles.D_Tube.String = rho;
handles.Cp_Tube.String = cp;
handles.Visc_Tube.String = visc;
handles.k_Tube.String = k;


function autoFillShellProps(T,handles)

%handles.D_Tube, handles.Cp_Tube, handles.Q_Tube, handles.Visc_Tube
    [rho, cp, visc, k] = getWaterProperties(T)
    handles.D_Shell.String = rho;
    handles.Cp_Shell.String = cp;
    handles.Visc_Shell.String = visc;
    handles.k_Shell.String = k;

function limpiarTube(handles)

handles.D_Tube.String = '';
handles.Cp_Tube.String = '';
handles.Visc_Tube.String = '';
handles.k_Tube.String = '';


function limpiarShell(handles)

    handles.D_Shell.String = '';
    handles.Cp_Shell.String = '';
    handles.Visc_Shell.String = '';
    handles.k_Shell.String = '';



function k_Tube_Callback(hObject, eventdata, handles)
% hObject    handle to k_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_Tube as text
%        str2double(get(hObject,'String')) returns contents of k_Tube as a double


% --- Executes during object creation, after setting all properties.
function k_Tube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_Tube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function operating_time_Callback(hObject, eventdata, handles)
% hObject    handle to operating_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of operating_time as text
%        str2double(get(hObject,'String')) returns contents of operating_time as a double


% --- Executes during object creation, after setting all properties.
function operating_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to operating_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_efficiency_Callback(hObject, eventdata, handles)
% hObject    handle to p_efficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p_efficiency as text
%        str2double(get(hObject,'String')) returns contents of p_efficiency as a double


% --- Executes during object creation, after setting all properties.
function p_efficiency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_efficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function energy_cost_Callback(hObject, eventdata, handles)
% hObject    handle to energy_cost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of energy_cost as text
%        str2double(get(hObject,'String')) returns contents of energy_cost as a double


% --- Executes during object creation, after setting all properties.
function energy_cost_CreateFcn(hObject, eventdata, handles)
% hObject    handle to energy_cost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k_Shell_Callback(hObject, eventdata, handles)
% hObject    handle to k_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_Shell as text
%        str2double(get(hObject,'String')) returns contents of k_Shell as a double


% --- Executes during object creation, after setting all properties.
function k_Shell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_Shell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
