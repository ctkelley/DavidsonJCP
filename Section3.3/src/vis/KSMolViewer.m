function varargout = KSMolViewer(varargin)
% KSMOLVIEWER2 MATLAB code for KSMolViewer.fig
%    KSMOLVIEWER2(mol) creates a new KSMOLVIEWER2 with displayment of bonds
%    and atoms.
%
%    KSMOLVIEWER2(mol,H) creates a new KSMOLVIEWER2 with displayment of bonds
%    and atoms together with the density given in Hamiltonian H.
%
%    KSMOLVIEWER2(mol,rho) creates a new KSMOLVIEWER2 with displayment of
%    bonds and atoms together with the density rho.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  Copyright (c) 2014 Chad A. Greene
%  This file is distributed under the terms of the MIT License.

% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KSMolViewer (see VARARGIN)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @KSMolViewer_OpeningFcn, ...
    'gui_OutputFcn',  @KSMolViewer_OutputFcn, ...
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

end
% End initialization code - DO NOT EDIT

% --- Executes just before KSMolViewer is made visible.
function KSMolViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% Choose default command line output for KSMolViewer
handles.output = hObject;

usrdat = [];

% Process input data
if length(varargin) < 1 
    usrdat.msg = 'error';
elseif ischar(varargin{1})
    usrdat.msg = '';
    [usrdat.mol,H] = KSload(varargin{1});
    usrdat.rho = H.rho;
elseif isa(varargin{1},'Molecule')
    usrdat.msg = '';
    usrdat.mol = varargin{1};
    if length(varargin) == 2
        if isa(varargin{2},'Ham') || isa(varargin{2},'BlochHam')
            usrdat.rho = varargin{2}.rho;
        else
            usrdat.rho = varargin{2};
        end
    else
        usrdat.rho = [];
    end
end

handles.output.UserData = usrdat;

% Update handles structure
guidata(hObject, handles);

end

% --- Outputs from this function are returned to the command line.
function varargout = KSMolViewer_OutputFcn(hObject, eventdata, handles)

usrdat = handles.output.UserData;
usrdat.fatoms = 1;
usrdat.fbonds = 1;
usrdat.fdensity = 0;

jColorBackground = java.awt.Color(hObject.Color(1), ...
    hObject.Color(2), hObject.Color(3));

usrdat.iso_slider = javax.swing.JSlider;
javacomponent(usrdat.iso_slider,[15 30 300 35]);
jslider = findjobj(usrdat.iso_slider);

jslider.Background = jColorBackground;
jslider.MajorTickSpacing = 20;
jslider.PaintTicks = 1;
jslider.StateChangedCallback = @(hO,ev)iso_slider_Callback(hO,ev,handles);

usrdat.alpha_slider = javax.swing.JSlider;
javacomponent(usrdat.alpha_slider,[335 25 300 35]);
aslider = findjobj(usrdat.alpha_slider);

aslider.Background = jColorBackground;
aslider.MajorTickSpacing = 20;
aslider.PaintTicks = 1;
aslider.PaintLabels = 1;
aslider.StateChangedCallback = @(hO,ev)alpha_slider_Callback(hO,ev,handles);

if strcmpi(usrdat.msg,'error')
    hf = errordlg('Input error.','Error');
    uiwait(hf);
    hf=findobj('Name','KSMolViewer');
    close(hf);
    return;
end

if isempty(usrdat.rho)
    handles.flag_density.Enable = 'off';
    handles.text_iso_val.Enable = 'off';
    handles.text_alpha_val.Enable = 'off';
    usrdat.iso_slider.disable();
    usrdat.isoval = 0;
    usrdat.alpha_slider.disable();
    usrdat.alphaval = 0;
    usrdat.fcenterdensity = 0;
    plotall(handles.axes_main,usrdat.mol,[],usrdat.isoval,0,0,0);
else
    minrho = ceil(10000*minel(usrdat.rho));
    maxrho = round(10000*maxel(usrdat.rho));
    usrdat.minrho = minrho;
    usrdat.maxrho = maxrho;
    usrdat.iso_slider.setMajorTickSpacing((maxrho-minrho)/5);
    usrdat.iso_slider.setMinimum(minrho);
    usrdat.iso_slider.setMaximum(maxrho);
    usrdat.iso_slider.setValue((minrho+maxrho)/2);
    usrdat.isoval = (minrho+maxrho)/20000;
    handles.edit_iso_val.String = num2str(usrdat.isoval);
    
    usrdat.alpha_slider.setValue(50);
    usrdat.alphaval = 0.5;
    plotall(handles.axes_main,usrdat.mol,usrdat.rho, ...
        usrdat.isoval,usrdat.alphaval,1,1);
    usrdat.iso_slider.disable();
    usrdat.alpha_slider.disable();
    
    usrdat.fcenterdensity = 1;
end

updateplot(handles.axes_main,usrdat.mol,usrdat.rho, ...
    usrdat.isoval,usrdat.alphaval, ...
    usrdat.fatoms,usrdat.fbonds,usrdat.fdensity,usrdat.fcenterdensity);
handles.output.UserData = usrdat;
% Get default command line output from handles structure
varargout{1} = handles.output;

end

function iso_slider_Callback(hObject, eventdata, handles)

handles.output.UserData.isoval = hObject.Value/10000;
usrdat = handles.output.UserData;
handles.edit_iso_val.String = num2str(usrdat.isoval);

if isfield(usrdat,'fatoms')
    updateplot(handles.axes_main,usrdat.mol,usrdat.rho, ...
        usrdat.isoval,usrdat.alphaval, ...
        usrdat.fatoms,usrdat.fbonds,usrdat.fdensity,usrdat.fcenterdensity);
end

end

function alpha_slider_Callback(hObject, eventdata, handles)

handles.output.UserData.alphaval = hObject.Value/100;
usrdat = handles.output.UserData;

if isfield(usrdat,'fatoms')
    updateplot(handles.axes_main,usrdat.mol,usrdat.rho, ...
        usrdat.isoval,usrdat.alphaval, ...
        usrdat.fatoms,usrdat.fbonds,usrdat.fdensity,usrdat.fcenterdensity);
end

end

% --- Executes on button press in flag_atoms.
function flag_atoms_Callback(hObject, eventdata, handles)

handles.output.UserData.fatoms = hObject.Value;
usrdat = handles.output.UserData;

updateplot(handles.axes_main,usrdat.mol,usrdat.rho, ...
    usrdat.isoval,usrdat.alphaval, ...
    usrdat.fatoms,usrdat.fbonds,usrdat.fdensity,usrdat.fcenterdensity);

end


% --- Executes on button press in flag_bonds.
function flag_bonds_Callback(hObject, eventdata, handles)

handles.output.UserData.fbonds = hObject.Value;
usrdat = handles.output.UserData;

updateplot(handles.axes_main,usrdat.mol,usrdat.rho, ...
    usrdat.isoval,usrdat.alphaval, ...
    usrdat.fatoms,usrdat.fbonds,usrdat.fdensity,usrdat.fcenterdensity);

end

% --- Executes on button press in flag_density.
function flag_density_Callback(hObject, eventdata, handles)

handles.output.UserData.fdensity = hObject.Value;
usrdat = handles.output.UserData;

if hObject.Value
    handles.text_iso_val.Enable = 'on';
    handles.edit_iso_val.Enable = 'on';
    handles.text_alpha_val.Enable = 'on';
    handles.flag_center_density.Enable = 'on';
    usrdat.iso_slider.enable();
    usrdat.alpha_slider.enable();
else
    handles.text_iso_val.Enable = 'off';
    handles.edit_iso_val.Enable = 'off';
    handles.text_alpha_val.Enable = 'off';
    handles.flag_center_density.Enable = 'off';
    usrdat.iso_slider.disable();
    usrdat.alpha_slider.disable();
end

updateplot(handles.axes_main,usrdat.mol,usrdat.rho, ...
    usrdat.isoval,usrdat.alphaval, ...
    usrdat.fatoms,usrdat.fbonds,usrdat.fdensity,usrdat.fcenterdensity);

end

function flag_center_density_Callback(hObject, eventdata, handles)

handles.output.UserData.fcenterdensity = hObject.Value;
usrdat = handles.output.UserData;

updateplot(handles.axes_main,usrdat.mol,usrdat.rho, ...
    usrdat.isoval,usrdat.alphaval, ...
    usrdat.fatoms,usrdat.fbonds,usrdat.fdensity,usrdat.fcenterdensity);

end

function edit_iso_val_Callback(hObject, eventdata, handles)

isoval = str2double(get(hObject,'String'));
if isempty(isoval) || isoval*10000 < handles.output.UserData.minrho
    warndlg('Invalid Isosurface Value.');
    hObject.String = num2str(handles.output.UserData.minrho/10000);
elseif isoval*10000 > handles.output.UserData.maxrho
    warndlg('Invalid Isosurface Value.');
    hObject.String = num2str(handles.output.UserData.maxrho/10000);
end

handles.output.UserData.isoval = str2double(get(hObject,'String'));
handles.output.UserData.iso_slider.setValue( ...
    handles.output.UserData.isoval*10000);
usrdat = handles.output.UserData;

updateplot(handles.axes_main,usrdat.mol,usrdat.rho, ...
    usrdat.isoval,usrdat.alphaval, ...
    usrdat.fatoms,usrdat.fbonds,usrdat.fdensity,usrdat.fcenterdensity);

end

% --- Executes during object creation, after setting all properties.
function edit_iso_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_iso_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --------------------------------------------------------------------
function uipushtool_save_figure_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname,~] = uiputfile( ...
    {'*.eps', 'EPS file (*.eps)'; ...
     '*.jpg', 'JPEG image (*.jpg)'; ...
     '*.bmp', 'Bitmap file (*.bmp)'; ...
     '*.pdf', 'Portable Document Format (*.pdf)'; ...
     '*.png', 'Portable Network Graphics file (*.png)'; ...
     '*.tif', 'TIFF image (*.tif)'}, ...
     'Save Image', ...
     handles.output.UserData.mol.name);

if isempty(filename)
    return;
end

export_fig(handles.axes_main,[pathname filename]);

end
