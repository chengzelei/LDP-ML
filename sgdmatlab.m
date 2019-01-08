function varargout = sgdmatlab(varargin)
% SGDMATLAB MATLAB code for sgdmatlab.fig
%      SGDMATLAB, by itself, creates a new SGDMATLAB or raises the existing
%      singleton*.
%
%      H = SGDMATLAB returns the handle to a new SGDMATLAB or the handle to
%      the existing singleton*.
%
%      SGDMATLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SGDMATLAB.M with the given input arguments.
%
%      SGDMATLAB('Property','Value',...) creates a new SGDMATLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sgdmatlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sgdmatlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sgdmatlab

% Last Modified by GUIDE v2.5 02-Jan-2019 20:18:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sgdmatlab_OpeningFcn, ...
                   'gui_OutputFcn',  @sgdmatlab_OutputFcn, ...
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


% --- Executes just before sgdmatlab is made visible.
function sgdmatlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sgdmatlab (see VARARGIN)

% Choose default command line output for sgdmatlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sgdmatlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sgdmatlab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
initialise();
global timerforstop;
handles.timer = timer('Period',1.5,'ExecutionMode','FixedRate', 'TimerFcn',{@printgif, handles}); 
timerforstop=handles.timer;
start( timerforstop);

function printgif(hOject,eventdata,handles)
global x;
global y;
global z;
global X;
global Y;
global Z;
global noisy_grad;
global true_grad;
i=getiternum();
cla(handles.axes1);
hold on;
surf(handles.axes1,X,Y,Z);
shading(handles.axes1,'interp');
colorbar(handles.axes1);

hold(handles.axes1,'on')
plot3(handles.axes1,x(i), y(i), z(i), 'ro-', 'Linewidth', 2);
rotate3d on;
rotate3d(handles.axes1)
hold on;

axes(handles.axes3);
A=imread('aggre1.png');
B=imresize(A,[256 256]);
imshow(B,'InitialMagnification','fit');

set(handles.text2,'String',sprintf('%s%d','Iteration num = ',i));   
set(handles.text5,'String',sprintf('%s%.3d','Nosiy average gradient = ',noisy_grad(i)));   
set(handles.text4,'String',sprintf('%s%.3d','True average gradient = ',true_grad(i)));   

function initialise()
global x;
global y;
global z;
global X;
global Y;
global Z;
global true_grad;
global noisy_grad;
x1 = load('betao.txt');
y1 = load('betat.txt');
z1 = load('loss.txt');
noisy_grad=load('noisy_grad.txt');
true_grad=load('true_grad.txt');
% 抽稀，以免内存不够
count    = 1;     % 新变量计数器
interval = 20; % 抽稀间隔
for i = 1 : interval : 10000
    x(count) = x1(i);
    y(count) = y1(i);
    z(count) = z1(i);
    count = count + 1;
end
[X,Y]=meshgrid(min(x):0.1:max(x),min(y):0.1:max(y)); 
%在网格点位置插值求Z，注意：不同的插值方法得到的曲线光滑度不同
Z=griddata(x,y,z,X,Y,'v4');
hold on;
figure('Visible','off')

function iternum=getiternum()
persistent a;
if isempty(a)
a=0;
end
a=a+1;
iternum=a;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global timerforstop;
stop(timerforstop );

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global timerforstop;
start(timerforstop );

% --- Executes on key press with focus on pushbutton2 and none of its controls.
function pushbutton2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
