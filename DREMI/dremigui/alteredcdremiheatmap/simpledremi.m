function varargout = simpledremi(varargin)
% SIMPLEDREMI MATLAB code for simpledremi.fig
%      SIMPLEDREMI, by itself, creates a new SIMPLEDREMI or raises the existing
%      singleton*.
%
%      H = SIMPLEDREMI returns the handle to a new SIMPLEDREMI or the handle to
%      the existing singleton*.
%
%      SIMPLEDREMI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMPLEDREMI.M with the given input arguments.
%
%      SIMPLEDREMI('Property','Value',...) creates a new SIMPLEDREMI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simpledremi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to simpledremi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help simpledremi

% Last Modified by GUIDE v2.5 18-Jan-2016 13:12:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simpledremi_OpeningFcn, ...
                   'gui_OutputFcn',  @simpledremi_OutputFcn, ...
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

  

end

% -- returns saved variable or empty matrix if the variable is not found.
function var = retr(name)
    hgui=getappdata(0,'hgui');  
    var=getappdata(hgui, name);
end

function handles=gethand
    hgui=getappdata(0,'hgui');
    handles=guihandles(hgui);
end

function put(name, what)
    hgui=getappdata(0,'hgui');
    setappdata(hgui, name, what);
    
    listener = retr([name '_listener']);
    if ~isempty(listener)
        listener();
    end
end


% --- Executes just before simpledremi is made visible.
function simpledremi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simpledremi (see VARARGIN)

% Choose default command line output for simpledremi
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

setappdata(0,'hgui',gcf);


% UIWAIT makes simpledremi wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = simpledremi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function placeIcon(hObject, filename);
    try 
        [image_pic, map] = imread(filename);
        if ~isempty( map ) 
            image_pic = ind2rgb( image_pic, map ); 
        end 
        set(hObject,'cdata',image_pic); 
    catch e
        fprintf('warning: failed reading icon image ''%s''', filename);
    end

end

% --- Executes on selection change in lst_Channels.
function lst_Channels_Callback(hObject, eventdata, handles)
% hObject    handle to lst_Channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lst_Channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lst_Channels
      
      sessionData = retr('sessionData');
      %get the channels that are selected
      selected_channels = get(handles.lst_Channels,'Value');
      channels = cell(1, length(selected_channels));
      
      for i=1:length(selected_channels)
          channels{i} = sessionData.marker_channels{selected_channels(i)};
      end
      
      num_selected_channels = length(selected_channels);
      sig_edges_string = cell(num_selected_channels * (num_selected_channels-1),1);
      current_loc = 1;
      for i=1:length(channels)
          
          for j=1:length(channels)
              
              if(i==j)
                  continue;
              end
              
              sig_edges_string{current_loc}= sprintf('%s -> %s', channels{i}, channels{j});
              current_loc = current_loc+1;
          end
      end
      set(handles.lst_sigedges,'String',sig_edges_string);
      

%       num_sig_edges = size(cdata.DREMI_sig_edges,1);
%       sig_edges_string = cell(num_sig_edges,1);
%       for i=1:num_sig_edges
%           sig_edges_string{i} = sprintf('%s->%s',cdata.DREMI_sig_edges{i,1},cdata.DREMI_sig_edges{i,2});
%       end
%           
%       
%       set(handles.lst_sigedges,'String', sig_edges_string);




end

% --- Executes during object creation, after setting all properties.
function lst_Channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lst_Channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end



function Edt_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Edt_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edt_threshold as text
%        str2double(get(hObject,'String')) returns contents of Edt_threshold as a double

    sessionData = retr('sessionData');
    
    threshold = str2num(get(handles.Edt_threshold,'String'));
    num_datasets = length(sessionData.datasets);
    for i=1:num_datasets
        
        sessionData.datasets{i}.DREMI_noise_threshold = threshold;
        
    end
    
    put('sessionData', sessionData);
end

% --- Executes during object creation, after setting all properties.
function Edt_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edt_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on selection change in popup_Xaxis.
function popup_Xaxis_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Xaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Xaxis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Xaxis


end


% --- Executes during object creation, after setting all properties.
function popup_Xaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Xaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on selection change in popup_Yaxis.
function popup_Yaxis_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Yaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Yaxis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Yaxis
    %get the new value 
 
 %can I get the string associated with that one? 

end

% --- Executes during object creation, after setting all properties.
function popup_Yaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Yaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on selection change in lst_FIles.
function lst_FIles_Callback(hObject, eventdata, handles)
% hObject    handle to lst_FIles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lst_FIles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lst_FIles

        %redraw it on these particular cells. 
        
end


% --- Executes during object creation, after setting all properties.
function lst_FIles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lst_FIles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on button press in btn_AddFile.
function btn_AddFile_Callback(hObject, eventdata, handles)
% hObject    handle to btn_AddFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    OpenMenuItem_Callback(hObject, eventdata,handles);
end



function OpenMenuItem_Callback(hObject,eventdata, handles)
   
    
    currentfolder = retr('currentfolder');
    if (isempty(currentfolder)) 
        currentfolder = '';
    end

    
    %/Users/mtadmor/Data/Cytof/AML_Ateam/
    files = uipickfiles('num',[1 inf],'out','cell', 'FilterSpec', [currentfolder '*.fcs'])
    if isequal(files,0) ~=0
        return 
    end
    
    

    [path, ~, ext] = fileparts(files{1});
    put('currentfolder', [path filesep]);
    
    tic
    
    [cytof_datasets]=cellfun(@cytof_data, files, 'UniformOutput', false);
    for i=1:length(cytof_datasets)
       
        cytof_datasets{i} = cytof_datasets{i};
    end
    
    disp(sprintf('Files loaded: %gs',toc));
    
   
        sessionData  = retr('sessionData');
        if (isempty(sessionData)) 
            sessionData = cytof_series(cytof_datasets, files,'cytof_series');
            
        else 
            sessionData = sessionData.add_datasets(cytof_datasets, files); 
        end
        
    put('sessionData', sessionData);
    set(handles.lst_FIles,'Value',[]); % Add this line so that the list can be changed
    set(handles.lst_FIles,'String', sessionData.datalabels);
    
    set(handles.lst_Channels,'Value',[]); % Add this line so that the list can be changed
    set(handles.lst_Channels, 'String',sessionData.marker_channels);
    
    
    %the popup menu strings should be changed here too 
    set(handles.popup_Xaxis,'Value',1);
    set(handles.popup_Xaxis, 'String',sessionData.marker_channels);
    
    set(handles.popup_Yaxis,'Value',1);
    set(handles.popup_Yaxis, 'String',sessionData.marker_channels);
    
%     set(handles.lst_FIles,'Value',1);
%     set(handles.lstChannels,'Value',[]);
%     
%     lstGates_Callback;
%     lstChannels_Callback;
      
   
end


% --- Executes on button press in btn_DREVI.
function btn_DREVI_Callback(hObject, eventdata, handles)
% hObject    handle to btn_DREVI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

        all_channels = get(handles.popup_Yaxis,'String');
        selected_Ychannel = get(handles.popup_Yaxis,'Value');
        Ychannel = all_channels{selected_Ychannel};
        selected_Xchannel = get(handles.popup_Xaxis,'Value');
        Xchannel = all_channels{selected_Xchannel};
 
        selected_gates = get(handles.lst_FIles, 'Value');
        sessionData = retr('sessionData');
        
        %axes(handles.axes1);
       
       cla;
       Yset = get(handles.Popup_autoY,'Value');
       if(Yset ==1)
            sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits');
       else
           
           maxyval = str2num(get(handles.Edt_maxy,'String'));
           sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits', 'maxy', maxyval);
           
       end
end

% --- Executes on button press in btn_DREMI.
function btn_DREMI_Callback(hObject, eventdata, handles)
% hObject    handle to btn_DREMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%for each subplot

        all_channels = get(handles.popup_Yaxis,'String');
        selected_Ychannel = get(handles.popup_Yaxis,'Value');
        Ychannel = all_channels{selected_Ychannel};
        selected_Xchannel = get(handles.popup_Xaxis,'Value');
        Xchannel = all_channels{selected_Xchannel};
        threshold = str2num(get(handles.Edt_threshold,'String'));

    selected_gates = get(handles.lst_FIles, 'Value');
    sessionData = retr('sessionData');
    Yset = get(handles.Popup_autoY,'Value');
    if(Yset ==1)
        [dremi_values] = sessionData.compute_series_mi( Xchannel, Ychannel, selected_gates,threshold);
    else
         maxyval = str2num(get(handles.Edt_maxy,'String'));
        [dremi_values] = sessionData.compute_series_mi( Xchannel, Ychannel, selected_gates,threshold,'maxy', maxyval);
    end
    for i=1:length(selected_gates)
       
        subplot(1,length(selected_gates),i);
        dremi_string = sprintf('dremi: %.2f',dremi_values(i));
        title(dremi_string);
    end

end


% --- Executes on selection change in popup_Response.
function popup_Response_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Response contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Response

%choose to plot with a response function 

     

     
     
     

end

% --- Executes during object creation, after setting all properties.
function popup_Response_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in btn_Response.
function btn_Response_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

      selected_response = get(handles.popup_Response,'Value')
      all_channels = get(handles.popup_Yaxis,'String');
        selected_Ychannel = get(handles.popup_Yaxis,'Value');
        Ychannel = all_channels{selected_Ychannel};
        selected_Xchannel = get(handles.popup_Xaxis,'Value');
        Xchannel = all_channels{selected_Xchannel};
 
        selected_gates = get(handles.lst_FIles, 'Value');
        sessionData = retr('sessionData');
     
     if(selected_response == 1)
       
        %hold on;  
       cla;
       
        sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits','sigmoid');
        %hold off; 
     elseif(selected_response==2)
        %hold on;
       cla;
      
        sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits','linear');
        %hold off; 
     elseif(selected_response==3)
        %hold on;
       cla;
        sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits','condmean');
        %hold off;
%         
%      elseif(selected_response==4)
%         %hold on;
%        cla;      
%         sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits','smooth_condmean');
%         %hold off; 
     end
end

% % --- Executes on button press in btn_Response.
% function btn_CondMean_Callback(hObject, eventdata, handles)
% % hObject    handle to btn_Response (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
%       selected_response = get(handles.popup_Response,'Value');
%       all_channels = get(handles.popup_Yaxis,'String');
%       selected_Ychannel = get(handles.popup_Yaxis,'Value');
%       Ychannel = all_channels{selected_Ychannel};
%       selected_Xchannel = get(handles.popup_Xaxis,'Value');
%       Xchannel = all_channels{selected_Xchannel};
%  
%       selected_gates = get(handles.lst_FIles, 'Value');
%       sessionData = retr('sessionData');
%      
%       if(selected_response == 1)
%        
%         hold on;  
%        %cla;
%        
%         sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits','cond_mean');
%         %hold off; 
% %      elseif(selected_response==2)
% %         %hold on;
% %        cla;
% %       
% %         sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits','linear');
% %         %hold off; 
% %      elseif(selected_response==3)
% %         %hold on;
% %        cla;
% %         sessionData.pairwise_correlation_visualization_datasets(Xchannel, Ychannel, selected_gates, 'fix_limits','condmean');
% %         %hold off;    
%      end
% end


% --- Executes on button press in btn_moveup.
function btn_moveup_Callback(hObject, eventdata, handles)
% hObject    handle to btn_moveup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     
      selected_gates = get(handles.lst_FIles, 'Value');
      if((length(selected_gates)>1)|(length(selected_gates)==0))
          return;
      end
      
     selected_gate = selected_gates(1);
      
     if(selected_gate==1)
         return;
     end
     sessionData = retr('sessionData');
     cdata1 = sessionData.datasets{selected_gate-1};
     sessionData.datasets{selected_gate-1} = sessionData.datasets{selected_gate};
     sessionData.datasets{selected_gate} = cdata1;
     
     %do the same with datalabels
     label1 = sessionData.datalabels{selected_gate-1};
     sessionData.datalabels{selected_gate-1} = sessionData.datalabels{selected_gate};
     sessionData.datalabels{selected_gate} = label1;
     
     set(handles.lst_FIles,'String', sessionData.datalabels);
     set(handles.lst_FIles,'Value',[selected_gate-1]);
     put('sessionData', sessionData);
     
end

function btn_moveup_CreateFcn(hObject, eventdata, handles)
    placeIcon(hObject, 'arrowup.gif');
end

function btn_movedown_CreateFcn(hObject, eventdata, handles)
    placeIcon(hObject, 'arrowup.gif');
end


% --- Executes on button press in btn_movedown.
function btn_movedown_Callback(hObject, eventdata, handles)
% hObject    handle to btn_movedown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
      selected_gates = get(handles.lst_FIles, 'Value');
      if((length(selected_gates)>1)|(length(selected_gates)==0))
          return;
      end
      
     selected_gate = selected_gates(1);
     sessionData = retr('sessionData');
      
     if(selected_gate==length(sessionData.datalabels))
         return;
     end
   
     cdata1 = sessionData.datasets{selected_gate+1};
     sessionData.datasets{selected_gate+1} = sessionData.datasets{selected_gate};
     sessionData.datasets{selected_gate} = cdata1;
     
     %do the same with datalabels
     label1 = sessionData.datalabels{selected_gate+1};
     sessionData.datalabels{selected_gate+1} = sessionData.datalabels{selected_gate};
     sessionData.datalabels{selected_gate} = label1;
     
     set(handles.lst_FIles,'String', sessionData.datalabels);
     set(handles.lst_FIles,'Value',[selected_gate+1]);
     put('sessionData', sessionData);
end


% --- Executes on button press in Btn_PairwiseDREMI.
function Btn_PairwiseDREMI_Callback(hObject, eventdata, handles)
% hObject    handle to Btn_PairwiseDREMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
      
       s= sprintf('in pairwise callback \n');
      %get the dataset
      selected_gates = get(handles.lst_FIles, 'Value');
      if((length(selected_gates)>1)|(length(selected_gates)==0))
	    matrix_empty=[]
	    for i=1:length(selected_gates)
		  %disp(class(selected_gates));
	      selected_gate = selected_gates(i);
		  %selected_gate= selected_gate(7:strlength(selected_gate));
          sessionData = retr('sessionData');
          cdata = sessionData.datasets{selected_gate};
      
          selected_channels = get(handles.lst_Channels,'Value');
          channels = cell(1, length(selected_channels));
      
          for i=1:length(selected_channels)
            channels{i} = sessionData.marker_channels{selected_channels(i)};
          end     
          noise_threshold = str2num(get(handles.Edt_threshold,'String'));
          DREMI_cutoff = str2num(get(handles.Etd_DREMICutoff,'String'));      
          cdatactrl = sessionData.make_combined_cytof_data_object();      
          [cdata, mi_matrix ] = cdata.pairwise_mi_compute(channels, noise_threshold, 'ctrl_data', cdatactrl);          
          load continuous_BuPu9.mat
          colormap(continuous_BuPu9)
          CLIM = [0 1];
          cla;
          subplot(1,1,1);
		  mi_matrix= reshape(mi_matrix.',1,[])
	      matrix_empty= vertcat(matrix_empty, mi_matrix)
		writematrix(matrix_empty,'MMM.csv');		
      	end
        else  
		  
      
        selected_gate = selected_gates(1);
        sessionData = retr('sessionData');
        cdata = sessionData.datasets{selected_gate};
      
      %get the channels that are selected
      
        selected_channels = get(handles.lst_Channels,'Value');
        channels = cell(1, length(selected_channels));
      
        for i=1:length(selected_channels)
		  disp(channels{i})
          channels{i} = sessionData.marker_channels{selected_channels(i)};
		  xyz= sessionData.marker_channels{selected_channels(i)};
		  xyz= xyz(7:strlength(xyz))
		  disp(xyz)
		  %channels{i} = xyz;
		  
		  %jk= char(channels{i});
		  %jk= jk(4:strlength(jk));
		  %disp(jk);
		  %channels{i}= jk;
		  %channels{i} = channels{i}(4:strlength(channels{i}))
		  %channels{i} = sessionData.marker_channels{selected_channels(i)}(4:strlength(sessionData.marker_channels{selected_channels(i)}));
        end
     
        disp(channels)
		disp(channels{1})
		%channels{1}= char("abc")
     
        noise_threshold = str2num(get(handles.Edt_threshold,'String'));
        DREMI_cutoff = str2num(get(handles.Etd_DREMICutoff,'String'));
      
        cdatactrl = sessionData.make_combined_cytof_data_object();
        
		
		
      
    
        [cdata, mi_matrix ] = cdata.pairwise_mi_compute(channels, noise_threshold, 'ctrl_data', cdatactrl); 
         
		 
		for i=1:length(selected_channels)
		  %disp(channels{i})
          %channels{i} = sessionData.marker_channels{selected_channels(i)};
		  xyz= sessionData.marker_channels{selected_channels(i)};
		  xyz= xyz(4:strlength(xyz))
		  xyz = upper(xyz)
		  channels{i} = xyz;
		  
		end
     
    
        load continuous_BuPu9.mat
        colormap(continuous_BuPu9)
        CLIM = [0 1];
        cla;
        subplot(1,1,1);
	    %TTT = array2table(mi_matrix);
	    %TTT.Properties.VariableNames = channels
	    %TTT.Properties.RowNames = channels;
	    %disp(TTT);
	    %matrix_new= table2array(TTT);
	    %Z= linkage(matrix_new)
	    %dendrogram(Z)
	    %heatmap_1= clustergram(matrix_new, 'Colormap',redbluecmap);
	    %set(heatmap_1,'RowLabels', channels,'ColumnLabels',channels);
	  
	  
	  
	    writematrix(mi_matrix,'MMM.csv');
        imagesc(mi_matrix);
		
        hold on;
        for i = 1:length(channels)
            plot([.5,length(channels)+0.5],[i-.5,i-.5],'k-');
            plot([i-.5,i-.5],[.5,length(channels)+0.5],'k-');
        end
        set(gca, 'Layer','top');
		
        set(gca,'fontsize',10);
        set(gca,'FontWeight','bold');
        set(gca,'ytick',1:length(channels));
        set(gca,'yticklabel',channels);
        disp(channels);
        set(gca,'xtick',1:length(channels));
        set(gca,'xticklabel',channels);
        set(gca,'XTickLabelRotation',90);
        set(gca,'TickLabelInterpreter','none');
        %mesh(1:length(channels));
        
        
        %xticklabel_rotate([1:length(channels)],45,channels);
        colorbar
        title('Concatenated DREMI Heatmap of 18 FCS Files from CLL Samples Analyzed with PLAYR');
		
      
        cdata = cdata.write_mi_graph(channels, DREMI_cutoff);
       
        cdata.draw_denovo_DREMI_graph();
        sessionData.datasets{selected_gate} = cdata;
		
		
%       set(handles.lst_sigedges,'Visible','on');
%       set(handles.txt_selectedge,'Visible','on');
%       set(handles.txt_selectedge,'Visible','on');
%       set(handles.btn_sigedgeplot,'Visible','on');
%       num_sig_edges = size(cdata.DREMI_sig_edges,1);
%       sig_edges_string = cell(num_sig_edges,1);
%       for i=1:num_sig_edges
%           sig_edges_string{i} = sprintf('%s->%s',cdata.DREMI_sig_edges{i,1},cdata.DREMI_sig_edges{i,2});
%       end
%           
%       
%       set(handles.lst_sigedges,'String', sig_edges_string);
      
       put('sessionData', sessionData);
	   end
      
      %next I should draw the edges and and make them clickable 
end


function Etd_DREMICutoff_Callback(hObject, eventdata, handles)
% hObject    handle to Etd_DREMICutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Etd_DREMICutoff as text
%        str2double(get(hObject,'String')) returns contents of Etd_DREMICutoff as a double
end

% --- Executes during object creation, after setting all properties.
function Etd_DREMICutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Etd_DREMICutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in lst_sigedges.
function lst_sigedges_Callback(hObject, eventdata, handles)
% hObject    handle to lst_sigedges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lst_sigedges contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lst_sigedges
end

% --- Executes during object creation, after setting all properties.
function lst_sigedges_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lst_sigedges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in btn_sigedgeplot.
function btn_sigedgeplot_Callback(hObject, eventdata, handles)
% hObject    handle to btn_sigedgeplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   
      selected_gates = get(handles.lst_FIles, 'Value')
      
      if(length(selected_gates)==0)
          return;
      end
      
      sessionData = retr('sessionData');
      
       
      
      selected_sigedges = get(handles.lst_sigedges,'Value');
      threshold = str2num(get(handles.Edt_threshold,'String'));
      edge_list = get(handles.lst_sigedges,'String');
      nodes = {};
      mi_matrix = [];
      selected_edge_xy_list = cell(length(selected_sigedges),2); 
      all_channels = get(handles.popup_Yaxis,'String');
      
      for i = 1:length(selected_sigedges)
          
          %decipher what edge value this is for x and y 
          %how do I decipher what it is? 
          edge_string = edge_list{selected_sigedges(i)};
          parts = textscan(edge_string, '%s -> %s');
          j = find(strcmp(all_channels,parts{1}))
          channel1_name = all_channels{j}
          
          a = find(strcmp(nodes,channel1_name));
          if(length(a)==0)
              nodes{end+1} = channel1_name;
          end
          j = find(strcmp(all_channels,parts{2}))
         
          channel2_name = all_channels{j}
        
          a = find(strcmp(nodes,channel2_name));
          if(length(a)==0)
              nodes{end+1} = channel2_name;
          end
          
          selected_edge_xy_list{i,1} = channel1_name;
          selected_edge_xy_list{i,2} = channel2_name;
          
          [mi_vector, ~] = sessionData.compute_series_mi(channel1_name, channel2_name, selected_gates, threshold);
          mi_matrix = [mi_matrix; mi_vector];
      end
      
      
      CLIM = [0 1];
      cla
      subplot(1,1,1);
      load continuous_BuPu9.mat
      colormap(continuous_BuPu9)
      imagesc(mi_matrix);
      set(gca,'fontsize',14);
      set(gca,'YTick',1:length(selected_sigedges));
      ylabels = edge_list(selected_sigedges)
	  
	  
	  %debug
	  %disp(ylabels)
      set(gca, 'YTickLabel',ylabels);
      xticks = 1:length(selected_gates)
	  xlabels = edge_list(selected_sigedges)
      set(gca,'XTick',1:length(selected_gates));
      %timepoint_labels = sessionData.datalabels(selected_gates)
      %xticklabel_rotate([1:length(timepoint_labels)],45,timepoint_labels);
      title('DREMI Heatmap2');
      
     
      colorbar
      %additionally draw the biograph. 
      
      if(length(selected_gates)==1)
           selected_gate = selected_gates(1);
           cdata = sessionData.datasets{selected_gate};
           %need to call the appropriate plotting function  
    
           cdata.compute_and_graph_render_edges_dremi(selected_edge_xy_list, nodes);
      end
      
      %plot the heat map. 
end


% --- Executes when pnl_axes is resized.
function pnl_axes_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to pnl_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --- Executes on selection change in Popup_colorscheme.
function Popup_colorscheme_Callback(hObject, eventdata, handles)
% hObject    handle to Popup_colorscheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Popup_colorscheme contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Popup_colorscheme

    sessionData = retr('sessionData');
    selected_colorscheme = get(handles.Popup_colorscheme,'Value');
    num_datasets = length(sessionData.datasets);
    for i=1:num_datasets
        
        sessionData.datasets{i}.DREVI_colormode = selected_colorscheme;
        
    end
    
    put('sessionData', sessionData);
    
end

% --- Executes during object creation, after setting all properties.
function Popup_colorscheme_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Popup_colorscheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in btn_compute_AUC.
function btn_compute_AUC_Callback(hObject, eventdata, handles)
% hObject    handle to btn_compute_AUC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    all_channels = get(handles.popup_Yaxis,'String');
    selected_Ychannel = get(handles.popup_Yaxis,'Value');
    Ychannel = all_channels{selected_Ychannel};
    selected_Xchannel = get(handles.popup_Xaxis,'Value');
    Xchannel = all_channels{selected_Xchannel};
    

    selected_gates = get(handles.lst_FIles, 'Value');
    sessionData = retr('sessionData');
    
    [auc_values] = sessionData.compute_series_auc(Xchannel, Ychannel, selected_gates);
    
    for i=1:length(selected_gates)
       
        subplot(1,length(selected_gates),i);
        auc_string = sprintf('auc: %.2f',auc_values(i));
        title(auc_string);
    end


end


function Edt_maxy_Callback(hObject, eventdata, handles)
% hObject    handle to Edt_maxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edt_maxy as text
%        str2double(get(hObject,'String')) returns contents of Edt_maxy as a double
end

% --- Executes during object creation, after setting all properties.
function Edt_maxy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edt_maxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in Popup_autoY.
function Popup_autoY_Callback(hObject, eventdata, handles)
% hObject    handle to Popup_autoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Popup_autoY contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Popup_autoY


end


% --- Executes during object creation, after setting all properties.
function Popup_autoY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Popup_autoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ha = gca;
    ha.Parent
    map = colormap;

    f_new = figure;
    movegui(f_new,'northwest'); 
    copyobj(ha.Parent.Children, f_new);
    colormap(map);
    drawnow;
    
    % get new filename
    [filename, pathname, ~] = uiputfile({'*.png;*.eps;*.pdf', 'Image files (.png, .eps,.pdf)';...
                                         '*.png', 'Portable Network Graphics (.png)';...
                                         '*.eps', 'Encapsulated Postscript Vector Graphics (.eps)';...
                                         '*.pdf', 'Color PDF file format (.pdf)'},...
                                        'Save Image');
 
     if isequal(filename,0) || isequal(pathname,0)
        return;
     end
    
    % save new figure
    if (endswith(filename, '.pdf'))
        screen2pdf(f_new, [pathname filename]);
    elseif (endswith(filename, '.eps'))
        screen2eps(f_new, [pathname filename]);
    else 
        screen2png(f_new, [pathname filename]);
    end

end