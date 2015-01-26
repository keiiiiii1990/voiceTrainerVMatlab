function training_callbacks(varargin)

% TRAINING_CALLBACKS   Callback functions for the training interface
%
% training_callbacks(fname,varargin) executes the callback function
% fname_callback with various parameters

fname=[varargin{1},'_callback'];
feval(fname,varargin{2:end});



%%%exiting and possibly starting the evaluation phase
function proceed_callback(handles)

if handles.run_all,
    %suggesting a break
    wfig=warndlg(['You have been working for ' int2str(round(etime(clock,handles.time)/60)) 'minutes. It is recommended that you take a break of at least the same duration before starting the evaluation phase. Click on OK when you are ready.'],'Warning');
    uiwait(wfig);
    %exiting and starting evaluation with the same experiment order
    close(gcbf);
    %mushram('evaluation',handles.expe_order);
    evaluationVowelSpeech('evaluation');
else,
    %exiting
    close(gcbf);
end



%%%playing sound files

function play_callback(handles,e,f)

if f,
    %%{
    if ~isempty(handles.files{handles.expe_order(e),handles.file_order(e,f)})
        myplay(handles.files{handles.expe_order(e),handles.file_order(e,f)});
    end
    %}
    %myplay(handles.files{handles.expe_order(e),handles.file_order(e,f)});
else,
    %%{
    if ~isempty(handles.files{handles.expe_order(e),1})
        myplay(handles.files{handles.expe_order(e),1});
    end
    %}
    %myplay(handles.files{handles.expe_order(e),1});
end

function myplay(file)
if isunix,
    %using system's play (from the sox package) on Unix (MATLAB's sound does not work)
    %[s,w]=unix(['play ' file]);
    [s,w]=unix(['afplay ' file]);
else,
    %using MATLAB's wavplay on Windows
    [y,fs]=wavread(file);
    wavplay(y,fs);
end
