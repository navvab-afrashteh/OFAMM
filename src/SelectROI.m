function handles = SelectROI(hObject, eventdata, handles)
[~, flag] = getimage(handles.MainAxes); % check if there is an image in the main axes
if flag
    handles.ROI.BW = roipoly;
    [r,c] = find(handles.ROI.BW);
    handles.ROI.xy = [c';r'];
    handles.ROI.selected = 1;
    guidata(gcbo,handles);
else
    handles.ROI.selected = 0;
end