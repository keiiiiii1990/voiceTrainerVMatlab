function customCursor(figureHandle)
cdata = ones(16);
cdata(3:14,3:14) = NaN;
cdata(1:2:16,8) = 1;
cdata(8,1:2:16) = 1;
cdata(2:2:16,9) = 1;
cdata(9,2:2:16) = 1;
set(figureHandle,'pointer','custom','PointerShapeCData',cdata,...
        'PointerShapeHotSpot',[9 9])
end