function setLineWidth(handle,ws)
%handle=plot()
%set each line from handle=plot() width as numer in ws.
%When length(ws)<length(handle),repeat ws.
ws=reshape(ws,[1,length(ws)]);
if length(ws)<length(handle)
    ws=repmat(ws,[length(handle)/length(ws),1]);
end
for index=1:length(handle)
    handle(index).LineWidth=ws(index);
end
end