function avgpath=mcxmeanpath_detID(id,detp)

 if(isfield(detp,'unitinmm'))
    unitinmm=detp.unitinmm;
else
    unitinmm=1;
 end

%Find idx
idx = find(detp.detid == id);

detw= mcxdetweight(detp,detp.prop,unitinmm);
avgpath=sum(detp.ppath(idx,:).*unitinmm.*repmat(detw(idx),1,size(detp.ppath,2))) / sum(detw(idx));

% detw=mcxdetweight_detID(idx,detp,detp.prop,unitinmm);
% avgpath=sum(detp.ppath(idx,:).*unitinmm.*repmat(detw(:),1,size(detp.ppath,2))) / sum(detw(:));
