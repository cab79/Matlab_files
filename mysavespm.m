
%==========================================================================
function mysavespm(action,fname) % function from spm_results_ui
%==========================================================================
xSPM = evalin('caller','xSPM;');
XYZ  = xSPM.XYZ;

switch lower(action)
    case 'thresh'
        Z = xSPM.Z;
        
    case 'binary'
        Z = ones(size(xSPM.Z));
        
    case 'n-ary'
        if ~isfield(xSPM,'G')
            Z       = spm_clusters(XYZ);
            num     = max(Z);
            [n, ni] = sort(histc(Z,1:num), 2, 'descend');
            n       = size(ni);
            n(ni)   = 1:num;
            Z       = n(Z);
        else
            C       = NaN(1,size(xSPM.G.vertices,1));
            C(xSPM.XYZ(1,:)) = ones(size(xSPM.Z));
            C       = spm_mesh_clusters(xSPM.G,C);
            Z       = C(xSPM.XYZ(1,:));
        end
        
    case 'current'
        [xyzmm,i] = spm_XYZreg('NearestXYZ',...
            spm_results_ui('GetCoords'),xSPM.XYZmm);
        spm_results_ui('SetCoords',xSPM.XYZmm(:,i));
        
        if ~isfield(xSPM,'G')
            A   = spm_clusters(XYZ);
            j   = find(A == A(i));
            Z   = ones(1,numel(j));
            XYZ = xSPM.XYZ(:,j);
        else
            C   = NaN(1,size(xSPM.G.vertices,1));
            C(xSPM.XYZ(1,:)) = ones(size(xSPM.Z));
            C   = spm_mesh_clusters(xSPM.G,C);
            C   = C==C(xSPM.XYZ(1,i));
            Z   = C(xSPM.XYZ(1,:));
        end
        
    otherwise
        error('Unknown action.');
end

if isfield(xSPM,'G')
    F     = spm_input('Output filename',1,'s');
    if isempty(spm_file(F,'ext'))
        F = spm_file(F,'ext','.gii');
    end
    F     = spm_file(F,'CPath');
    M     = gifti(xSPM.G);
    C     = zeros(1,size(xSPM.G.vertices,1));
    C(xSPM.XYZ(1,:)) = Z; % or use NODE_INDEX
    M.cdata = C;
    save(M,F);
    cmd   = 'spm_mesh_render(''Disp'',''%s'')';
else
    V   = spm_write_filtered(Z, XYZ, xSPM.DIM, xSPM.M,...
    sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k),fname);
    cmd = 'spm_image(''display'',''%s'')';
    F   = V.fname;
end
fprintf('Written %s\n',spm_file(F,'link',cmd));                         %-#