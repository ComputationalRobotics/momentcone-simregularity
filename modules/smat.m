%%*****************************************************************
%% smat: compute the matrix smat(x).
%%
%%   M = smat(blk,x,isspM); 
%%*****************************************************************
%% SDPT3: version 4.0
%% Copyright (c) 1997 by
%% Kim-Chuan Toh, Michael J. Todd, Reha H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************

   function M = smat(blk,xvec,isspM)

   if (nargin < 3); isspM = zeros(size(blk,1),1); end 
%%
   if ~iscell(xvec) 
      if strcmp(blk{1},'s')
         if exist('mexsmat', 'file') == 3
            M = mexsmat(blk,xvec,isspM);
         else
            M = local_smat_single_block(blk, xvec);
         end
      else
         M = xvec; 
      end   
   else 
      M = cell(size(blk,1),1);     
      if (length(isspM)==1)
         isspM = isspM*ones(size(blk,1),1); 
      end
      for p=1:size(blk,1)
         pblk = blk(p,:);
         if strcmp(pblk{1},'s');
            if exist('mexsmat', 'file') == 3
               M{p} = mexsmat(pblk,xvec{p},isspM(p));
            else
               M{p} = local_smat_single_block(pblk, xvec{p});
            end
         else
            M{p} = xvec{p}; 
         end   
      end
   end
%%*********************************************************
end

function M = local_smat_single_block(pblk, xvec)
   dims = pblk{2};
   if numel(dims) ~= 1
      error('local smat fallback supports a single SDP block only.');
   end
   n = dims(1);
   if numel(xvec) ~= n * (n + 1) / 2
      error('local smat fallback got inconsistent vector length.');
   end
   M = zeros(n, n);
   idx = 1;
   for j = 1:n
      for i = j:n
         v = xvec(idx);
         M(i, j) = v;
         M(j, i) = v;
         idx = idx + 1;
      end
   end
end

