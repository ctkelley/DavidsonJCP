function cry = finalize(cry)
% CRYSTAL/FINALIZE Finalize function for crystal class
%    cry = FINALIZE(cry) returns a crystal class of with finalized fields.
%
%    See also Crystal.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

cry = finalize@Molecule(cry);

if isempty(cry.kpts)
    cry.kpts = [0,0,0];
else
    % NOTE: the transpose is very important.
    Creci = 2*pi*inv(cry.supercell)';
    cry.kpts = cry.kpts*Creci;
end

cry.nkpts = size(cry.kpts,1);

if numel(cry.wks) ~= cry.nkpts
    cry.wks = ones(cry.nkpts,1)/cry.nkpts;
end

end
