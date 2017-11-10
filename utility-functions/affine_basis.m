function [B, rind] = affine_basis(type, flat)
% FORMAT B = affine_basis(type, ('2d'))
% type - * 'translation'
%        * 'rotation'
%        * 'rigid'      or 6
%        * 'similitude' or 7
%        * 'affine'     or 12 [default]
% 
% B    - 4x4xQ array.
% rind - Indices of basis that shoudl be reularised (all but tr/rot)
%
% Returns a Lie algebra basis system encoding for one of the above
% transformation types.

    if nargin < 2
        flat = '3d';
        if nargin < 1
            type = 'affine';
        end
    end
    flat = ischar(flat) && strcmpi(flat, '2d');
    if ~ischar(type)
        type = num2str(type);
    end
    type = deblank(lower(type));
    
    % --- Define basis vectors
    
    Bt = [ 0 0 0 1   0 0 0 0   0 0 0 0 ;
           0 0 0 0   0 0 0 1   0 0 0 0 ;
           0 0 0 0   0 0 0 0   0 0 0 1 ;
           0 0 0 0   0 0 0 0   0 0 0 0 ];
     
    Br = [ 0 -1 0 0   0 0 -1 0   0 0  0 0 ;
           1  0 0 0   0 0  0 0   0 0 -1 0 ;
           0  0 0 0   1 0  0 0   0 1  0 0 ;
           0  0 0 0   0 0  0 0   0 0  0 0 ];
      
    Bsim = [ 1 0 0 0 ;
             0 1 0 0 ;
             0 0 1 0 ;
             0 0 0 0 ];
    
    Bscl = [ 1 0 0 0   0 0 0 0   0 0 0 0 ;
             0 0 0 0   0 1 0 0   0 0 0 0 ;
             0 0 0 0   0 0 0 0   0 0 1 0 ;
             0 0 0 0   0 0 0 0   0 0 0 0 ];
         
    Bshr = [ 0 1 0 0   0 0 1 0   0 0 0 0 ;
             1 0 0 0   0 0 0 0   0 0 1 0 ;
             0 0 0 0   1 0 0 0   0 1 0 0 ;
             0 0 0 0   0 0 0 0   0 0 0 0 ];
    
    % --- Remove 3D basis if 2D
    if flat
        Bt = Bt(:,1:8);
        Br = Br(:,1:4);
        Bsim = [ 1 0 0 0 ;
                 0 1 0 0 ;
                 0 0 0 0 ;
                 0 0 0 0 ];
        Bscl = Bscl(:,1:8);
        Bshr = Bshr(:,1:4);
    end
         
    % --- Build complete basis
    switch type
        case 'translation'
            B = Bt;
            rind = [];
        case 'rotation'
            B = Br;
            rind = [];
        case {'rigid', '6'}
            B = [Bt Br];
            rind = [];
        case {'similitude', '7'}
            B = [Bt Br Bsim];
            if flat,  rind = [4 5];
            else      rind = [7 8 9];  end
        case {'affine', '12'}
            B = [Bt Br Bscl Bshr];
            if flat,  rind = [4 5 6];
            else      rind = [7 8 9 10 11 12];  end
    end
    B = reshape(B, 4, 4, []);
    
end