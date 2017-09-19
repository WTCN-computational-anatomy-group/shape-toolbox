function B = affine_basis(type, flat)
% FORMAT B = affine_basis(type, ('2d'))
% type - * 'translation'
%        * 'rotation'
%        * 'rigid'      or 6
%        * 'similitude' or 7
%        * 'affine'     or 12 [default]
% 
% B    - 4x4xQ array.
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
        case 'rotation'
            B = Br;
        case {'rigid', '6'}
            B = [Bt Br];
        case {'similitude', '7'}
            B = [Bt Br Bsim];
        case {'affine', '12'}
            B = [Bt Br Bscl Bshr];
    end
    B = reshape(B, 4, 4, []);
    
end