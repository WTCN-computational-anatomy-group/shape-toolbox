function dtype = fileArray2FileArrayType(dtype, be)
% FORMAT mtype = fileArray2MatlabType(fatype)
% fatype - file_array data type namw (ex: 'FLOAT32-LE') or code (ex: 16)
% mtype  - Matlab datatype (ex: 'single')
    if nargin < 2
        be = 2;
    end
    
    if ischar(dtype)
        if be > 1
            if endsWith(dtype, '-le', 'IgnoreCase', true)
                be = 0;
            elseif endsWith(dtype, '-be', 'IgnoreCase', true)
                be = 1;
            else
                be = 2;
            end
        end
        if startsWith(dtype, 'UNKNOWN', 'IgnoreCase', true)
            dtype = 0;
        elseif startsWith(dtype, 'BINARY', 'IgnoreCase', true)
            dtype = 1;
        elseif startsWith(dtype, 'INT8', 'IgnoreCase', true)
            dtype = 256;
        elseif startsWith(dtype, 'UINT8', 'IgnoreCase', true)
            dtype = 2;
        elseif startsWith(dtype, 'INT16', 'IgnoreCase', true)
            dtype = 4;
        elseif startsWith(dtype, 'UINT16', 'IgnoreCase', true)
            dtype = 512;
        elseif startsWith(dtype, 'INT32', 'IgnoreCase', true)
            dtype = 8;
        elseif startsWith(dtype, 'UINT32', 'IgnoreCase', true)
            dtype = 768;
        elseif startsWith(dtype, 'INT64', 'IgnoreCase', true)
            dtype = 1024;
        elseif startsWith(dtype, 'UINT64', 'IgnoreCase', true)
            dtype = 1280;
        elseif startsWith(dtype, 'FLOAT32', 'IgnoreCase', true)
            dtype = 16;
        elseif startsWith(dtype, 'FLOAT64', 'IgnoreCase', true)
            dtype = 64;
        elseif startsWith(dtype, 'FLOAT128', 'IgnoreCase', true)
            dtype = 1536;
        elseif startsWith(dtype, 'COMPLEX64', 'IgnoreCase', true)
            dtype = 32;
        elseif startsWith(dtype, 'COMPLEX128', 'IgnoreCase', true)
            dtype = 1792;
        elseif startsWith(dtype, 'COMPLEX256', 'IgnoreCase', true)
            dtype = 2048;
        elseif startsWith(dtype, 'RGB24', 'IgnoreCase', true)
            dtype = 128;
        else
            dtype = 0;
        end
        dtype = double([dtype be]);
    else
        if length(dtype) == 2
            be = dtype(2);
            dtype = dtype(1);
        end
        switch dtype
            case 0
                dtype = 'UNKNOWN';
            case 1
                dtype = 'BINARY';
            case 256
                dtype = 'INT8';
            case 2
                dtype = 'UINT8';
            case 4
                dtype = 'INT16';
            case 512
                dtype = 'UINT16';
            case 8
                dtype = 'INT32';
            case 768
                dtype = 'UINT32';
            case 1024
                dtype = 'INT64';
            case 1280
                dtype = 'UINT64';
            case 16
                dtype = 'FLOAT32';
            case 64
                dtype = 'FLOAT64';
            case 1536
                dtype = 'FLOAT128';
            case 32
                dtype = 'COMPLEX64';
            case 1792
                dtype = 'COMPLEX128';
            case 2048
                dtype = 'COMPLEX256';
            case 128
                dtype = 'RGB24';
            otherwise
                dtype = 'UNKNOWN';
        end
        if be < 2
            if be
                dtype = [dtype '-BE'];
            else
                dtype = [dtype '-LE'];
            end
        end
    end
end