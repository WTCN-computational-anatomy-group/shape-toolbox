function dtype = fileArray2MatlabType(dtype)
% FORMAT mtype = fileArray2MatlabType(fatype)
% fatype - file_array data type namw (ex: 'FLOAT32-LE') or code (ex: 16)
% mtype  - Matlab datatype (ex: 'single')
    if ischar(dtype)
        if startsWith(dtype, 'BINARY', 'IgnoreCase', true)
            dtype = 'logical';
        elseif startsWith(dtype, 'INT8', 'IgnoreCase', true)
            dtype = 'int8';
        elseif startsWith(dtype, 'UINT8', 'IgnoreCase', true)
            dtype = 'uint8';
        elseif startsWith(dtype, 'INT16', 'IgnoreCase', true)
            dtype = 'int16';
        elseif startsWith(dtype, 'UINT16', 'IgnoreCase', true)
            dtype = 'uint16';
        elseif startsWith(dtype, 'INT32', 'IgnoreCase', true)
            dtype = 'int32';
        elseif startsWith(dtype, 'UINT32', 'IgnoreCase', true)
            dtype = 'uint32';
        elseif startsWith(dtype, 'INT64', 'IgnoreCase', true)
            dtype = 'int64';
        elseif startsWith(dtype, 'UINT64', 'IgnoreCase', true)
            dtype = 'uint64';
        elseif startsWith(dtype, 'FLOAT32', 'IgnoreCase', true)
            dtype = 'single';
        elseif startsWith(dtype, 'FLOAT64', 'IgnoreCase', true)
            dtype = 'double';
        elseif startsWith(dtype, 'FLOAT128', 'IgnoreCase', true)
            dtype = 'double';
        elseif startsWith(dtype, 'COMPLEX64', 'IgnoreCase', true)
            dtype = 'single';
        elseif startsWith(dtype, 'COMPLEX128', 'IgnoreCase', true)
            dtype = 'double';
        elseif startsWith(dtype, 'COMPLEX256', 'IgnoreCase', true)
            dtype = 'double';
        else
            dtype = 'double';
        end
    else
        switch dtype
            case 1
                dtype = 'logical';
            case 256
                dtype = 'int8';
            case 2
                dtype = 'uint8';
            case 4
                dtype = 'int16';
            case 512
                dtype = 'uint16';
            case 8
                dtype = 'int32';
            case 768
                dtype = 'uint32';
            case 1024
                dtype = 'int64';
            case 1280
                dtype = 'uint64';
            case 16
                dtype = 'single';
            case 64
                dtype = 'double';
            case 1536
                dtype = 'double';
            case 32
                dtype = 'single';
            case 1792
                dtype = 'double';
            case 2048
                dtype = 'double';
            otherwise
                dtype = 'double';
        end
    end
end