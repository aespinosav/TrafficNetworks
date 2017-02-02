#Script to import .tntp files from Bar-Gera's test problems
using DataFrames

function remove_from_split_string(splitted_array, strs...)
    for s in strs
        for i in 1:length(splitted_array)
            splitted_array[i] = splitted_array[i][find(x -> x != s, splitted_array[i])]
        end
    end
    splitted_array
end

function load_tntp_to_dataframe(filename)

    open(filename, "r") do f
        lines = readlines(f)
    end

    splits = map(x -> split(x, '\t'), lines) #separate according to separator
    splits = remove_from_split_string(splits, "\n", "", ";\n", "~", "~ ") # get rid of unwanted characters
    splits = splits[find(x -> length(x) != 0, splits)] # get rid of empty lines

    zeroth_line = find(x -> x[1] == "<END OF METADATA>", splits)[1]
    header_line = zeroth_line + 1
    begining_line = zeroth_line + 2
    end_line = length(splits)

    n_edges = end_line - header_line

    data = DataFrame()

    for i in 1:length(splits[header_line])
        if i == 1 || 2
            typ = Int
        else
            typ = Float64
        end 
        data[Symbol(splits[header_line][i])] = @data [parse(typ, row[i]) for row in splits[begining_line:end] ]
    end
    data
end

function convert_to_tn(data_frame_from_tntp)
    
end