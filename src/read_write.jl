function write_chromosomes(input_prefix; output="chromosomes.txt")
    chrs = read_bim(string(input_prefix, ".bim")).CHR_CODE |> unique |> sort
    open(output, "w") do io
        for chr in chrs
            println(io, chr)
        end
    end
end