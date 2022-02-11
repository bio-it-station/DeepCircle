input_bed = commandArgs(trailingOnly = TRUE)[1]
genome_dir = commandArgs(trailingOnly = TRUE)[2]

get_1000_window = function(input_bed, genome_dir){
    
    out_dir = paste((dirname(getwd())), "/output/", sep = "")    
    bed_name = tail(unlist(strsplit(input_bed, "/")), 1)
    
    bed = read.csv(input_bed, sep = "\t", header = FALSE)

    number_records = nrow(bed)
    message(paste("Number of lines before processing:", number_records))

    #obtain 1000bp window centered on eccDNA midpoint
    smaller_1000 = bed[((bed[,3] - bed[,2]) <= 1000),]
    larger_1000 = bed[((bed[,3] - bed[,2]) > 1000),]
    window_1000 = bed[((bed[,3] - bed[,2]) <= 1000),]
    shift = (1000 - (window_1000[,3] - window_1000[,2]))/2
    window_1000[,2] = window_1000[,2] - floor(shift)
    window_1000[,3] = window_1000[,3] + ceiling(shift)
    message(paste("Number of records with 1000 bp window:", nrow(window_1000)))

    hg38 = read.csv(genome_dir, sep = "\t", header = FALSE)

    #adjust coordinate if < 0 or > chromosome size
    #if coordinate < 0 then change to 0
    window_1000[(window_1000[,2] < 0) ,2] = 0
    window_1000[(window_1000[,3] < 0) ,3] = 0

    #if coordinate > chromosome size then change to chromosome size
    size_of_matched_chr_vec = hg38[match(window_1000[,1],hg38[,1]),2]
    window_1000[window_1000[,2] > size_of_matched_chr_vec, 2] = size_of_matched_chr_vec[window_1000[,2] > size_of_matched_chr_vec]
    window_1000[window_1000[,3] > size_of_matched_chr_vec, 3] = size_of_matched_chr_vec[window_1000[,3] > size_of_matched_chr_vec]

    lower_than_0 = sum(window_1000[,2] < 0 | window_1000[,3] < 0)
    number_of_overbound = sum(window_1000[,2] > size_of_matched_chr_vec | window_1000[,3] > size_of_matched_chr_vec)
    message(paste("Number of records with exceeding boundary: ", number_of_overbound))
    message(paste("Number of records coordinate < 0:", lower_than_0))

    pre = unlist(strsplit(bed_name, "\\."))[1]
    out_path = paste(out_dir, pre, "_1000_window.bed", sep = "")
    smaller_out_path = paste(out_dir, pre, "_smaller_1000.bed", sep = "")
    larger_out_path = paste(out_dir, pre, "_larger_1000.bed", sep = "")

    write.table(window_1000, file = out_path, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(smaller_1000, file = smaller_out_path, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(larger_1000, file = larger_out_path, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

    message(paste("Window written to:", out_path))
    message(paste("Smaller than 1000 written to:", smaller_out_path))
    message(paste("Larger than 1000 written to:", larger_out_path))
}
get_1000_window(input_bed, genome_dir)
