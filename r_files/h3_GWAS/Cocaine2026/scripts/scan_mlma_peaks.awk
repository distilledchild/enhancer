BEGIN {
    FS = "[ \t]+"; OFS = "\t"
    number = "^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$"
}
NF == 0 { next }
{
    sub(/\r$/, "", $NF)
    if ($1 == "Chr") {
        headers++
        for (i = 1; i <= NF; i++) col[$i] = i
        if (!col["SNP"] || !col["p"] || !col["b"] || !col["se"] || !col["bp"]) {
            print "Missing required MLMA columns" > "/dev/stderr"; bad_header = 1; exit 2
        }
        next
    }
    if (!headers) { bad_header = 1; exit 2 }
    chr = $(col["Chr"]); seen[chr] = 1; rows[chr]++
    raw = $(col["p"])
    if (raw == "NA" || raw == "NaN" || raw == "nan" || raw == "-nan" || raw == ".") { missing[chr]++; next }
    if (raw !~ number || raw + 0 < 0 || raw + 0 > 1) { invalid[chr]++; next }
    p = raw + 0
    if (p == 0) { zero[chr]++; next }
    valid[chr]++
    if (p <= 1e-7) passing[chr]++
    if (!(chr in best) || p < best[chr]) {
        best[chr] = p; best_raw[chr] = raw; snp[chr] = $(col["SNP"])
        bp[chr] = $(col["bp"]); beta[chr] = $(col["b"]); se[chr] = $(col["se"])
        af[chr] = $(col["Freq"]); ties[chr] = 1
    } else if (p == best[chr]) ties[chr]++
}
END {
    if (bad_header) exit 2
    print "trait", "Chr", "n_rows", "n_valid", "n_missing_P", "n_invalid_P", "n_zero_P", "n_P_le_1e_minus7", "lead_SNP", "lead_bp", "min_P", "lead_beta", "lead_se", "lead_AF", "lead_ties", "headers_read"
    for (chr in seen) {
        print trait, chr, rows[chr]+0, valid[chr]+0, missing[chr]+0, invalid[chr]+0, zero[chr]+0, passing[chr]+0,
              (valid[chr] ? snp[chr] : "NA"), (valid[chr] ? bp[chr] : "NA"),
              (valid[chr] ? best_raw[chr] : "NA"), (valid[chr] ? beta[chr] : "NA"),
              (valid[chr] ? se[chr] : "NA"), (valid[chr] ? af[chr] : "NA"), ties[chr]+0, headers
    }
}
