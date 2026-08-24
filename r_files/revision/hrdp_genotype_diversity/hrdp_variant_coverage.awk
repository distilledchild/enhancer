BEGIN {
  FS = OFS = "\t"
  split(panel_names, panel_name, ",")
  split(study_names, study_name, ",")
}

function count_alleles(field, count, gt, allele, n, i, total) {
  gt = field
  sub(/:.*/, "", gt)
  gsub(/\|/, "/", gt)
  n = split(gt, allele, "/")
  total = 0
  for (i = 1; i <= n; i++) {
    if (allele[i] ~ /^[0-9]+$/) {
      count[allele[i]]++
      total++
    }
  }
  return total
}

/^#CHROM/ {
  for (i = 10; i <= NF; i++) sample_column[$i] = i

  for (i in panel_name) {
    if (!(panel_name[i] in sample_column)) {
      print "Missing HRDP specimen: " panel_name[i] > "/dev/stderr"
      exit 10
    }
    panel_column[i] = sample_column[panel_name[i]]
  }
  for (i in study_name) {
    if (!(study_name[i] in sample_column)) {
      print "Missing study specimen: " study_name[i] > "/dev/stderr"
      exit 11
    }
    study_column[i] = sample_column[study_name[i]]
  }
  next
}

/^#/ { next }

$7 == "." || $7 == "PASS" {
  quality = $6 + 0
  if (quality < 30) next

  delete panel_allele_count
  delete study_allele_count
  panel_total_alleles = 0
  for (i in panel_column) {
    panel_total_alleles += count_alleles($(panel_column[i]), panel_allele_count)
  }
  for (i in study_column) {
    count_alleles($(study_column[i]), study_allele_count)
  }

  n_alt = split($5, alt, ",")
  delete panel_site_type
  delete study_site_type

  for (j = 1; j <= n_alt; j++) {
    if (!(j in panel_allele_count) || alt[j] == "*" || alt[j] ~ /^</) next

    allele_frequency = panel_allele_count[j] / panel_total_alleles
    minor_allele_frequency = allele_frequency
    if (minor_allele_frequency > 0.5) {
      minor_allele_frequency = 1 - minor_allele_frequency
    }

    if (length($4) == 1 && length(alt[j]) == 1) {
      type = "SNP"
    } else if (length($4) != length(alt[j])) {
      type = "INDEL"
    } else {
      type = "OTHER"
    }

    for (m = 1; m <= 3; m++) {
      if (m == 2 && minor_allele_frequency <= 0.1) continue
      if (m == 3 && minor_allele_frequency <= 0.2) continue

      panel_allele[type, 30, m]++
      panel_site_type[type, m] = 1
      if (j in study_allele_count) {
        study_allele[type, 30, m]++
        study_site_type[type, m] = 1
      }

      if (quality >= 40) {
        panel_allele[type, 40, m]++
        if (j in study_allele_count) study_allele[type, 40, m]++
      }
    }
  }

  for (site_key in panel_site_type) {
    split(site_key, site_key_part, SUBSEP)
    type = site_key_part[1]
    m = site_key_part[2]
    panel_site[type, 30, m]++
    if ((type, m) in study_site_type) study_site[type, 30, m]++
    if (quality >= 40) {
      panel_site[type, 40, m]++
      if ((type, m) in study_site_type) study_site[type, 40, m]++
    }
  }
}

END {
  print "measurement_unit", "variant_type", "qual_threshold", \
    "maf_threshold", \
    "n_hrdp_panel_variants", "n_study_background_variants"

  type_name[1] = "SNP"
  type_name[2] = "INDEL"
  type_name[3] = "OTHER"
  threshold_value[1] = 30
  threshold_value[2] = 40
  maf_value[1] = 0
  maf_value[2] = 0.1
  maf_value[3] = 0.2
  unit_name[1] = "site"
  unit_name[2] = "allele"

  for (u = 1; u <= 2; u++) {
    unit = unit_name[u]
    for (i = 1; i <= 3; i++) {
      for (j = 1; j <= 2; j++) {
        for (m = 1; m <= 3; m++) {
          type = type_name[i]
          threshold = threshold_value[j]
          maf_threshold = maf_value[m]
          if (unit == "site") {
            print unit, type, threshold, maf_threshold, \
              panel_site[type, threshold, m] + 0, \
              study_site[type, threshold, m] + 0
          } else {
            print unit, type, threshold, maf_threshold, \
              panel_allele[type, threshold, m] + 0, \
              study_allele[type, threshold, m] + 0
          }
        }
      }
    }
  }
}
