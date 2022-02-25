filenames <- c(
  "DTIY00000000.fna",
  "JPGV00000000.fna",
  "NC_007775.fna",
  "NC_007776.fna",
  "NC_008536.fna",
  "NC_009523.fna",
  "NC_009767.fna",
  "NC_009828.fna",
  "NC_010175.fna",
  "NC_011831.fna",
  "NC_012032.fna",
  "NC_013946.fna",
  "NC_014212.fna",
  "NC_016024_NC_016025.fna",
  "NC_017464.fna",
  "NC_019701.fna",
  "NC_022795.fna",
  "NZ_ADVR00000000.fna",
  "NZ_AP014945.fna",
  "NZ_AP017560.fna",
  "NZ_ARKI00000000.fna",
  "NZ_AUIU00000000.fna",
  "NZ_AZXV00000000.fna",
  "NZ_CP007028.fna",
  "NZ_CP010822.fna",
  "NZ_CP016312.fna",
  "NZ_CP017675.fna",
  "NZ_CP021130.fna",
  "NZ_JPIM00000000.fna",
  "NZ_LWQS00000000.fna",
  "NZ_MAVV00000000.fna"
)

dir.in <- "data/amanda/ref/whole"

for (filename in filenames) {
  cat(filename,"\n")
  cat(names(readDNAStringSet(file.path(dir.in,filename))),sep = "\n")
  cat("\n")
}


