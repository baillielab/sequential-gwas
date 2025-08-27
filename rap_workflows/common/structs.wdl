version 1.0

struct PGENFileset {
    String chr
    File pgen
    File psam
    File pvar
}

struct BGENFileset {
    String chr
    File bgen
    File bgi
    File sample
    File vcf_info
    File vcf_info_index
}

struct PLINKFileset {
    String chr
    File bed
    File bim
    File fam
}

struct BCFFileset {
    String chr
    File bcf
    File csi
}