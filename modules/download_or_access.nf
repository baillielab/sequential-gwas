process DownloadOrAccessStoredResource {
    storeDir "${store_dir}"

    input:
        tuple val(ftp_dir), val(filename)
        val store_dir

    output:
        path(filename)

    script:
        """
        wget -O ${filename} ${ftp_dir}/${filename}
        """
}
