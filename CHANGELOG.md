# v1.5.5.dev
- Added code in `prepare_regions.sh` which outputs `metadata.txt` file which contains metadata about when and what version of `prepare_regions.sh` script was executed to generate the geneset data. 

# v1.5.4
Note that starting from this version, the version of the repository tag has diverged from container version. v1.5.4 is the uploaded container version, but the `burden_testing` code from v1.5.3 is used.

- Fixed `GMMAT` package version inside the container to  commit ID: `4063452` (See [GMMAT #34](https://github.com/hanchenphd/GMMAT/pull/34))
- Default container build file changed from Singularity to Dockerfile


# v1.5.3
- Relatedness matrices from GEMMA > 0.4 are now supported.
- Added `--method-optim` option for `step2`


# v1.5.2
1. Deleted singularity-deploy related files since remote building fails.
2. Changed singularity definition file.
    - Changed to compile HTSlib and BCFTools within /opt instead of /root, since this might have the chance of saving unwanted files in the host root's directory.
    - Changed git cloning of htslib and bcftools to downloading the recommended source tarball to compile the two libraries, since the container was missing the htslib (tabix, bgzip) and bcftools executables so prepare_regions.sh was failing.
3. Changes in prepare_regions.sh.
    - Refactored some code for improved readability.
    - Used absolute paths in some cd commands during VEP build to fix VEP build.
    - Fixed incorrect URL for downloading VEP's homo_sapiens_vep_98_GRCh38.tar.gz cache file.
    - Added code to download GENCODE data from http:// when download fails with ftp://
    - Changed fetching GENCODE data MD5 hash value with FTP to HTTP as fetching from FTP appears to hang in some environments (e.g. suanpan). Fetching with HTTP appeared more reliable across environments.
    - Added code to fetch APPRIS file's MD5 hash and compare with the downloaded APPRIS file's MD5 hash
    - Added Debug mode (-m), although not many debug messages added yet.

This has only been tested on our internal server, and not yet on other environments.

# v1.5.1
First tracked release.