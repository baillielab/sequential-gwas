using PackageCompiler
using SequentialGWAS

sysimage_path = isempty(ARGS) ? joinpath(pkgdir(SequentialGWAS), "SequentialGWAS.so") : ARGS[1]
create_sysimage(["SequentialGWAS"]; sysimage_path=sysimage_path)