using PackageCompiler
using GenomiccWorkflows

sysimage_path = isempty(ARGS) ? joinpath(pkgdir(GenomiccWorkflows), "GenomiccWorkflows.so") : ARGS[1]
create_sysimage(["GenomiccWorkflows"]; sysimage_path=sysimage_path)