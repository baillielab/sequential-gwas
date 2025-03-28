using Pkg
Pkg.add("PackageCompiler")
using PackageCompiler

# create_app(
#     ".", 
#     "build";
#     executables=["seq-gwas" => "julia_main"],
# )

create_sysimage(["SequentialGWAS"], 
    sysimage_path="FlowOMMIC.so", 
    include_transitive_dependencies=false, 
    cpu_target=PackageCompiler.default_app_cpu_target()
)