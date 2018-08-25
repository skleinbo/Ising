for package in ["GLM","Plots","LaTeXStrings", "StatPlots", "DataFrames", "DataFramesMeta","Interact","CSV"]
    try
        sp = Symbol(package)
        @eval import $sp
    catch err
        info("Package $package not installed. Trying to...")
        Pkg.add(package)
    end
end
