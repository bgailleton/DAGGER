module pork
  using CxxWrap
  @wrapmodule("./build/lib/libtestlib.so", :define_julia_module)

  function __init__()
    @initcxx
  end
end
