module dagger
  using CxxWrap
  @wrapmodule("./libjudagger.so", :define_julia_module)

  function __init__()
    @initcxx
  end
end
